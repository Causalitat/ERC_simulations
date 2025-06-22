# Code used for simulation functions
library(dplyr)
library(purrr)
library(MASS)
# library(pscl) # For ZINB, if using pscl::rzinb. Or use custom function.
# library(countreg) # For ZINB using gamlss::rzinb

# Function to generate Zero-Inflated Negative Binomial (ZINB) random variables
# n: number of observations
# mu: mean of the negative binomial part
# theta: dispersion parameter of the negative binomial part (size)
# pi: zero-inflation probability
rzinb_custom <- function(n, mu, theta, pi) {
  # Generate from Negative Binomial
  nb_counts <- rnbinom(n, size = theta, mu = mu)
  # Generate zero-inflation component
  is_zero <- rbinom(n, size = 1, prob = pi)
  # Combine: if is_zero is 1, then count is 0, otherwise it's from NB
  counts <- ifelse(is_zero == 1, 0, nb_counts)
  return(counts)
}


# Function to generate synthetic data for simulations
# Sample size is sample for simulation 
# GPS mode refers to what assumptions are being made when generating exposure as a function of covariates
# Exposure response relationship refers to the exposure-outcome true relationship
# Outcome interaction refers to if there is an interaction between confounders and exposure in the outcome model
# Outcome sd is the sd of the outcome model
# Data application is whether to use covariates and exposure from real data or not
sim_data_generate <- function(sample_size = 1000, 
                              gps_mod = NA,
                              exposure_response_relationship = NA, 
                              outcome_interaction = FALSE,
                              outcome_sd = 10, 
                              data_application = FALSE,
                              data_app_sample = NULL,
                              zinb_mu = 5,      # Default mean for ZINB
                              zinb_theta = 1,   # Default dispersion for ZINB
                              zinb_pi = 0.3) {  # Default zero-inflation probability for ZINB
  
  # Make sure input parameters are correct
  if (!gps_mod %in% 1:5) { # Updated to include 5 for ZINB
    stop("Invalid value for gps_mod. It should be in 1:5.")
  }
  
  if (!exposure_response_relationship %in% c("linear", "sublinear", "threshold")) {
    stop("Invalid value for exposure_relationship. It should be 'linear', 'sublinear', or 'threshold'.")
  }
  
  # Define functions used for scaling exposure and gamma function
  cov_function <- function(confounders) as.vector(-0.8 + matrix(c(0.1, 0.1, -0.1, 0.2, 0.1, 0.1), nrow = 1) %*% t(confounders))
  scale_exposure <- function(x){20 * (x-min(x))/(max(x)-min(x))}
  
  if (!data_application) {
    # Generate the confounders
    cf <- mvrnorm(n = 2 * sample_size,
                  mu = rep(0, 4),
                  Sigma = diag(4))
    
    cf5 <- sample((-2):2, 2 * sample_size, replace = TRUE)
    cf6 <- runif(2 * sample_size, min = -3, max = 3)
    
    # Combine confounders and set names
    confounders_large <- data.frame(cf, cf5, cf6) %>%
      rename_with(~ paste0("cf", 1:6))
    
    # Generate appropriate exposure based on gps_mod
    # Note: Using 2*sample_size for initial generation, then filtering and sampling.
    
    common_linear_term <- cov_function(confounders_large) # This is applied to 2*sample_size rows

    if (gps_mod %in% 1:4) {
      exposure_values_raw <-
        case_when(
          # Ensure rhs of ~ matches length of confounders_large (i.e. 2*sample_size)
          gps_mod == 1 ~ 9 * common_linear_term + 18 + rnorm(2 * sample_size, 0, sqrt(10)),
          gps_mod == 2 ~ 9 * common_linear_term + 18 + sqrt(5) * rt(2 * sample_size, df = 3),
          gps_mod == 3 ~ 9 * common_linear_term + 15 + 2 * (confounders_large[, "cf3"])^2 + rnorm(2 * sample_size, 0, sqrt(10)),
          gps_mod == 4 ~ 9 * common_linear_term + 2 * confounders_large[, "cf3"]^2 + 2 * confounders_large[, "cf1"] * confounders_large[, "cf4"] + 15 + rnorm(2 * sample_size, 0, sqrt(10))
        )

      exposure_df <-
        tibble(confounders_large, exposure_raw = exposure_values_raw) %>%
        filter(exposure_raw > 0) %>%
        slice_sample(n = sample_size, replace = ifelse(nrow(.) < sample_size, TRUE, FALSE)) %>%
        rename(exposure = exposure_raw) %>%
        dplyr::select(exposure, starts_with("cf"))

    } else if (gps_mod == 5) { # ZINB model
      log_mu_base <- log(zinb_mu) # Use the zinb_mu parameter from function arguments

      # common_linear_term is already calculated based on confounders_large (2*sample_size)
      log_mu_nb_values <- log_mu_base + (common_linear_term / 2)
      mu_nb_values <- exp(log_mu_nb_values)

      mu_nb_values <- pmin(mu_nb_values, 50)
      mu_nb_values[mu_nb_values <= 0.01] <- 0.01

      exposure_values_raw <- rzinb_custom(n = 2 * sample_size, # Matches length of mu_nb_values
                                          mu = mu_nb_values,
                                          theta = zinb_theta, # zinb_theta from function arguments
                                          pi = zinb_pi)       # zinb_pi from function arguments

      exposure_df <-
        tibble(confounders_large, exposure_raw = exposure_values_raw) %>%
        slice_sample(n = sample_size, replace = FALSE) %>%
        rename(exposure = exposure_raw) %>%
        dplyr::select(exposure, starts_with("cf"))
    }

    # Now separate out exposure and confounders, name cf for brevity
    exposure <- exposure_df$exposure
    cf <- as.matrix(exposure_df %>% dplyr::select(-exposure))
  } else {
    # Use the data_app_sample argument
    exposure <- data_app_sample$exposure
    # Confounders are everything minus the exposure
    cf <- data_app_sample %>% dplyr::select(-exposure)
    cf <- as.matrix(cf)
  }
  
  # Add on exposure effect depending on relationship
  transformed_exp <- 
    case_when(
      exposure_response_relationship == "linear" ~ exposure,
      exposure_response_relationship == "sublinear" ~ 5 * log(exposure + 1),
      exposure_response_relationship == "threshold" ~ {
        1.5 * ifelse(exposure <= 5, 0, exposure - 5)
      }
    )
  
  # Add effect of the confounders and create outcome
  if (!outcome_interaction) {
    Y <- transformed_exp + as.numeric(20 - c(2, 2, 3, -1, -2, -2) %*% t(cf) + rnorm(sample_size, 0, outcome_sd))
  } else {
    Y <- transformed_exp + as.numeric(
      20 - c(2, 2, 3, -1, -2, -2) %*% t(cf) - 
        transformed_exp * (-0.1 * cf[, 1] + 0.1 * cf[, 3] ^ 2 + 0.1 * cf[, 4] + 0.1 * cf[, 5]) +
        rnorm(sample_size, 0, outcome_sd))
  }
  
  # Combine into one simulated dataset
  simulated_data <- tibble(Y, exposure, as_tibble(cf))
  
  return(simulated_data)
}




# This function fits all estimators to the data and extracts the absolute bias, RMSE
# Sim data is output from first function sim data generate
# Exposure response relationship determines if the relationship is linear, sublinear or threshold
# Outcome interaction is a logical that specifies if an interaction is present in the outcome model
metrics_from_data <- function(sim_data = sim_data,
                              exposure_response_relationship,
                              outcome_interaction) {
  
  # Rename for simplicity
  data_example <- sim_data
  
  # Make sure inputs are formatted correctly
  if (!exposure_response_relationship %in% c("linear", "sublinear", "threshold")) {
    stop("Invalid value for exposure_relationship. It should be 'linear', 'sublinear', or 'threshold'.")
  }
  
  # Now fit a linear model
  linear_fit <- glm(Y ~ exposure + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, 
                    data = data_example, 
                    family = "gaussian")
  
  # fit GAM model
  gam_fit <- mgcv::gam(Y ~ s(exposure, bs = 'cr', k = 4) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                       data = data_example, family = "gaussian")
  
  # Fit threshold model
  change_model <- chngptm(Y ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6, ~ exposure, 
                          family = "gaussian", type = "segmented", 
                          var.type = "bootstrap", save.boot = T, data = data_example)
  
  # Fit weighting models --------------------------------------
  max_attempts <- 5
  
  # Initialize a data frame to store convergence information
  convergence_info <- data.frame(
    method = c("ebal", "ebal_moments"),
    converged = NA
  )
  
  # List to store models
  weightit_models <- list(
    ebal = NULL,
    ebal_moments = NULL
  )
  
  # Define the model formula
  model_formula <- exposure ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6
  
  # List of methods
  methods <- list(
    ebal = list(method = "ebal", stabilize = TRUE),
    ebal_moments = list(method = "ebal", stabilize = TRUE, moments = 2)
  )
  
  # Loop over methods
  for (method_name in names(methods)) {
    
    converged <- FALSE  # Reset the convergence flag at the start of each method
    
    # Loop for max_attempts
    for (i in 1:max_attempts) {
      
      # Try to fit the model 
      weightit_models[[method_name]] <- do.call(weightit, c(list(formula = model_formula, data = data_example), methods[[method_name]]))
      
      # Check convergence (if weights are all very small this indicates it did not converge)
      if (all(weightit_models[[method_name]]$weights < 0.001) || any(weightit_models[[method_name]]$weights < 0)) {
        converged <- FALSE
      } else {
        converged <- TRUE
        break
      }
    }
    
    # Record the convergence status in the data frame
    convergence_info$converged[convergence_info$method == method_name] <- converged
    
    # If not converged after max_attempts, set to NA and log a message
    if (!converged) {
      warning(paste("Method", method_name, "did not converge after", max_attempts, "attempts"))
      # Convert to 1's for now, but will be discarded from simulation run
      weightit_models[[method_name]]$weights <- rep(1, nrow(data_example))
    }
  }
  
  
  # Add correct weights and truncate to 99.5th percentile
  data_example <- 
    data_example %>% 
    mutate(
      ent = weightit_models$ebal$weights,
      ent2 = weightit_models$ebal_moments$weights
    ) %>%
    mutate(
      ent = ifelse(ent > quantile(ent, 0.995), quantile(ent, 0.995), ent),
      ent2 = ifelse(ent2 > quantile(ent2, 0.995), quantile(ent2, 0.995), ent2)
    )
  
  correlation_table <- 
    purrr::map_dfr(c("ent", "ent2"), function(entropy_type) {
      # Grab weight of interest
      correct_weight <- data_example %>% dplyr::pull(!!rlang::sym(entropy_type))
      
      # Check if the weights are all zero or contain negative values
      if (all(correct_weight == 0) || any(correct_weight < 0)) {
        warning(paste("Skipping covariance calculation for", entropy_type, "due to non-positive weights"))
        return(data.frame(covariate = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"), 
                          pre_cor = NA, 
                          post_cor = NA, 
                          method = entropy_type))
      }
      
      post_cor <- cov.wt(data_example %>% dplyr::select(exposure, cf1:cf6), 
                         wt = correct_weight, cor = TRUE)$cor
      post_cor <- abs(post_cor[-1, 1])
      
      pre_cor <- cov.wt(data_example %>% dplyr::select(exposure, cf1:cf6), cor = TRUE)$cor
      pre_cor <- abs(pre_cor[-1, 1])
      
      tibble(covariate = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"), 
             pre_cor = pre_cor, 
             post_cor = post_cor, 
             method = entropy_type)
    })
  
  
  # Now fit the entropy weighting models
  entropy_lm <- glm(Y ~ exposure + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, 
                    data = data_example, weights = ent)
  entropy_gam <- mgcv::gam(Y ~ s(exposure, bs = 'cr', k = 4) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                           data = data_example, weights = ent)
  change_model_ent <- chngptm(Y ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6, ~ exposure, 
                              family = "gaussian", type = "segmented", var.type = "bootstrap", 
                              save.boot = T, data = data_example, weights = data_example$ent)
  
  # Now fit second moment entropy weighting
  entropy_lm2 <- glm(Y ~ exposure + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, 
                     data = data_example, weights = ent2)
  entropy_gam2 <- mgcv::gam(Y ~ s(exposure, bs = 'cr', k = 4) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, 
                            data = data_example, weights = ent2)
  change_model_ent2 <- chngptm(Y ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6, ~ exposure, 
                              family = "gaussian", type = "segmented", var.type = "bootstrap", 
                              save.boot = T, data = data_example, weights = data_example$ent2)
    
  # Now fit CausalGPS package ---------------------------------------------
  
  # Add ID column for CausalGPS package
  data_example <- mutate(data_example, id = row_number())
  
  # Set grid to tune over for causalGPS
  tune_grid <- tidyr::expand_grid(
    nrounds = 100,
    eta = c(0.05, 0.1, 0.2, 0.3),
    max_depth = c(2, 3),
    delta = seq(1, 3, by = 0.1)
  )
  
  # Make list to pass to run in parallel using mclapply
  tune_grid_list <- as.list(as.data.frame(t(tune_grid)))
  
  # Wrapper function for running causalGPS
  wrapper_func <- function(tune_param, return_type = "results") {
    # Run causalGPS with specified parameters
    suppressMessages({
      pseudo_pop_tune <- generate_pseudo_pop(Y = data.frame(dplyr::select(data_example, Y, id)),
                                             w = data.frame(dplyr::select(data_example, exposure, id)),
                                             c = data.frame(dplyr::select(data_example, cf1, cf2, cf3, cf4, cf5, cf6, id)),
                                             ci_appr = "matching",
                                             pred_model = "sl",
                                             use_cov_transform = FALSE,
                                             exposure_trim_qtls = c(0.01, 0.99),
                                             optimized_compile = T,
                                             sl_lib = c("m_xgboost"),
                                             params = list(xgb_rounds = tune_param[[1]],
                                                           xgb_eta = tune_param[[2]],
                                                           xgb_max_depth = tune_param[[3]]),
                                             covar_bl_method = "absolute",
                                             covar_bl_trs = 0.1,
                                             covar_bl_trs_type = "mean",
                                             dist_measure = "l1",
                                             max_attempt = 1,
                                             matching_fun = "matching_l1",
                                             delta_n = tune_param[[4]],
                                             scale = 1,
                                             nthread = 1)
    })
    
    matched_pop_tune <- pseudo_pop_tune$pseudo_pop
    
    # Now generate covariate balance
    post_cor <- cov.wt(matched_pop_tune %>% dplyr::select(exposure, cf1:cf6), 
                       wt = matched_pop_tune$counter_weight, cor = TRUE)$cor
    post_cor <- abs(post_cor[-1, 1])
    mean_post_cor <- mean(post_cor)
    
    # Return metrics for results
    results <- 
      data.frame(nrounds = tune_param[[1]], 
                 eta = tune_param[[2]],
                 max_depth = tune_param[[3]],
                 delta = tune_param[[4]],
                 scale = 1,
                 mean_corr = mean_post_cor,
                 max_corr = max(post_cor))
    
    if (return_type == "results") {
      return(results)
    } else if (return_type == "matched_pop") {
      return(matched_pop_tune)
    } else {
      stop("Invalid return_type specified.")
    }
  }
  
  # Now iterate through simulation function, bind together results
  # For now using lapply to reduce memory load
  pseudo_pop_list <- mclapply(tune_grid_list, mc.cores = 10, wrapper_func)
  corr_search <- do.call("rbind", pseudo_pop_list)
  min_corr = corr_search %>% filter(mean_corr == min(mean_corr)) %>% slice_head()

  # Extract minimum as your result
  pseudo_pop_tuned <- wrapper_func(as.list(min_corr[1:4]), return_type = "matched_pop")
  
  # Check that correlation
  post_cor <- cov.wt(pseudo_pop_tuned[, c("exposure", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6")], 
                     wt = (pseudo_pop_tuned$counter_weight), cor = T)$cor
  post_cor <- abs(post_cor[-1, 1])
  
  # Add to correlation table
  cor_causalgps <- 
    data.frame(covariate = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"),
               pre_cor = correlation_table %>% filter(method == "ent") %>% arrange(covariate) %>% pull(pre_cor),
               post_cor = post_cor,
               method = "causal_gps_tuned",
               delta = min_corr[["delta"]])
  
  # Delta is not relevant for other models but add to correlation table
  correlation_table$delta = NA
  correlation_table <- rbind(correlation_table, cor_causalgps)
  
  # Now outcome model with weights
  causal_gps_tuned <- 
    mgcv::gam(formula = Y ~ s(exposure, bs = 'cr', k = 4) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
              family = "gaussian",
              data = data.frame(pseudo_pop_tuned),
              weights = counter_weight)
  
  # Now many points to evaluate the ERC 
  discrete_points = 100
  
  # Predict the ERC for every point
  data_prediction <- 
    map_dfr(seq(0, 20, length.out = discrete_points), function(pot_exp) {
      
      # Create potential outcome population at exposure
      # Create potential outcome population at exposure
      potential_data <- 
        data_example %>% 
        dplyr::select(cf1:cf6) %>%
        mutate(exposure = pot_exp,
               transformed_exp = case_when(
                 exposure_response_relationship == "linear" ~ exposure,
                 exposure_response_relationship == "sublinear" ~ 5 * log(exposure + 1),
                 exposure_response_relationship == "threshold" ~ 1.5 * ifelse(exposure <= 5, 0, exposure - 5),
                 TRUE ~ NA_real_
               ))
      
      confounder_matrix <- as.matrix(dplyr::select(potential_data, cf1:cf6))
      
      # Add true fit 
      Y <- 
        potential_data$transformed_exp + 
        if (outcome_interaction) {
          as.numeric(20 - c(2, 2, 3, -1, -2, -2) %*% t(confounder_matrix)) - 
            potential_data$transformed_exp * (-0.1 * potential_data$cf1 + 0.1 * potential_data$cf3^2 + 0.1 * potential_data$cf4 + 0.1 * potential_data$cf5)
        } else {
          as.numeric(20 - c(2, 2, 3, -1, -2, -2) %*% t(confounder_matrix))
        }
      
      # Fit each model and take mean for ERC
      potential_outcome <- 
        potential_data %>% 
        mutate(linear_model = predict(linear_fit, potential_data, type = "response"),
               gam_model = predict(gam_fit, newdata = potential_data, type = "response"),
               change_model = predict(change_model, newdata = potential_data, type = "response"),
               ent_linear = predict(entropy_lm, newdata = potential_data, type = "response"),
               ent_gam = predict(entropy_gam, newdata = potential_data, type = "response"),
               ent_change = predict(change_model_ent, newdata = potential_data, type = "response"),
               ent_linear2 = predict(entropy_lm2, newdata = potential_data, type = "response"), # Adjust for second moment
               ent_gam2 = predict(entropy_gam2, newdata = potential_data, type = "response"),
               ent_change2 = predict(change_model_ent2, newdata = potential_data, type = "response"),
               causal_gps_tuned = predict(causal_gps_tuned, newdata = potential_data, type = "response"),
               true_fit = Y) %>% 
        dplyr::select(exposure, linear_model, gam_model, change_model, ent_linear, ent_gam, ent_change,
                      ent_linear2, ent_gam2, ent_change2, causal_gps_tuned, true_fit) %>% 
        summarize_all(mean) 
      return(potential_outcome)
    })
  
  model_types <- c("linear_model", "gam_model", "change_model", "ent_linear", "ent_gam", 
                   "ent_change","ent_linear2", "ent_gam2", "ent_change2","causal_gps_tuned")
  
  # Separate true_fit from the prediction data
  true_fit_data <- data_prediction %>% dplyr::select(exposure, true_fit)
  
  # Diagnostic plot of model fit
  plot_fits <- 
    data_prediction %>% 
    dplyr::select(-true_fit) %>%
    pivot_longer(all_of(model_types), names_to = "model", values_to = "prediction") %>%
    mutate(model = factor(model, levels = model_types)) %>% # Convert model to a factor with specified levels order
    left_join(true_fit_data, by = "exposure") %>%
    arrange(exposure) %>%
    ggplot(aes(x = exposure, y = prediction, color = model, linetype = "Prediction")) +
    geom_line() +
    geom_line(aes(y = true_fit, linetype = "True Fit"), data = true_fit_data, color = "black", linetype = "dashed") +
    labs(x = "Exposure concentration", y = "Relative risk of death") +
    facet_wrap(~ model) +
    theme_bw() +
    scale_linetype_manual(values = c("Prediction" = "solid", "True Fit" = "dashed"))
  
  # Save Bias and MSE by model type and exposure value
  data_metrics <- 
    data_prediction %>% 
    tidyr::pivot_longer(cols = all_of(model_types), names_to = "model", values_to = "prediction") %>%
    dplyr::group_by(model, exposure) %>%
    dplyr::summarize(
      bias = mean(prediction - true_fit, na.rm = TRUE),
      mse = mean((prediction - true_fit)^2, na.rm = TRUE),
      .groups = 'drop' # this will ungroup the data after summarizing
    )
  
  return_columns <- c("exposure", model_types, "true_fit")
  return(list(metrics = data_metrics, 
              predictions = dplyr::select(data_prediction, !!return_columns), 
              cor_table = correlation_table,
              convergence_info = convergence_info))
}

