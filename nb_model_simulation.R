# Negative Binomial (NB) Model Simulation Script
#
# This script provides an example of how to run a single simulation scenario
# with a Negative Binomial (NB) distributed treatment.
# It utilizes the existing ZINB functionality (gps_mod = 5) by setting
# the zero-inflation probability (zinb_pi) to 0.
#
# Its purpose is to illustrate how to achieve NB treatment generation
# using the functions in functions/simulation_functions.R.

# Load necessary libraries
# Ensure they are installed in your R environment.
# install.packages(c("purrr", "dplyr", "tidyr", "ggplot2", "MASS", "parallel", "xgboost", "SuperLearner", "WeightIt", "chngpt", "cobalt", "CausalGPS"))
.libPaths(new = c("~/R/x86_64-pc-linux-gnu-library/4.2", .libPaths())) # Specific to FASSE cluster

library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
# library(parallel) # Not strictly needed for a single run
library(xgboost)
# library(SuperLearner)
library(WeightIt)
library(chngpt)
library(cobalt)
library(CausalGPS)

# --- Configuration ---
# Define the directory where the repository is located.
# Adjust this path to match your local setup.
# REPO_DIR <- "~/Desktop/Francesca_research/Simulation_studies/" # Example local path
REPO_DIR <- "." # Assumes running from the root of the repository

# Source the simulation functions
source(paste0(REPO_DIR, "/functions/simulation_functions.R"))

# --- Simulation Parameters for NB Treatment ---
set.seed(456) # Use a different seed for variety

SAMPLE_SIZE <- 300 # A different sample size for this example
GPS_MODEL_SPEC <- 5  # Use gps_mod = 5, which is the ZINB mechanism
EXPOSURE_RESPONSE <- "sublinear" # Example: sublinear relationship
OUTCOME_INTERACTION_PRESENT <- FALSE
OUTCOME_SD <- 10

# Parameters for Negative Binomial (achieved via ZINB settings)
NB_MEAN_PARAM <- 8       # This will be used for `zinb_mu` (base mean for NB component)
NB_DISPERSION_PARAM <- 1.5 # This will be used for `zinb_theta` (dispersion for NB component)
NB_ZERO_INFLATION_PROB <- 0 # CRITICAL: Set to 0 for pure Negative Binomial

message(paste("--- Negative Binomial (NB) Treatment Simulation ---"))
message(paste("Simulating NB by using ZINB mechanism (gps_mod=5) with zero-inflation probability = 0."))

# --- 1. Data Generation for NB Treatment ---
message(paste("\nGenerating synthetic data with NB treatment (Sample Size:", SAMPLE_SIZE,
              ", Base NB Mean:", NB_MEAN_PARAM, ", NB Dispersion:", NB_DISPERSION_PARAM, ")"))

sim_data_nb <- sim_data_generate(
  sample_size = SAMPLE_SIZE,
  gps_mod = GPS_MODEL_SPEC, # Should be 5
  exposure_response_relationship = EXPOSURE_RESPONSE,
  outcome_interaction = OUTCOME_INTERACTION_PRESENT,
  outcome_sd = OUTCOME_SD,
  data_application = FALSE,
  zinb_mu = NB_MEAN_PARAM,
  zinb_theta = NB_DISPERSION_PARAM,
  zinb_pi = NB_ZERO_INFLATION_PROB # Key for NB
)

# Explore the generated data (optional)
message("\nGenerated NB data summary:")
print(head(sim_data_nb))
message(paste("Number of rows in generated data:", nrow(sim_data_nb)))
summary(sim_data_nb$exposure)
message(paste("Proportion of actual zeros in exposure:", mean(sim_data_nb$exposure == 0)))
# For NB (pi=0), any zeros are from the NB distribution itself, not from inflation.
# If mu is high and theta not too small, actual zeros from NB might be rare.

# --- 2. Model Fitting and Metric Calculation ---
message("\nFitting models and calculating metrics for NB data...")

metrics_and_predictions_nb <- metrics_from_data(
  sim_data = sim_data_nb,
  exposure_response_relationship = EXPOSURE_RESPONSE,
  outcome_interaction = OUTCOME_INTERACTION_PRESENT
)

model_metrics_nb <- metrics_and_predictions_nb$metrics
model_predictions_nb <- metrics_and_predictions_nb$predictions
correlation_table_nb <- metrics_and_predictions_nb$cor_table
convergence_status_nb <- metrics_and_predictions_nb$convergence_info

# --- 3. Reviewing Results ---
message("\n--- NB Model Simulation Results ---")

message("\nConvergence Status of Weighting Methods (NB):")
print(convergence_status_nb)

message("\nCovariate Balance (Absolute Correlation with Exposure - NB):")
print(head(correlation_table_nb %>% filter(method == "ent")))

message("\nModel Performance Metrics (Bias and MSE - NB):")
print(head(model_metrics_nb %>% filter(model == "gam_model", exposure < 5))) # Example for GAM

message("\nPredicted Exposure-Response Curves (NB - First few points):")
print(head(model_predictions_nb %>% dplyr::select(exposure, true_fit, gam_model, linear_model)))

# --- 4. Visualization (Optional) ---
plot_data_nb <- model_predictions_nb %>%
  select(exposure, true_fit, gam_model, linear_model, ent_gam) %>%
  pivot_longer(cols = -c(exposure, true_fit), names_to = "model_type", values_to = "predicted_outcome")

erc_plot_nb <- ggplot(plot_data_nb, aes(x = exposure)) +
  geom_line(aes(y = predicted_outcome, color = model_type, linetype = "Estimated ERC")) +
  geom_line(aes(y = true_fit, linetype = "True ERC"), color = "black") +
  labs(
    title = "NB Model: Estimated vs. True Exposure-Response Curve",
    subtitle = paste("Scenario: ", EXPOSURE_RESPONSE, " ER, N=", SAMPLE_SIZE, ", Base NB Mean=", NB_MEAN_PARAM, sep=""),
    x = "Exposure Level",
    y = "Outcome",
    color = "Model Type",
    linetype = "Curve Type"
  ) +
  scale_linetype_manual(values = c("Estimated ERC" = "solid", "True ERC" = "dashed")) +
  theme_minimal() +
  theme(legend.position = "top")

# print(erc_plot_nb) # Display the plot
# To save:
# dir.create(paste0(REPO_DIR, "/figures"), showWarnings = FALSE) # Ensure directory exists
# ggsave(paste0(REPO_DIR, "/figures/nb_model_erc_plot.png"), erc_plot_nb, width = 8, height = 6)
# message(paste0("\nPlot saved to ", REPO_DIR, "/figures/nb_model_erc_plot.png (if figures directory exists)"))

message("\nNB model simulation finished.")
message("This script uses gps_mod=5 with zinb_pi=0 to simulate Negative Binomial treatment.")
message("Modify parameters as needed to explore other NB scenarios.")
