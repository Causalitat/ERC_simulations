# Function for running ERC sim on FASSE cluster 

# Change library path for running on cluster with updates
.libPaths(new = c("~/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))

# Load libraries needed 
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(parallel)
library(xgboost)
library(SuperLearner)
library(WeightIt)
library(chngpt)
library(cobalt)
library(CausalGPS)

# File paths
# Set directory (based on username on cluster or local computer)
if (Sys.getenv("USER") == "mcork") {
  repo_dir <- "~/nsaph_projects/ERC_simulation/Simulation_studies/"
} else if (Sys.getenv("USER") == "michaelcork") {
  repo_dir <- "~/Desktop/Francesca_research/Simulation_studies/"
}

# Define simulation number
sim.num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# pass argument for model run
args <- commandArgs(T)
exposure_response_relationship = as.character(args[[1]])
run_title = as.character(args[[2]])

# Create out directory
out_dir <- paste0(repo_dir, "results/", exposure_response_relationship, "_TRUE")
# if (!dir.exists(out_dir)) dir.create(out_dir)

# Create correct output directory (previously using time, now using whatever character you want)
dir.create(paste0(out_dir, "/", run_title))
out_dir <- paste0(out_dir, "/", run_title)

# Source the appropriate functions needed for simulation 
source(paste0(repo_dir, "/functions/simulation_functions.R"))

# Start the simulations ------------------------------------------------------------------------------------------
set.seed(sim.num)

# Store output from for loop
sim_metrics <- tibble()
sim_predictions <- tibble()
sim_cor_table <- tibble()
sim_convergence <- tibble()


# loop through all sample sizes included
for (sample_size in c(200, 1000, 10000)) {
  message(paste("Running with sample size", sample_size))
  # loop through all gps model specifications (now includes 5 for ZINB)
  for (gps_mod in 1:5) {
    message(paste("Running with gps model", gps_mod))

    # Specific parameters for ZINB (gps_mod == 5).
    # These are currently set to the defaults in sim_data_generate.
    # They can be customized here if needed for specific simulation runs.
    # e.g., zinb_params <- list(zinb_mu = 10, zinb_theta = 2, zinb_pi = 0.25)
    # And then pass these to sim_data_generate:
    # sim_data_generate(..., zinb_mu = zinb_params$zinb_mu, ...)
    # For now, we rely on the defaults in sim_data_generate.

    # Loop through two different outcome relationships
    for (outcome_interaction in c("T", "F")) {
      
      if (outcome_interaction) {
        message("Running with interaction in outcome model")
      } else {
        message("Running without interaction in outcome model")
      }
      
      # Generate synthetic data for simulation study
      sim_data <- sim_data_generate(sample_size = sample_size, 
                                    gps_mod = gps_mod,
                                    exposure_response_relationship = exposure_response_relationship, 
                                    outcome_interaction = as.logical(outcome_interaction),
                                    outcome_sd = 10, 
                                    data_application = FALSE)
        
        # Get metrics and predictions from sample
        metrics_predictions <- 
          metrics_from_data(sim_data = sim_data,
                            exposure_response_relationship = exposure_response_relationship, 
                            outcome_interaction = as.logical(outcome_interaction))
        
        metrics <- metrics_predictions$metrics
        predictions <-  metrics_predictions$predictions
        cor_table <- metrics_predictions$cor_table
        convergence_info <- metrics_predictions$convergence_info
        
        # Add appropriate columns to metrics and predictions
        metrics <- 
          metrics %>% 
          mutate(gps_mod = gps_mod,
                 sample_size = sample_size,
                 sim = sim.num,
                 outcome_interaction = outcome_interaction)
        
        predictions <- 
          predictions %>% 
          mutate(gps_mod = gps_mod, 
                 sample_size = sample_size, 
                 sim = sim.num,
                 outcome_interaction = outcome_interaction)
        
        cor_table <- 
          cor_table %>% 
          mutate(gps_mod = gps_mod, 
                 sample_size = sample_size, 
                 sim = sim.num,
                 outcome_interaction = outcome_interaction)
        
        # Add on convergence info for entropy weighting
        convergence_info <- 
          convergence_info %>% 
          mutate(gps_mod = gps_mod, 
                 sample_size = sample_size, 
                 sim = sim.num,
                 outcome_interaction = outcome_interaction)
        
        
        # Bind together data table with all results
        sim_metrics <- bind_rows(sim_metrics, metrics)
        sim_predictions <- bind_rows(sim_predictions, predictions)
        sim_cor_table <- bind_rows(sim_cor_table, cor_table)
        sim_convergence <- bind_rows(sim_convergence, convergence_info)
    }
  }
}

message("Saving output")

# Save the final sim results in list and write to output directory
sim_final <- list(sim_metrics, sim_predictions, sim_cor_table, sim_convergence)
saveRDS(sim_final, file = paste0(out_dir, "/sim_final_", sim.num, ".RDS"))
