# Toy Model Simulation Script
#
# This script provides a simplified example of how to run a single simulation scenario
# using the functions defined in functions/simulation_functions.R.
# Its purpose is to illustrate the core mechanics of the simulation process
# in an accessible way for users new to the codebase.

# Load necessary libraries
# These libraries are commonly used in the main simulation scripts.
# Ensure they are installed in your R environment.
# install.packages(c("purrr", "dplyr", "tidyr", "ggplot2", "MASS", "parallel", "xgboost", "SuperLearner", "WeightIt", "chngpt", "cobalt", "CausalGPS"))
.libPaths(new = c("~/R/x86_64-pc-linux-gnu-library/4.2", .libPaths())) # Specific to FASSE cluster

library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
# library(parallel) # Not strictly needed for a single run, but often used
library(xgboost)
# library(SuperLearner) # May not be directly used in the simplest toy model path
library(WeightIt)
library(chngpt)
library(cobalt)
library(CausalGPS)

# --- Configuration ---
# Define the directory where the repository is located.
# Adjust this path to match your local setup.
# REPO_DIR <- "~/Desktop/Francesca_research/Simulation_studies/" # Example local path
# Or, use a relative path if running from within the repo's root directory
REPO_DIR <- "." # Assumes you are running this script from the root of the repository

# Source the simulation functions
# This file contains the core logic for data generation and analysis.
source(paste0(REPO_DIR, "/functions/simulation_functions.R"))

# --- Simulation Parameters ---
# These parameters define the specific scenario for our toy model.
# Refer to `functions/simulation_functions.R` and the main simulation scripts
# for more details on these parameters and their possible values.

# Set a seed for reproducibility
set.seed(123) # Allows you to get the same results every time you run the script

# Sample size for the simulated dataset
SAMPLE_SIZE <- 200 # A small sample size for quick execution

# GPS model specification (defines how treatment is generated based on covariates)
# 1: Linear relationship with normal error
# 2: Linear relationship with t-distributed error
# 3: Non-linear relationship (quadratic term) with normal error
# 4: Non-linear relationship (quadratic and interaction terms) with normal error
GPS_MODEL_SPEC <- 1 # Using the simplest linear model

# Exposure-response relationship (defines the true relationship between treatment and outcome)
# "linear": Outcome = Treatment + Confounders
# "sublinear": Outcome = log(Treatment) + Confounders
# "threshold": Outcome = I(Treatment > T_val) * (Treatment - T_val) + Confounders
EXPOSURE_RESPONSE <- "linear"

# Outcome interaction (defines if there's an interaction between treatment and confounders in the outcome model)
# TRUE: Outcome = Treatment + Confounders + Treatment * Confounders
# FALSE: Outcome = Treatment + Confounders
OUTCOME_INTERACTION_PRESENT <- FALSE # No interaction for simplicity

# Standard deviation for the outcome model's error term
OUTCOME_SD <- 10

# --- 1. Data Generation ---
# Generate synthetic data based on the specified parameters.
message(paste("Generating synthetic data with sample size:", SAMPLE_SIZE,
              "GPS model:", GPS_MODEL_SPEC,
              "Exposure-response:", EXPOSURE_RESPONSE,
              "Outcome interaction:", OUTCOME_INTERACTION_PRESENT))

sim_data <- sim_data_generate(
  sample_size = SAMPLE_SIZE,
  gps_mod = GPS_MODEL_SPEC,
  exposure_response_relationship = EXPOSURE_RESPONSE,
  outcome_interaction = OUTCOME_INTERACTION_PRESENT,
  outcome_sd = OUTCOME_SD,
  data_application = FALSE # We are generating synthetic data, not using real data
)

# Explore the generated data (optional)
message("Generated data summary:")
print(head(sim_data))
message(paste("Number of rows in generated data:", nrow(sim_data)))
message(paste("Number of columns in generated data:", ncol(sim_data)))

# --- 2. Model Fitting and Metric Calculation ---
# Fit various statistical models to the generated data and calculate performance metrics.
# This step mimics how different methods are evaluated in the main simulation.
message("Fitting models and calculating metrics...")

metrics_and_predictions <- metrics_from_data(
  sim_data = sim_data,
  exposure_response_relationship = EXPOSURE_RESPONSE,
  outcome_interaction = OUTCOME_INTERACTION_PRESENT
)

# The function 'metrics_from_data' returns a list containing:
# - metrics: Bias and MSE for each model at different exposure levels
# - predictions: Predicted outcome values for each model across a range of exposures
# - cor_table: Covariate balance (correlation) before and after weighting/matching
# - convergence_info: Information about the convergence of weighting methods

# Extract the results
model_metrics <- metrics_and_predictions$metrics
model_predictions <- metrics_and_predictions$predictions
correlation_table <- metrics_and_predictions$cor_table
convergence_status <- metrics_and_predictions$convergence_info

# --- 3. Reviewing Results ---
# Display some of the key results from the toy model run.

message("\n--- Toy Model Simulation Results ---")

message("\nConvergence Status of Weighting Methods:")
print(convergence_status)

message("\nCovariate Balance (Absolute Correlation with Exposure):")
# Show a snippet of the correlation table, e.g., for one method
print(head(correlation_table %>% filter(method == "ent"))) # Example for entropy balancing

message("\nModel Performance Metrics (Bias and MSE):")
# Show a snippet of the metrics, e.g., for one model at a specific exposure level
print(head(model_metrics %>% filter(model == "gam_model", exposure < 1))) # Example for GAM model

message("\nPredicted Exposure-Response Curves (First few points):")
print(head(model_predictions %>% dplyr::select(exposure, true_fit, gam_model, linear_model))) # Show true, GAM, and linear model predictions

# --- 4. Visualization (Optional) ---
# Create a simple plot of the estimated ERC for one or more models.
# This uses ggplot2, similar to the main analysis scripts.

# Prepare data for plotting (similar to what's done in metrics_from_data for its internal plot)
plot_data <- model_predictions %>%
  select(exposure, true_fit, gam_model, linear_model, ent_gam) %>% # Select a few models to plot
  pivot_longer(cols = -c(exposure, true_fit), names_to = "model_type", values_to = "predicted_outcome")

erc_plot <- ggplot(plot_data, aes(x = exposure)) +
  geom_line(aes(y = predicted_outcome, color = model_type, linetype = "Estimated ERC")) +
  geom_line(aes(y = true_fit, linetype = "True ERC"), color = "black") +
  labs(
    title = paste("Toy Model: Estimated vs. True Exposure-Response Curve"),
    subtitle = paste("Scenario: Linear ER, N=", SAMPLE_SIZE, ", GPS Mod=", GPS_MODEL_SPEC, sep=""),
    x = "Exposure Level",
    y = "Outcome",
    color = "Model Type",
    linetype = "Curve Type"
  ) +
  scale_linetype_manual(values = c("Estimated ERC" = "solid", "True ERC" = "dashed")) +
  theme_minimal() +
  theme(legend.position = "top")

# print(erc_plot) # Display the plot
# To save the plot:
# ggsave(paste0(REPO_DIR, "/figures/toy_model_erc_plot.png"), erc_plot, width = 8, height = 6)
# message(paste0("\nPlot saved to ", REPO_DIR, "/figures/toy_model_erc_plot.png (if figures directory exists)"))
# Ensure 'figures' directory exists or specify a different path.

# --- End of Toy Model Script ---
message("\nToy model simulation finished. Review the printed outputs and generated plot (if enabled).")
message("Modify parameters in the 'Simulation Parameters' section to explore other scenarios.")

# Note: This toy model runs a single iteration. The main simulation script `1_run_simulation.R`
# typically runs many iterations (e.g., 100 or 1000) across various parameter combinations
# and often utilizes parallel processing for speed, saving results to disk for later aggregation.
# This script is designed for understanding and experimentation.
