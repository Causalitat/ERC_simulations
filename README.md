# ERC simulations
Code for running analysis in: Methods for Estimating the Exposure-Response Curve to Inform the New Safety Standards for Fine Particulate Matter

This respository contains code used to conduct simulation studies comparing methods for exposure-response curves. It also contains a data application to the Medicare data set. 

The main code scripts are as follows:
* `0_launch.jobs.sh` is script to launch jobs on the FASSE cluster
* `1_run_simulation.R` is script to launch simulation for one exposure-outcome model, which launches 100 jobs
* `2_aggregate_simulation.R` is script bind together 100 simulation results to calculate predictions and metrics
* `3_paper_figures.R` is script to produce figures used in the paper
* `/functions/` contains scripts for running simulation
* `/figures/` saves figures used in publication
* `/markdown_files/` contains scripts used for exploratory data analysis and making additional plots

## Toy Model for Understanding Simulations
A new script `toy_model_simulation.R` has been added to help users understand the simulation process.
* **Purpose**: This script runs a single, simplified simulation scenario. It is heavily commented to explain each step, from data generation to model fitting and results evaluation. It's an excellent starting point for new users to grasp the core mechanics before diving into the more complex, large-scale simulations.
* **How to Run**:
    1. Open `toy_model_simulation.R` in RStudio or your preferred R environment.
    2. Modify the `REPO_DIR` variable at the top of the script if your working directory is not the root of the repository.
    3. Adjust simulation parameters (e.g., `SAMPLE_SIZE`, `GPS_MODEL_SPEC`, `EXPOSURE_RESPONSE`) as desired to explore different scenarios.
    4. Run the script. It will print output to the console and can be configured to save a plot of the estimated exposure-response curve.
* **Key Features**:
    * Uses the same core functions (`sim_data_generate`, `metrics_from_data`) as the main simulation.
    * Focuses on clarity and step-by-step execution.
    * Allows easy experimentation with different simulation settings.

## Simulation Details

### Treatment Generation (`gps_mod`)
The `gps_mod` parameter in `sim_data_generate` (within `functions/simulation_functions.R`) controls how the exposure (treatment) is generated based on covariates. The following options are available:

*   `gps_mod = 1`: Linear relationship with normal error. Exposure ~ β0 + βX + N(0, σ^2).
*   `gps_mod = 2`: Linear relationship with t-distributed error. Exposure ~ β0 + βX + t(df).
*   `gps_mod = 3`: Non-linear relationship (quadratic term) with normal error. Exposure ~ β0 + βX + β_k*X_k^2 + N(0, σ^2).
*   `gps_mod = 4`: Non-linear relationship (quadratic and interaction terms) with normal error. Exposure ~ β0 + βX + β_k*X_k^2 + β_ij*X_i*X_j + N(0, σ^2).
*   `gps_mod = 5`: Zero-Inflated Negative Binomial (ZINB) distribution.
    *   This option generates count data with excess zeros, which can be useful for modeling treatments that are often zero but positive otherwise.
    *   The mean of the Negative Binomial component (μ_NB) is modeled as a function of covariates: `μ_NB = exp(log(zinb_mu_base) + common_linear_term / 2)`.
    *   The ZINB generation uses the following parameters (which can be passed to `sim_data_generate`):
        *   `zinb_mu` (Default: 5): The base mean for the Negative Binomial component. The actual mean for each observation will vary based on its covariates.
        *   `zinb_theta` (Default: 1): The dispersion parameter (size) for the Negative Binomial component. Higher values mean less dispersion.
        *   `zinb_pi` (Default: 0.3): The zero-inflation probability. This is the probability that an observation is an "excess" zero, regardless of the NB component.

### Outcome Generation
The outcome `Y` is generated based on the exposure, confounders, and specified error distribution. Key parameters include:
* `exposure_response_relationship`: Defines the true functional form between exposure and outcome (e.g., "linear", "sublinear", "threshold").
* `outcome_interaction`: A boolean (`TRUE`/`FALSE`) indicating whether there's an interaction term between the (transformed) exposure and some confounders in the true outcome model.

## Data application
Code related to fitting model to Medicare database is found in `/data_application/`

* `0_data_prep.R` is script to load data and prepare it for analysis. This is run each time a new input dataset is prepared. This also prepares the bootstrap samples
* `0_tune_CausalGPS.R` is script that tunes the XGBoost and caliper hyperparameters for the CausalGPS package, this only needs to be run once for each input dataset to find optimal parameters
* `1_model_fit.R` is script to fit all ERC estimators to the data and generate the ERC for the Medicare dataset
* `1_model_fit_boot.R` is script to fit all ERC estimators to the bootstrap data and generate the ERC for the Medicare dataset
* `2_post_processing.R` is script to generate plots of ERC as well as covariate balance plots with included uncertainty from bootstrap
* `/functions/` contains two functions used for data application:
  * `/functions/entropy_wt_functions.R` is a more efficient implementation of entropy weighting for a large dataset
  * `/functions/fit_causalGPS_by_year.R` fits the CausalGPS package by year for Medicare dataset
* `bash_*.sh` all these scripts are used to launch jobs on the FASSE cluster 
