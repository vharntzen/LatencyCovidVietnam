
# Run small simulation on the server to check if the way we address right truncation
# in JAGS for the doubly interval censored data is correct.
# The script uses paralellization.
#
# The data generation code can be found in the help file of the estimation 
# function in our R package 'doublIn' as well.

date_of_today <- "230324"
n_runs <- 500

# Load software
source("fun_estimate_dic.R")
require(dplyr)
require(parallel)
require(patchwork)

# Pick up the number of cores; HASH WHEN RUNNING LOCALLY
cores <- system("nproc", intern = TRUE)

# Pick up the setting for the number of test moments (1-5)
# HASH WHEN RUNNING LOCALLY
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
scenario <- as.numeric(taskid)

# The simulation function
Simulate_truncation <- function(i, scenario = scenario){

N <- 1000
  
# 1. Generate the data
  
  # Draw an exposure window width 1, 2, 3, 4, 5
  L0 <- rep(0, N)
  L1 <- sample(1:5, N, replace = T)
  
  # Draw the infection moment from a uniform distribution on (L0, L1)
  L <- runif(N, 0, L1)
  
  # Draw latency times
  # As estimated by Xin et al., 2022 (doi: 10.1093/cid/ciab746)
  times <- rgamma(N, shape = 4.05, rate = 0.74)
  R <- L + times 
  
  # Create empty vectors
  R0 <- R1 <- Trunc <- rep(NA, N)
  
  # Define the data set
  mydat <- data.frame(L, L0, L1, R , R0, R1, Trunc)
  
  # (i) draw test moments
  if(scenario == 1){last_test <- sample( c(5, 10, 15, 20, 25), N, replace = T)}
  if(scenario == 2){last_test <- sample( c(14, 21), N, replace = T)}
  if(scenario == 3){last_test <- rep(7, N)}
  
  # For each individual
  for(r in 1:nrow(mydat)){

    
    # (ii) define a window containing the endpoint
    mydat$R1[r] <- floor(mydat$R[r]) + 1
    mydat$R0[r] <- floor(mydat$R[r])
    
    # (iii) apply truncation
    if( mydat$R[r] > (mydat$L1[r] + last_test[r]) ){
      
      mydat$Trunc[r] <- NA
      
    } else {
      
      # (iv) set the truncation time equal to the last test moment
      mydat$Trunc[r] <-  mydat$L1[r] + last_test[r]
      
    }
    
  }
  
  # We set end of exposure back to first positive test day but in this
  # exercise we don't want to introduce any dependency.
  mydat$L1 <- ifelse(mydat$L1 > mydat$R1, mydat$R1, mydat$L1)
  
# 2. Truncate the data
  truncated_data <- mydat %>% filter(!is.na(Trunc))
  
# 3. Fit model addressing truncation
addressing_truncation <- 
    Estimate_doubly_interval_censored(truncated_data, 
                              infection_risk_distribution = "constant", 
                              method = "gamma", 
                              percentiles = c(0.5, 0.9, 0.95, 0.99),
                              right_truncation = T, iters = 5500, 
                              burnin_period = 500, thin = 1, 
                              further_thin_plots = 1)
  
# 4. Fit model NOT addressing truncation
no_addressing_truncation <- Estimate_doubly_interval_censored(
    truncated_data,
    infection_risk_distribution = "constant",
    method = "gamma",
    percentiles = c(0.5, 0.9, 0.95, 0.99),
    right_truncation = F, iters = 5500,
    burnin_period = 500, thin = 1,
    further_thin_plots = 1)
  
# 5. Fit the 'truth'

# Process output - truncation addressed
out_addressed <- as.data.frame(addressing_truncation$estimates[c(6,8), 1:3])
out_addressed$truncation <- "addressed"
out_addressed$truth_dat <- quantile(mydat$R - mydat$L, c(0.5, 0.95) )[1:2]
out_addressed$truth <- qgamma(shape = 4.05, rate = 0.74, p = c(0.5, 0.95))
out_addressed$deviation <- out_addressed$est - out_addressed$truth
out_addressed$true_in_CI <- ifelse(out_addressed$truth >= out_addressed$lower_CI
                                & out_addressed$truth <= out_addressed$upper_CI,
                                   1, 0)
out_addressed$quantile <- c(0.5, 0.95)

# Process output - truncation unaddressed
out_unaddressed <- as.data.frame(no_addressing_truncation$estimates[c(6,8), 1:3])
out_unaddressed$truncation <- "unaddressed"
out_unaddressed$truth_dat <- quantile(mydat$R - mydat$L, c(0.5, 0.95) )[1:2]
out_unaddressed$truth <- qgamma(shape = 4.05, rate = 0.74, p = c(0.5, 0.95))
out_unaddressed$deviation <- out_unaddressed$est - out_unaddressed$truth
out_unaddressed$true_in_CI <- ifelse(out_unaddressed$truth >= 
                                       out_unaddressed$lower_CI
                                   & out_unaddressed$truth <= 
                                     out_unaddressed$upper_CI,
                                   1, 0)
out_unaddressed$quantile <- c(0.5, 0.95)

out <- rbind(out_addressed, out_unaddressed)
out$scenario <- scenario
out$frac_included_data <- nrow(truncated_data)/N

# Save diagnostic plots for five iterations only
if(i < 6){
  
  # Save the plots
  addressing_truncation$plot_running_quantiles / 
                 addressing_truncation$plot_parameters
  ggsave(paste("results_sim/Diagnostic_plots_addr_trunc_FinalFinalsim_scenario_", 
               scenario, "_Iteration", i, "_", "_on_", date_of_today, 
               ".jpeg", sep = ""))
}

# Return the output
return(out)

}

# Run the simulation
result <- mclapply(1:n_runs, 
              function(i) Simulate_truncation(i, scenario =
                                                scenario), mc.cores = cores)

saveRDS( result, paste("results_sim/FinalFinalSim_scenario_", scenario, 
                      "_on_", 
                             date_of_today, ".RDS", sep = "") )

