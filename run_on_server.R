# Run on the server
run_name <- "estimates_vndata"
date_of_today <- "150324"

# Load functions
source("fun_estimate_dic.R")

# Load packages
require(dplyr)
require(tidyr)
require(ggplot2)
require(magrittr)
require(flexsurv)
require(coarseDataTools)
require(rjags)
require(coda)
require(optimParallel)
require(pracma)
require(tibble)
require(patchwork)

# Select the data
dat_lt1 <- readRDS("data/dat_lt1.rds")
dat_lt2 <- readRDS("data/dat_lt2.rds")

# Subset of narrow exposure windows
dat_lt1_narrow <- dat_lt1 %>% filter(L1 - L0 < 5)
dat_lt2_narrow <- dat_lt2 %>% filter(L1 - L0 < 5)

list_data_sets_strict <- list(full = dat_lt1,
                              narrow_only = dat_lt1_narrow)
list_data_sets_loose <- list(full = dat_lt2,
                             narrow_only = dat_lt2_narrow)

# . . . . . Settings: run locally and save . . . . . . . . . . . . . . . . . . .

#JAGS-model specific settings
#    (i) Less thinning is needed when there is almost no autocorrelation (i.e. no
#        pattern in the parameter plot)
#    (i) For the final estimates, try to have 5000 iterations left after thinning.
#    (i) Burnin_period. Note that burnin has to be 
#        shorter than iter/thin.
# 
# burnin_period = 10000
# iterations = 500000
# thinning = 10
# further_thin_plots = 1
# 
# # Data set
# data_set = c("full", "narrow_only")
# 
# # Parametric assumptions
# exposure_window_definition = c("strict", "loose")
# infection_risk_distribution = c("constant", "exp_growth")
# r_SE = 0.004466976
# r_mean = 0.1058827
# 
# method = c("gamma", "GenGamma", "Weibull")
# 
# # Truncation
# right_truncation = c(T, F)
# 
# opt <- expand_grid(burnin_period, iterations, thinning,
#                    exposure_window_definition, data_set,
#                    infection_risk_distribution, r_mean, r_SE,
#                    method, right_truncation, further_thin_plots)
# # Remove the dupicates for constant (r does not matter)
# 
# setwd("~/Library/CloudStorage/OneDrive-UniversiteitLeiden/To_ALICE/Run_files")
# saveRDS(opt, paste("results/Options_", run_name, "_on_", date_of_today, ".RDS", sep = "" ))

# . . . . . Run the model . . . . . . . . . . . . . . . . . . . . . . . . . . .

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
run <- as.numeric(taskid)

opt <- readRDS(paste("results/Options_", run_name, "_on_", date_of_today, ".RDS", sep = "" ))

temp_val <- opt$iterations[run]/opt$thinning[run]

if( (opt$burnin_period[run] >= temp_val) == T){
  cat("Burnin period is longer than the ratio iterations/thinning, 
      which is not allowed. Please adjust.")
  break
  }
  
  cat("Starting iteration number ", run, "/", nrow(opt), ".")
  
  set.seed(1004)
  
tmp <- try(Estimate_doubly_interval_censored(
  
  # Data set
  dat = ifelse(opt$exposure_window_definition[run] == "strict", 
               list_data_sets_strict[ opt$data_set[run] ], 
               list_data_sets_loose[ opt$data_set[run] ]) %>% as.data.frame(), 
  
  # Parametric assumptions
  infection_risk_distribution = opt$infection_risk_distribution[run], 
  exp_growth_rate = opt$r_mean[run], 
  exp_growth_rate_SE =  opt$r_SE[run],
  method =  opt$method[run], 
  
  # Truncation
  right_truncation =  opt$right_truncation[run], 

  # Other settings
  percentiles = c(0.5, 0.9, 0.95, 0.99), 
  iters =  opt$iterations[run], 
  burnin_period =  opt$burnin_period[run], 
  thin =  opt$thinning[run], 
  further_thin_plots =  opt$further_thin_plots[run]
  
), silent = T)

if(class(tmp) == "try-error"){
cat("This run gave an error")
break
}

# . . . . . Save the output. . . . . . . . . . . . . . . . . . . . . . . . . . .

if(class(tmp)!= "try-error"){
# Save the estimates
saveRDS(tmp$estimates, paste("results/Estimates_Run", run, "_", run_name, "_on_", 
                             date_of_today, ".RDS", sep = ""))

# Save the plots
tmp$plot_running_quantiles / tmp$plot_parameters
ggsave(paste("results/Diagnostic_plots_Run", run, "_", run_name, "_on_", date_of_today, 
           ".jpeg", sep = ""))

}

