###
## overview:  this script provides a power detection of MR methods on large-scale simulation GWAS data
## date:      02/2025
## by:        Jiyuan Jiao (jiyuanj@umich.edu)
###

### load required packages and set parameters ###
setwd("~/MR_GWAS_Simulation")


## Libraries
# install.packages("devtools")
library(devtools)
# devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
# devtools::install_github("tye27/mr.divw")
library(mr.divw)
# devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
source("Power_Fun.R")

## Hyper Parameters setup
M <- 461    # Number of SNPs
N <- 75000  # Sample size
gamma <- 0.1    # True Causal effect of X on Y
f_jX <- rbeta(461, 2, 5)     # Effect allele frequencies for Trait X
f_jY <- runif(M, min = 0.05, max = 0.8)  # Effect allele frequencies for Trait Y
match <- TRUE # Matched sign pleiotropy
num_sim <- 5000  # Number of Simulations


## Data Loading and Setting Up
sbp_dat <- read.csv("sbp_data.csv")
beta_Xj <- sbp_dat$beta.exposure
s_Xj <- sbp_dat$se.exposure
beta_Xj_hat <- rnorm(M, mean = beta_Xj, sd = s_Xj)  # beta.exp with Gaussian Noise
s_Yj <- sqrt(1 / (2 * f_jY * (1 - f_jY) * N))  # Standard error of beta.out
ml <- mr_method_list() |> dplyr::filter(heterogeneity_test == TRUE)
summary((beta_Xj_hat/(s_Xj*sqrt(sbp_dat$samplesize.exposure)))^2)

ple_mean <- c(0, 0.001, 0.005, 0.01, 0.05, 0.1) # Pleiotropy mean
ple_var <- c(0, 0.001, 0.003, 0.005) # Pleiotropy variance

type <- "Q"

### Power Calculation
## Q Style Pleiotropy Data
type <- "Q"

# Variance Scenario 1 (Mean 0, Variance 0.000)
s1_result <- power(num_sim, paste0("rs", 1:M), f_jY, beta_Xj_hat, s_Xj, s_Yj, 
                   gamma = gamma, ple_mean[1], ple_var[1], match, type) 

# Variance Scenario 2 (Mean 0, Variance 0.001)
s2_result <- power(num_sim, paste0("rs", 1:M), f_jY, beta_Xj_hat, s_Xj, s_Yj, 
                   gamma = gamma, ple_mean[1], ple_var[2], match, type)

# Variance Scenario 3 (Mean 0, Variance 0.003)
s3_result <- power(num_sim, paste0("rs", 1:M), f_jY, beta_Xj_hat, s_Xj, s_Yj, 
                   gamma = gamma, ple_mean[1], ple_var[3], match, type)

# Variance Scenario 4 (Mean 0, Variance 0.005)
s4_result <- power(num_sim, paste0("rs", 1:M), f_jY, beta_Xj_hat, s_Xj, s_Yj, 
                   gamma = gamma, ple_mean[1], ple_var[4], match, type)

# Save result for each cycle of simulation
for (i in 1:length(ple_var)) {
  df_name <- paste0("s", i, "_result")
  var_name <- paste0("s", i, "_result_Qdata")
  file_name <- paste0("./results/", var_name, ".csv")
  write.csv(get(df_name), file_name)
}

## PRESSO Style Pleiotropy Data
type <- "PRESSO"

# Variance Scenario 1 (Mean 0, Variance 0.000)
v1_result <- power(num_sim, paste0("rs", 1:M), f_jY, beta_Xj_hat, s_Xj, s_Yj, 
                   gamma = gamma, ple_mean[1], ple_var[1], match, type) 

# Variance Scenario 2 (Mean 0, Variance 0.001)
v2_result <- power(num_sim, paste0("rs", 1:M), f_jY, beta_Xj_hat, s_Xj, s_Yj, 
                   gamma = gamma, ple_mean[1], ple_var[2], match, type)

# Variance Scenario 3 (Mean 0, Variance 0.003)
v3_result <- power(num_sim, paste0("rs", 1:M), f_jY, beta_Xj_hat, s_Xj, s_Yj, 
                   gamma = gamma, ple_mean[1], ple_var[3], match, type)

# Variance Scenario 4 (Mean 0, Variance 0.005)
v4_result <- power(num_sim, paste0("rs", 1:M), f_jY, beta_Xj_hat, s_Xj, s_Yj, 
                   gamma = gamma, ple_mean[1], ple_var[4], match, type)

# Save result for each cycle of simulation
for (i in 1:length(ple_var)) {
  df_name <- paste0("v", i, "_result")
  var_name <- paste0("v", i, "_result_Qdata")
  file_name <- paste0("./results/", var_name, ".csv")
  write.csv(get(df_name), file_name)
}

