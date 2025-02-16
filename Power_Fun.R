q_gwas_data_fun <- function(SNP, f_j, beta.exposure, se.exposure, se.outcome, gamma=0.1, 
                            ple.mean, ple.se, match_sign = TRUE){
  
  u_j <- rnorm(length(beta.exposure), mean=ple.mean, sd=ple.se)
  if(match_sign){
    u_j <- abs(u_j)*sign(beta.exposure)
  }
  beta_Yj <- gamma * beta.exposure + u_j
  beta.outcome <- rnorm(length(beta.exposure), mean = beta_Yj, sd = se.outcome)
  
  temp_df <- data.frame(
    SNP = SNP,
    f_j = f_j,
    beta.exposure = beta.exposure,
    se.exposure = se.exposure,
    beta.outcome = beta.outcome,
    se.outcome = se.outcome,
    id.exposure = "Trait X",  
    exposure = "X",
    id.outcome = "Trait Y",
    outcome = "Y",
    mr_keep = TRUE)
  
  return(temp_df)
}

presso_gwas_data_fun <- function(SNP, f_j, beta.exposure, se.exposure, se.outcome, gamma, random_indice, 
                                 mean, var, match_sign = TRUE){
  
  u_j <- rep(0, length(beta.exposure))
  u_j[random_indice] <- rnorm(length(random_indice), mean = mean, sd = var)
  if(match_sign){
    u_j <- abs(u_j)*sign(beta.exposure)
  }
  beta_Yj <- gamma * beta_Xj + u_j
  beta.outcome <- rnorm(M, mean = beta_Yj, sd = se.outcome)
  
  temp_df <- data.frame(
    SNP = SNP,
    f_j = f_j,
    beta.exposure = beta.exposure,
    se.exposure = se.exposure,
    beta.outcome = beta.outcome,
    se.outcome = se.outcome,
    id.exposure = "Trait X",  
    exposure = "X",
    id.outcome = "Trait Y",
    outcome = "Y",
    pl_effect = u_j,
    mr_keep = TRUE)
  
  return(temp_df)
}

power <- function(num_sim, SNP, f_j, beta.exposure, se.exposure, se.outcome, gamma=0.1, 
                  ple.mean, ple.se, match_sign = TRUE, type){
  
  
  q_test_pvalues <- numeric(num_sim)
  egger_test_pvalues <- numeric(num_sim)
  qivwr_test_pvalues <- numeric(num_sim)
  presso_test_pvalues <- numeric(num_sim)
  
  ivw_est <- numeric(num_sim)
  weighted_median_est <- numeric(num_sim)
  egger_est <- numeric(num_sim)
  weighted_mode_est <- numeric(num_sim)
  raps_est <- numeric(num_sim)
  ivwr_est <- numeric(num_sim)
  
  ivw_bias <- numeric(num_sim)
  weighted_median_bias <- numeric(num_sim)
  egger_bias <- numeric(num_sim)
  weighted_mode_bias <- numeric(num_sim)
  raps_bias <- numeric(num_sim)
  ivwr_bias <- numeric(num_sim)
  
  for (i in 1:num_sim) {
    
    if(type == "Q"){
      df <- q_gwas_data_fun(SNP, f_j, beta.exposure, se.exposure, se.outcome, gamma = gamma, 
                            ple.mean, ple.se, match_sign)
      
    }else if(type == "PRESSO"){
      random_indice <- sample(1:M, n_out)
      df <- presso_gwas_data_fun(SNP, f_j, beta.exposure, se.exposure, se.outcome, gamma = gamma, 
                                 random_indice, ple.mean, ple.se, match_sign)
      
    }else{
      stop("Must enter a valid type of pleiotropy style.")
      
    }
    
    q_test_pvalues[i] <- mr_heterogeneity(df)$Q_pval[2]
    egger_test_pvalues[i] <- mr_pleiotropy_test(df)$pval
    qivwr_test_pvalues[i] <- mr_heterogeneity(df, method_list = ml$obj[4])$Q_pval
    invisible(capture.output(suppressMessages(suppressWarnings({est <- mr(df, method_list = c("mr_ivw", "mr_weighted_median", 
                                                                                              "mr_egger_regression", "mr_weighted_mode", "mr_ivw_radial"))}))))
    
    ivw_est[i] <- est$b[est$method=="Inverse variance weighted"]
    ivwr_est[i] <- est$b[est$method=="IVW radial"]
    weighted_median_est[i] <- est$b[est$method=="Weighted median"]
    egger_est[i] <- est$b[est$method=="MR Egger"]
    weighted_mode_est[i] <- est$b[est$method=="Weighted mode"]
    # stop("error encounterd")
    # raps_est[i] <- mr.raps::mr.raps(df)$beta.hat
    raps_est[i] <- mr.raps::mr.raps(b_exp=df$beta.exposure, b_out=df$beta.outcome,
                                    se_exp=df$se.exposure, se_out=df$se.outcome)$beta.hat
    # raps_est[i] <- est$b[est$method=="Robust adjusted profile score (RAPS)"]
    
    ivw_bias[i] <- abs(est$b[est$method=="Inverse variance weighted"] - gamma)
    ivwr_bias[i] <- abs(est$b[est$method=="IVW radial"] - gamma)
    weighted_median_bias[i] <- abs(est$b[est$method=="Weighted median"] - gamma)
    egger_bias[i] <- abs(est$b[est$method=="MR Egger"] - gamma)
    weighted_mode_bias[i] <- abs(est$b[est$method=="Weighted mode"] - gamma)
    raps_bias[i] <- abs(raps_est[i] - gamma)
    # presso_test_pvalues[i] <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
    #                                     SdOutcome = "se.outcome", SdExposure = "se.exposure", data = df, 
    #                                     OUTLIERtest = T, DISTORTIONtest = T)[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
  }
  
  power_Q_test <- mean(q_test_pvalues < 0.05)
  power_EGGER_test <- mean(egger_test_pvalues < 0.05)
  power_Qivwr_test <- mean(qivwr_test_pvalues < 0.05)
  avg_ivw_bias <- mean(ivw_bias)
  avg_weighted_median_bias <- mean(weighted_median_bias)
  avg_egger_bias <- mean(egger_bias)
  avg_weighted_mode_bias <- mean(weighted_mode_bias)
  avg_raps_bias <- mean(raps_bias)
  avg_ivwr_bias <- mean(ivwr_bias) 
  
  df_results <- data.frame(
    q_test_pvalue = q_test_pvalues,
    egger_test_pvalue = egger_test_pvalues,
    qivwr_test_pvalue = qivwr_test_pvalues,
    ivw_est = ivw_est,
    ivw_radial_est = ivwr_est,
    weighted_median_est = weighted_median_est,
    egger_est = egger_est,
    weighted_mode_est = weighted_mode_est,
    raps_est = raps_est
  )
  
  cat("Variance", ple.se, "Estimated power of Q test:", power_Q_test, "\n")
  cat("Variance", ple.se, "Estimated power of EGGER:", power_EGGER_test, "\n")
  cat("Variance", ple.se, "Estimated power of Q test radial IVW:", power_Qivwr_test, "\n")
  cat("Variance", ple.se, "Avg Bias for IVW:", avg_ivw_bias, "\n")
  cat("Variance", ple.se, "Avg Bias for IVW Radial:", avg_ivwr_bias, "\n")
  cat("Variance", ple.se, "Avg Bias for Weighted Median:", avg_weighted_median_bias, "\n")
  cat("Variance", ple.se, "Avg Bias for MR-EGGER:", avg_egger_bias, "\n")
  cat("Variance", ple.se, "Avg Bias for Weighted Mode:", avg_weighted_mode_bias, "\n")
  cat("Variance", ple.se, "Avg Bias for RAPS:", avg_raps_bias, "\n")
  # cat("Variance", ple.se, "PRESSO pvals are:", presso_test_pvalues, "\n")
  
  return(df_results)
}
