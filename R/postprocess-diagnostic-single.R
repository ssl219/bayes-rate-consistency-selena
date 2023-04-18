#' Make convergence diagnostic statistics
#'
#' @param fit A fitted CmdStanModel
#' @param outdir Directory to save outputs
#'
#' @return A summary table of estimates and diagnostics
make_convergence_diagnostic_stats <- function(fit, outdir = NA) {

  # Rhat and effective sample size
  fit_summary <- fit$summary(variables = NULL, "rhat", "ess_bulk")

  cat("\n The minimum and maximum effective sample size are ", range(fit_summary$ess_bulk, na.rm = TRUE))
  cat("\n The minimum and maximum Rhat are ", range(fit_summary$rhat, na.rm = TRUE))
  if(min(fit_summary$ess_bulk, na.rm = TRUE) < 500) cat('\n Minimum effective sample size smaller than 500')

  # Diagnostics
  sampler_diagnostics <- fit$diagnostic_summary()

  # # Compute WAIC and LOO
  # tryCatch({
  #   log_lik <- fit$draws("log_lik", format = "matrix")
  # 
  #   n_inf_log_lik <- sum(is.infinite(log_lik))
  #   if(n_inf_log_lik > 0){
  #     .message <- paste("Detected", n_inf_log_lik, "Inf values in log_lik. Removing those iterations.")
  #     warning(.message)
  #     log_lik[is.infinite(log_lik)] <- NA
  #   }
  #   log_lik <- na.omit(log_lik)
  #   WAIC <- loo::waic(log_lik)
  #   LOO <- loo::loo(log_lik)
  # }, error = function(e) e)

  # Time of execution
  time <- fit$time()

  # save
  if(!is.na(outdir)){
    saveRDS(fit_summary, file = file.path(outdir, "fit_summary.rds"))
    # saveRDS(WAIC, file = file.path(outdir, "WAIC.rds"))
    # saveRDS(LOO, file = file.path(outdir, "LOO.rds"))
    saveRDS(sampler_diagnostics, file = file.path(outdir, "sampler_diagnostics.rds"))
    saveRDS(time, file = file.path(outdir, "time_elapsed.rds"))
  } else {
    warning("\n outdir is not given. Results were not saved.")
  }

  return(fit_summary)
}


#' Extract posterior predicted quantities
#'
#' @param fit A fitted CmdStanModel
#' @param dt The dataset used to train the model (data.table)
#' @param predict_type Type of prediction to extract [yhat, yhat_strata]
#'
#' @return A data.table with predictions
#' @export
#'
#' @examples
#' \dontrun{
#'  fit <- readRDS("path/to/model.rds")
#'  dt.po <- extract_posterior_predictions(fit, predict_type = "stratified")
#' }
extract_posterior_predictions <- function(fit, dt){
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')

  # Predicted contacts
  po <- fit$draws("yhat", inc_warmup = FALSE, format = "draws_matrix")
  dt.po <- as.data.table(reshape2::melt(po))

  # Extract indices
  .pattern <- "yhat\\[([0-9]+),([0-9]+),([0-9]+)\\]"

  dt.po[, comb_idx := as.numeric(gsub(.pattern, "\\1", variable))]
  dt.po[, age_idx := as.numeric(gsub(.pattern, "\\2", variable))]
  dt.po[, alter_age_idx := as.numeric(gsub(.pattern, "\\3", variable))]

  # Calculate quantiles
  dt.po <- dt.po[, list( q=quantile(value, prob=ps, na.rm=T), q_label = p_labs),
                 by=list(comb_idx, age_idx, alter_age_idx)]
  dt.po <- as.data.table(
    dcast(dt.po, comb_idx + age_idx + alter_age_idx ~ q_label,
          value.var = "q")
  )

  # Recover gender and alter gender
  dt.po[, gender := fcase(
    comb_idx %in% c(1,3), "Male",
    comb_idx %in% c(2,4), "Female",
    default = NA)]

  dt.po[, alter_gender := fcase(
    comb_idx %in% c(1,4), "Male",
    comb_idx %in% c(2,3), "Female",
    default = NA)]

  # Recover age
  dt.po[, age := unique(dt$age)[age_idx]]
  dt.po[, alter_age := unique(dt$alter_age)[alter_age_idx]]
  dt <- merge(dt, dt.po[, list(age, alter_age, gender, alter_gender, CL, CU, M)],
              by=c("age", "alter_age", "gender", "alter_gender"))

  # Calculate contact intensities and rates
  dt[, cntct_intensity_predict := M/part]
  dt[, cntct_intensity_predict_CL := CL/part]
  dt[, cntct_intensity_predict_CU := CU/part]

  dt[, cntct_rate_predict := M/part/pop]
  dt[, cntct_rate_predict_CL := CL/part/pop]
  dt[, cntct_rate_predict_CU := CU/part/pop]

  return(dt)
}

make_ppd_check_covimod <- function(po, dt.survey, outdir=NA, new_hh=FALSE){
  if (new_hh){
    ps <- c(0.5, 0.025, 0.975)
    p_labs <- c('M','CL','CU')
    
    po_MM <- subset(po, "yhat_strata_MM")
    po_FF <- subset(po, "yhat_strata_FF")
    po_MF <- subset(po, "yhat_strata_MF")
    po_FM <- subset(po, "yhat_strata_FM")
  
    dt.po_MM <- as.data.table(reshape2::melt(po_MM))
    dt.po_FF <- as.data.table(reshape2::melt(po_FF))
    dt.po_MF <- as.data.table(reshape2::melt(po_MF))
    dt.po_FM <- as.data.table(reshape2::melt(po_FM))
  
    # Extract indices  
    .patternMM <- "yhat_strata_MM\\[([0-9]+),([0-9]+)\\]"
    .patternFF <- "yhat_strata_FF\\[([0-9]+),([0-9]+)\\]"
    .patternMF <- "yhat_strata_MF\\[([0-9]+),([0-9]+)\\]"
    .patternFM <- "yhat_strata_FM\\[([0-9]+),([0-9]+)\\]"
    
    dt.po_MM[, comb_idx := 1L]
    dt.po_MM[, part_idx := as.numeric(gsub(.patternMM, "\\1", variable))]
    dt.po_MM[, alter_age_strata_idx := as.numeric(gsub(.patternMM, "\\2", variable))]
    
    dt.po_FF[, comb_idx := 2L]
    dt.po_FF[, part_idx := as.numeric(gsub(.patternFF, "\\1", variable))]
    dt.po_FF[, alter_age_strata_idx := as.numeric(gsub(.patternFF, "\\2", variable))]
    
    dt.po_MF[, comb_idx := 3L]
    dt.po_MF[, part_idx := as.numeric(gsub(.patternMF, "\\1", variable))]
    dt.po_MF[, alter_age_strata_idx := as.numeric(gsub(.patternMF, "\\2", variable))]
    
    dt.po_FM[, comb_idx := 4L]
    dt.po_FM[, part_idx := as.numeric(gsub(.patternFM, "\\1", variable))]
    dt.po_FM[, alter_age_strata_idx := as.numeric(gsub(.patternFM, "\\2", variable))]
    
    # not needed
    # # Recover participant sizes
    # P_MM <- max(dt.po_MM$part_idx)
    # P_FF <- max(dt.po_FF$part_idx)
    # P_MF <- max(dt.po_MF$part_idx)
    # P_FM <- max(dt.po_FM$part_idx)
    # 
    # # Adjust participant index part_idx to merge everything in one dataset
    # dt.po_FF[, part_idx := part_idx + P_MM]
    # dt.po_MF[, part_idx := part_idx + P_MM + P_FF]
    # dt.po_FM[, part_idx := part_idx + P_MM + P_FF + P_MF]
    
    # Bind everything
    dt.po <- rbind(dt.po_MM, dt.po_FF, dt.po_MF, dt.po_FM)
    
    # Calculate quantiles
    dt.po <- dt.po[, list( q=quantile(value, prob=ps, na.rm=T), q_label = p_labs),
                   by=list(comb_idx, part_idx, alter_age_strata_idx)]
    dt.po <- data.table::dcast(dt.po, comb_idx + part_idx + alter_age_strata_idx ~ q_label, value.var = "q")
    
    # Recover gender and alter gender
    dt.po[, gender := fcase(
      comb_idx == 1, "Male",
      comb_idx == 2, "Female",
      comb_idx == 3, "Male",
      comb_idx == 4, "Female",
      default = NA)]
    
    dt.po[, alter_gender := fcase(
      comb_idx == 1, "Male",
      comb_idx == 2, "Female",
      comb_idx == 3, "Female",
      comb_idx == 4, "Male",
      default = NA)]
    
    # order survey in same ordering as data preprocessing
    covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
    # making sure order of factors in alter_age_strata is ascending instead of decreasing
    # note alter_age_strata_idx and age_strata_idx are the same, but they serve different purposes
    dt.survey[, age_strata_idx:=as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
    dt.survey <- dt.survey[order(age, new_id, age_strata_idx, gender, alter_gender)]
    
    # Recover new_id and alter_age_strata (have to make sure dt.survey is well ordered as done before)
    dt.po[, alter_age_strata := sort(unique(dt.survey$alter_age_strata))[alter_age_strata_idx]]
    dt.po[, new_id := fcase(
      comb_idx == 1, unique(dt.survey[gender=="Male" & alter_gender=="Male"]$new_id)[part_idx],
      comb_idx == 2, unique(dt.survey[gender=="Female" & alter_gender=="Female"]$new_id)[part_idx],
      comb_idx == 3, unique(dt.survey[gender=="Male" & alter_gender=="Female"]$new_id)[part_idx],
      comb_idx == 4, unique(dt.survey[gender=="Female" & alter_gender=="Male"]$new_id)[part_idx],
      default = NA
    )]
    
    # merge with contact data in dt.survey
    dt <- merge(dt.survey, dt.po[, list(new_id, alter_age_strata, gender, alter_gender, CL, CU, M)],
                by=c("new_id", "alter_age_strata", "gender", "alter_gender"))
    
    dt[, inside.CI := y >= CL & y <= CU]
    proportion_ppd <- mean(dt$inside.CI, na.rm=T)
    cat(" Proportion of points within posterior predictive 95% CI: ", proportion_ppd, "\n")
    
    if(!is.na(outdir)){
      saveRDS(dt, file.path(outdir, "ppd_check.rds"))
      saveRDS(proportion_ppd, file.path(outdir, "ppd_check_proportion.rds"))
    } else {
      warning("\n outdir is not specified. Results were not saved.")
    }
  }
  else{
    ps <- c(0.5, 0.025, 0.975)
    p_labs <- c('M','CL','CU')
    
    po <- subset(po, "yhat_strata")
    
    dt.po <- as.data.table(reshape2::melt(po))
    
    # Extract indices
    .pattern <- "yhat_strata\\[([0-9]+),([0-9]+), ([0-9]+)\\]"
    
    dt.po[, comb_idx := as.numeric(gsub(.pattern, "\\1", variable))]
    dt.po[, age_idx := as.numeric(gsub(.pattern, "\\2", variable))]
    dt.po[, alter_age_strata_idx := as.numeric(gsub(.pattern, "\\3", variable))]
    
    # Calculate quantiles
    dt.po <- dt.po[, list( q=quantile(value, prob=ps, na.rm=T), q_label = p_labs),
                   by=list(comb_idx, age_idx, alter_age_strata_idx)]
    dt.po <- data.table::dcast(dt.po, comb_idx + age_idx + alter_age_strata_idx ~ q_label, value.var = "q")
    
    # Recover gender and alter gender
    dt.po[, gender := fcase(
      comb_idx %in% c(1,3), "Male",
      comb_idx %in% c(2,4), "Female",
      default = NA)]
    
    dt.po[, alter_gender := fcase(
      comb_idx %in% c(1,4), "Male",
      comb_idx %in% c(2,3), "Female",
      default = NA)]
    
    # Recover age
    dt.po[, age := unique(dt.survey$age)[age_idx]]
    dt.po[, alter_age_strata := sort(unique(dt.survey$alter_age_strata))[alter_age_strata_idx]]
    dt <- merge(dt.survey, dt.po[, list(age, alter_age_strata, gender, alter_gender, CL, CU, M)],
                by=c("age", "alter_age_strata", "gender", "alter_gender"))
    
    dt[, inside.CI := y >= CL & y <= CU]
    proportion_ppd <- mean(dt$inside.CI, na.rm=T)
    cat(" Proportion of points within posterior predictive 95% CI: ", proportion_ppd, "\n")
    
    if(!is.na(outdir)){
      saveRDS(dt, file.path(outdir, "ppd_check.rds"))
      saveRDS(proportion_ppd, file.path(outdir, "ppd_check_proportion.rds"))
    } else {
      warning("\n outdir is not specified. Results were not saved.")
    }
  }

  return(dt)
}

extract_posterior_rates <- function(po){
  po <- subset(po, "log_cnt_rate")
  dt.po <- as.data.table(reshape2::melt(po))

  # Extract indices
  .pattern <- "log_cnt_rate\\[([0-9]+),([0-9]+),([0-9]+)\\]"

  dt.po[, comb_idx := as.numeric(gsub(.pattern, "\\1", variable))]
  dt.po[, age_idx := as.numeric(gsub(.pattern, "\\2", variable))]
  dt.po[, alter_age_idx := as.numeric(gsub(.pattern, "\\3", variable))]
  
  # Calculate posterior contact rate
  # dt.po[, value := exp(value)]

  return(dt.po)
}

extract_posterior_alpha <- function(po, gender_comb="MM"){
  # po <- subset(po, "alpha_age")
  dt.po <- as.data.table(reshape2::melt(po))
  
  # Extract indices
  if (gender_comb=="MM"){
    .pattern <- "alpha_age_MM\\[([0-9]+),([0-9]+)\\]"
  }
  
  if (gender_comb=="FF"){
    .pattern <- "alpha_age_FF\\[([0-9]+),([0-9]+)\\]"
  }
  
  if (gender_comb=="MF"){
    .pattern <- "alpha_age_MF\\[([0-9]+),([0-9]+)\\]"
  }
  
  if (gender_comb=="FM"){
    .pattern <- "alpha_age_FM\\[([0-9]+),([0-9]+)\\]"
  }
  
  
  dt.po[, part_idx := as.numeric(gsub(.pattern, "\\1", variable))]
  dt.po[, alter_age_idx := as.numeric(gsub(.pattern, "\\2", variable))]
  
  return(dt.po)
}

posterior_alpha <- function(fit, dt.po, type="matrix", outdir=NA, gender_comb="MM"){
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(type=="matrix"){ # Full contact intensity matrix
    # Calculate quantiles
    dt.po <- dt.po[, list(q=quantile(value, prob=ps, na.rm=T), q_label = p_labs),
                   by=list(part_idx, alter_age_idx)]
    dt.po <- data.table::dcast(dt.po, part_idx + alter_age_idx ~ q_label, value.var = "q")
    
    # If COVIMOD data
    dt.po[, alter_age := alter_age_idx - 1]
    
    # add participant ages
    
    # unique number of participants
    N_part <- max(dt.po$part_idx)
    # extract participant ages and add gender and alter_gender
    if (gender_comb=="MM"){
      dt.po[, gender:="Male"]
      dt.po[, alter_gender:="Male"]
      dt.part_age <- as.data.table(reshape2::melt(fit$draws(c("part_age_MM"),inc_warmup = FALSE, format="draws_matrix")[1, c(1:N_part)]))
    }
    
    if (gender_comb=="FF"){
      dt.po[, gender:="Female"]
      dt.po[, alter_gender:="Female"]
      dt.part_age <- as.data.table(reshape2::melt(fit$draws(c("part_age_FF"),inc_warmup = FALSE, format="draws_matrix")[1, c(1:N_part)]))
    }
    
    if (gender_comb=="MF"){
      dt.po[, gender:="Male"]
      dt.po[, alter_gender:="Female"]
      dt.part_age <- as.data.table(reshape2::melt(fit$draws(c("part_age_MF"),inc_warmup = FALSE, format="draws_matrix")[1, c(1:N_part)]))
    }
    
    if (gender_comb=="FM"){
      dt.po[, gender:="Female"]
      dt.po[, alter_gender:="Male"]
      dt.part_age <- as.data.table(reshape2::melt(fit$draws(c("part_age_FM"),inc_warmup = FALSE, format="draws_matrix")[1, c(1:N_part)]))
    }
    
    # create smaller dataset with id and age
    id_age_alpha <- data.table(part_idx=1:N_part, age=dt.part_age$value)
    
    # add to dataset
    dt.po <- merge(dt.po, id_age_alpha, by=c("part_idx"))
    
    # aggregate by age and alter_age
    group_var <- c("age", "alter_age")
    dt.po[, alpha_agg := mean(M), by=group_var]
    
    # save
    
    if(!is.na(outdir)){
      if (gender_comb=="MM"){
        saveRDS(dt.po, file.path(outdir, "alphaMM_matrix.rds"))
      }
      else if (gender_comb=="FF"){
        saveRDS(dt.po, file.path(outdir, "alphaFF_matrix.rds"))
      }
      else if (gender_comb=="MF"){
        saveRDS(dt.po, file.path(outdir, "alphaMF_matrix.rds"))
      }
      else if (gender_comb=="FM"){
        saveRDS(dt.po, file.path(outdir, "alphaFM_matrix.rds"))
      }
    }
    
    else{
      warning("\n outdir is not specified. Results were not saved.")
    }
    
    return(dt.po)
  }
  
}  


posterior_contact_intensity <- function(dt.po, dt.pop, type="matrix", simulation=FALSE, outdir=NA, new_hh=FALSE){
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')

  if(type=="matrix"){ # Full contact intensity matrix
    # Calculate quantiles
    dt.po <- dt.po[, list(q=quantile(value, prob=ps, na.rm=T), q_label = p_labs),
                   by=list(comb_idx, age_idx, alter_age_idx)]
    dt.po <- data.table::dcast(dt.po, comb_idx + age_idx + alter_age_idx ~ q_label, value.var = "q")

    # Convert back to rates
    dt.po[, M := exp(M)]
    dt.po[, CL := exp(CL)]
    dt.po[, CU := exp(CU)]

    # Recover age
    if(simulation){ # If simulation data
      dt.po[, age := age_idx + 5]
      dt.po[, alter_age := alter_age_idx + 5]
    } else { # If COVIMOD data
      dt.po[, age := age_idx - 1]
      dt.po[, alter_age := alter_age_idx - 1]
    }

    # Recover gender and alter gender
    dt.po[, gender := fcase(comb_idx %in% c(1,3), "Male",
                            comb_idx %in% c(2,4), "Female", default = NA)]
    dt.po[, alter_gender := fcase(comb_idx %in% c(1,4), "Male",
                                  comb_idx %in% c(2,3), "Female", default = NA)]
    
    if (new_hh){
      if(!is.na(outdir)){
        saveRDS(dt.po, file.path(outdir, "rate_matrix.rds"))
      } else {
        warning("\n outdir is not specified. Results were not saved.")
      }
    }
    
    else{
      # Load datasets
      dtp <- copy(dt.pop)
      setnames(dtp, c("age", "gender"), c("alter_age", "alter_gender"))
      dt.po <- merge(dt.po, dtp, by=c("alter_age", "alter_gender"), all.x = TRUE)
      
      dt.po[, intensity_M := M * pop]
      dt.po[, intensity_CL := CL * pop]
      dt.po[, intensity_CU := CU * pop]
      
      if(!is.na(outdir)){
        saveRDS(dt.po, file.path(outdir, "intensity_matrix.rds"))
      } 
      else {
        warning("\n outdir is not specified. Results were not saved.")
      }
    }
    
    return(dt.po)
  } else { # Marginal contact intensity by gender

    dtm <- dt.po[comb_idx == 1 | comb_idx == 3]
    if (simulation){
      dtm[, age := age_idx + 5]
    } else {
      dtm[, age := age_idx - 1]
    }

    dtm <- merge(dtm, dt.pop[gender == "Male"], by="age")
    dtm <- dtm[, .(value = sum(exp(value + log(pop)))), by=c("draw", "age")]
    dtm <- dtm[, list( q=quantile(value, prob=ps, na.rm=T), q_label = p_labs), by="age"]
    dtm[, gender := "Male"]
    gc()

    dtf <- dt.po[comb_idx == 2 | comb_idx == 4]
    if (simulation){
      dtf[, age := age_idx + 5]
    } else {
      dtf[, age := age_idx - 1]
    }

    dtf <- merge(dtf, dt.pop[gender == "Female"], by="age")
    dtf <- dtf[, .(value = sum(exp(value + log(pop)))), by=c("draw", "age")]
    dtf <- dtf[, list( q=quantile(value, prob=ps, na.rm=T), q_label = p_labs), by="age"]
    dtf[, gender := "Female"]
    gc()

    dt <- rbind(dtm, dtf)
    dt <- data.table::dcast(dt, age + gender ~ q_label, value.var = "q")
    setnames(dt, c("M", "CL", "CU"), c("intensity_M", "intensity_CL", "intensity_CU"))

    if(!is.na(outdir)){
      saveRDS(dt, file.path(outdir, "intensity_marginal.rds"))
    } else {
      warning("\n outdir is not specified. Results were not saved.")
    }

    return(dt)
  }
}

#' Makes a table for posterior predictive checks
#'
#' @param dt Output of `extract_posterior_predictions()`
#' @param predict_type Type of prediction to extract [yhat, yhat_strata]
#' @param outdir
#'
#' @return A data.table with indications for whether the predicted estimates lie in the posterior predictive 95% CI
#' @export
#'
#' @examples
#' \dontrun{
#' dt.po <- extract_posterior_predictions(fit, predict_type = "yhat")
#'
#' make_posterior_predictive_check(dt.po, predict_type="yhat")
#' }
make_posterior_predictive_check <- function(dt, outdir=NA){

  dt[, inside.CI := cntct_intensity <= cntct_intensity_predict_CU & cntct_intensity >= cntct_intensity_predict_CL]
  cat("\n Proportion of points within posterior predictive 95% CI: ", mean(dt$inside.CI, na.rm=T))

  if(!is.na(outdir)){
    saveRDS(dt, file.path(outdir, "ppd_check.rds"))
  } else {
    warning("\n outdir is not specified. Results were not saved.")
  }

  return(dt)
}

#' Make MSE summary stats table
#'
#' @param dt Output of `extract_posterior_predictions()`
#'
#' @return data.frame with squared bias and MSE values
#' @export
#'
#' @examples
#' \dontrun{
#' dt.po <- extract_posterior_predictions(fit)
#'
#' make_mse_table(dt.po)
#' }
make_error_table <- function(dt, outdir=NA){
  mse <- function(y, y_pred) mean( (y - y_pred)**2, na.rm=T )
  sbias <- function(y, y_pred) mean( y - y_pred, na.rm=T )^2

  df <- data.frame(
    metric = c("bias", "bias", "mse", "mse"),
    name = c("intensity", "rate", "intensity", "rate"),
    value = c(sbias(dt$cntct_intensity ,dt$cntct_intensity_predict),
              sbias(dt$cntct_rate, dt$cntct_rate_predict),
              mse(dt$cntct_intensity ,dt$cntct_intensity_predict),
              mse(dt$cntct_rate, dt$cntct_rate_predict))
  )

  if(!is.na(outdir)){
    saveRDS(df, file = file.path(outdir, "mse.rds"))
  } else {
    return(df)
  }
}
