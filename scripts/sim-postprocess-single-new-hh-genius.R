# Preamble: Generates diagnostic statistics and result plots

cat("\n ---------- Begin Post-processing ---------- \n")

library(optparse)
library(data.table)
library(cmdstanr)
library(bayesplot)
library(loo)
library(reshape2)
library(stringr)
library(ggplot2)
library(viridis)
library(pammtools)

##### ---------- I/O ---------- #####
option_list <- list(
  optparse::make_option("--repo_path", type = "character", default = "/rds/general/user/ssl219/home/bayes-rate-consistency-selena",
                       help = "Absolute file path to repository directory, used as long we don t build an R package [default]",
                       dest = "repo.path"),
  optparse::make_option("--hhsize", type = "integer", default = 4,
                        help = "Household size [default %default]",
                        dest = "hhsize"),
  optparse::make_option("--sample_size", type = "integer", default = 55,
                        help = "Number of participants with random age in the survey [default \"%default\"]",
                        dest = 'sample_size'),
  optparse::make_option("--scenario", type = "character", default = "flat",
                        help = "Scenario [default %default]",
                        dest = "scenario"),
  optparse::make_option("--data_path", type = "character", default = "/rds/general/user/ssl219/home",
                       help = "Absolute file path to data directory, used as long we don t build an R package [default]",
                       dest = 'data.path'),
  # optparse::make_option("--repo_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate/bayes-rate-consistency-selena",
  #                        help = "Absolute file path to repository directory, used as long we don t build an R package [default]",
  #                        dest = "repo.path"),
  # optparse::make_option("--data_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate",
  #                        help = "Absolute file path to data directory, used as long we don t build an R package [default]",
  #                        dest = 'data.path'),
  optparse::make_option("--wave", type="integer", default = 1,
                        help = "COVIMOD wave",
                        dest = "wave"),
  optparse::make_option("--model", type = "character", default = NA_character_,
                        help = "Name of the model",
                        dest = "model.name"),
  optparse::make_option("--mixing", type = "logical", default = TRUE,
                        help = "Whether to assess mixing",
                        dest = "mixing"),
  optparse::make_option("--ppc", type = "logical", default = TRUE,
                        help = "Whether to run posterior predictive checks",
                        dest = "ppc"),
  optparse::make_option("--plot", type = "logical", default = TRUE,
                        help = "Whether to plot posterior distributions",
                        dest = "plot"),
  optparse::make_option("--sim.no", type = "integer", default = 1,
                        help = "Simulated Dataset Number [default %default]",
                        dest = "sim.no")
)

cat("\n before args")

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

cat("\n after args")

# args$repo.path = "/Users/mac/Documents/M4R/code/bayes_consistency_rate/bayes-rate-consistency-selena"
# args$data.path = "/Users/mac/Documents/M4R/code/bayes_consistency_rate"
# args$model.name = "hsgp-eq-rd-new-hh-dropping-all-zeros-symmetric-poisson-1-sim-flat-new-hh-hh2-55"

model.path <- file.path(args$repo.path, "stan_fits", paste0(args$model.name, ".rds"))
data.path <- file.path(args$data.path, "data/simulations/datasets/new-hh-both", paste0("hh", args$hhsize, "-", args$sample_size), paste0("dataset", args$sim.no), paste0("data-hh", args$hhsize, "-", args$scenario, "-", args$sample_size, "-amended-genius.rds"))


# Error handling
if(!file.exists(model.path)) {
  cat("\n Model: ", model.path)
  stop("The specified model does not exist!")
}
if(!file.exists(data.path)) {
  stop("The specified dataset does not exist!")
}

# Output directories
export.path <- file.path(args$repo.path, "results-sim", args$model.name)
export.fig.path <- file.path(export.path, "figures")
if(!dir.exists(export.path)){
  dir.create(export.path, recursive = TRUE)
  dir.create(export.fig.path)
} else {
  if(!dir.exists(export.fig.path)){
    dir.create(export.fig.path)
  }
}

##### ---------- Setup ---------- #####
cat(paste("\n Model path:", model.path))
cat(paste("\n Data path:", data.path))

fit <- readRDS(model.path)
data <- readRDS(data.path)
dt.cnt <- data$contacts[wave == args$wave]
dt.offsets <- data$offsets[wave == args$wave]
dt.pop <- data$pop

source(file.path(args$repo.path, "R/covimod-utility.R"))
source(file.path(args$repo.path, "R/stan-utility.R"))
source(file.path(args$repo.path, "R/postprocess-diagnostic-single.R"))
source(file.path(args$repo.path, "R/postprocess-plotting-single.R"))
source(file.path(args$repo.path, "R", "sim-dataset-utility.R"))

#### ------------------- True contact rates --------------------- #####

if (args$scenario == "flat"){
  dt.sim.true.cntct <- as.data.table(readRDS( file.path(args$data.path, "data/simulations/intensity/new-hh/flat-data.rds") ))
}

if (args$scenario == "boarding_school"){
  dt.sim.true.cntct <- as.data.table(readRDS( file.path(args$data.path, "data/simulations/intensity/new-hh/boarding-data.rds") ))
}

#### ------------------- Stan data --------------------- #####

## Configure Stan data
cat(" Configuring Stan data ...")

# Initialize
stan_data <- init_stan_data(A=55, C=8)

# Add contact counts
stan_data <- add_contact_vector(stan_data, dt.cnt, survey="COVIMOD", single = TRUE, new_hh = TRUE)

# Add obs counts
stan_data <- add_N(stan_data, dt.cnt, survey = "COVIMOD", new_hh=TRUE)

# # Add missing u index
# dt.cnt[, u := fcase(wave == 1, 1)]

# Add flattened list of ages of contacts and corresponding cumulative list
stan_data <- add_ages_contacts(stan_data, dt.offsets)

# Add row major index
stan_data <- add_row_major_idx(stan_data, dt.cnt, survey="simulation")

# Add household offsets
stan_data <- add_household_offsets(stan_data, dt.offsets, no_log=FALSE, genius=TRUE)

# Map age to age strata
stan_data <- add_map_age_to_strata(stan_data, survey="simulation")

# Map individual to age for each gender combination
stan_data <- add_map_indiv_to_age(stan_data, dt.cnt, dt.offsets)

# Add Non-nuisance index
stan_data <- add_nn_idx(stan_data)

# Add standardized age indexes
stan_data <- add_std_age_idx(stan_data)

# Add HSGP parameters
stan_data <- add_hsgp_parms(stan_data, args$hsgp_c, args$hsgp_m , args$hsgp_m)
cat(" DONE!\n")

#### ---------- Assess convergence and mixing ---------- #####
if(args$mixing){
  cat(" Assess convergence and mixing\n")
  
  # Make convergence diagnostic tables
  fit_summary <- make_convergence_diagnostic_stats(fit, outdir=export.path)
  
  # Make trace plots
  cat("\n Making trace plots")
  bayesplot::color_scheme_set(scheme = "mix-blue-pink")
  
  pars <- c('gp_alpha', 'gp_rho_1', 'gp_rho_2')
  
  pars_po <- fit$draws(pars)
  p <- bayesplot::mcmc_trace(pars_po)
  ggsave(file = file.path(export.fig.path, 'mcmc_trace_parameters.png'), plot = p, h = 20, w = 20, limitsize = F)
  
  # Make pairs plots
  cat(" Making pairs plots\n")
  p <- bayesplot::mcmc_pairs(pars_po, off_diag_args=list(size=0.3, alpha=0.3))
  ggsave(file = file.path(export.fig.path, 'mcmc_pairs_parameters.png'), plot = p, h = 20, w = 20, limitsize = F)
  
  cat("\n DONE!\n")
}

##### ---------- Posterior predictive checks ---------- #####
if(args$ppc){
  cat(" Extracting posterior\n")
  po <- fit$draws(c("yhat_strata_MM", "yhat_strata_FF", "yhat_strata_MF", "yhat_strata_FM"), 
                  inc_warmup = FALSE, format="draws_matrix")
  
  cat(" Making posterior predictive checks\n")
  dt.ppc <- make_ppd_check_covimod(po, dt.offsets, stan_data, outdir=export.path, fig.outdir=export.fig.path, new_hh=TRUE)
  
  cat("\n DONE.\n")
}

##### ---------- Error Table ---------- #####

dt.ppc[, cntct_intensity_predict := M]
dt.ppc[, cntct_intensity := y]
dt.ppc[, cntct_rate:=y/Hic_b]
dt.ppc[, cntct_rate_predict:=M/Hic_b]
dt.ppc <- unique(dt.ppc, by=c("age", "gender", "alter_gender", "alter_age_strata", "new_id"))
error_table <- make_error_table(dt.ppc, rate=TRUE)

saveRDS(dt.ppc, file.path(outdir=export.path, "error_dt.rds"))
saveRDS(error_table, file.path(outdir=export.path, "error_table.rds"))

p3 <- ggplot(dt.ppc) +
  geom_tile(aes(x=age, y=alter_age_strata, fill = cntct_intensity)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  coord_equal() +
  viridis::scale_fill_viridis(na.value = "white", option="H") +
  labs(x="Participants' age", y="Contacts' age ", fill="Counts") +
  facet_grid( paste(alter_gender, "(Contacts)") ~ paste(gender, "(Participants)") ) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "bottom", text = element_text(size = 5),
        strip.background = element_rect(color=NA, fill = "transparent"))

p4 <- ggplot(dt.ppc) +
  geom_tile(aes(x=age, y=alter_age_strata, fill = cntct_intensity_predict)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  coord_equal() +
  viridis::scale_fill_viridis(na.value = "white", option="H") +
  labs(x="Participants' age", y="Contacts' age ", fill="Predicted counts") +
  facet_grid( paste(alter_gender, "(Contacts)") ~ paste(gender, "(Participants)") ) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "bottom", text = element_text(size = 5),
        strip.background = element_rect(color=NA, fill = "transparent"))

ggsave(file.path(export.fig.path, "empirical_cntct_intensity.png"), plot = p3, height = 3, width = 7)
ggsave(file.path(export.fig.path, "estimated_cntct_intensity.png"), plot = p4, height = 3, width = 7)

cat("\n DONE.\n")

##### ---------- Plotting ---------- #####
if(args$plot){
  cat(" Extracting posterior contact intensities\n")
  po <- fit$draws(c("log_cnt_rate"), inc_warmup = FALSE, format="draws_matrix")
  dt.po <- extract_posterior_rates(po)
  dt.matrix <- posterior_contact_intensity(dt.po, dt.pop, dt.sim.true.cntct=dt.sim.true.cntct, type="matrix", outdir=export.path, new_hh=TRUE, sim=TRUE)
  dt.margin <- posterior_contact_intensity(dt.po, dt.pop, type="marginal", outdir=export.path, new_hh=TRUE)
  
  rm(dt.po); suppressMessages(gc()); # Ease memory
  rm(dt.ppc); suppressMessages(gc()); # Ease memory
  
  po.alphaMM <- fit$draws(c("alpha_age_MM"), inc_warmup = FALSE, format="draws_matrix")
  po.alphaFF <- fit$draws(c("alpha_age_FF"), inc_warmup = FALSE, format="draws_matrix")
  po.alphaMF <- fit$draws(c("alpha_age_MF"), inc_warmup = FALSE, format="draws_matrix")
  po.alphaFM <- fit$draws(c("alpha_age_FM"), inc_warmup = FALSE, format="draws_matrix")
  
  dt.po.alphaMM <- extract_posterior_alpha(po.alphaMM, gender_comb = "MM")
  dt.matrix.alphaMM <- posterior_alpha(fit, dt.po.alphaMM, stan_data, type="matrix", outdir=export.path, gender_comb="MM")
  dt.margin.alphaMM <- posterior_alpha(fit, dt.po.alphaMM, stan_data, type="marginal", outdir=export.path, gender_comb="MM")
  
  
  dt.po.alphaFF <- extract_posterior_alpha(po.alphaFF, gender_comb = "FF")
  dt.matrix.alphaFF <- posterior_alpha(fit, dt.po.alphaFF, stan_data, type="matrix", outdir=export.path, gender_comb="FF")
  dt.margin.alphaFF <- posterior_alpha(fit, dt.po.alphaFF, stan_data, type="marginal", outdir=export.path, gender_comb="FF")
  
  
  dt.po.alphaMF <- extract_posterior_alpha(po.alphaMF, gender_comb = "MF")
  dt.matrix.alphaMF <- posterior_alpha(fit, dt.po.alphaMF, stan_data, type="matrix", outdir=export.path, gender_comb="MF")
  dt.margin.alphaMF <- posterior_alpha(fit, dt.po.alphaMF, stan_data, type="marginal", outdir=export.path, gender_comb="MF")
  
  
  dt.po.alphaFM <- extract_posterior_alpha(po.alphaFM, gender_comb = "FM")
  dt.matrix.alphaFM <- posterior_alpha(fit, dt.po.alphaFM, stan_data, type="matrix", outdir=export.path, gender_comb="FM")
  dt.margin.alphaFM <- posterior_alpha(fit, dt.po.alphaFM, stan_data, type="marginal", outdir=export.path, gender_comb="FM")
  
  # combine datasets with all gender combinations
  dt.matrix.alpha <- rbind(dt.matrix.alphaMM, dt.matrix.alphaFF, dt.matrix.alphaMF, dt.matrix.alphaFM)
  dt.margin.alpha <- rbind(dt.margin.alphaMM, dt.margin.alphaFF, dt.margin.alphaMF, dt.margin.alphaFM)
  
  saveRDS(dt.matrix.alpha, file.path(outdir=export.path, "alpha_matrix.rds"))
  saveRDS(dt.margin.alpha, file.path(outdir=export.path, "alpha_margin.rds"))
  
  # Ease memory
  rm(dt.po.alphaMM); suppressMessages(gc());
  rm(dt.matrix.alphaMM); suppressMessages(gc());
  rm(dt.po.alphaFF); suppressMessages(gc());
  rm(dt.matrix.alphaFF); suppressMessages(gc());
  rm(dt.po.alphaMF); suppressMessages(gc());
  rm(dt.matrix.alphaMF); suppressMessages(gc());
  rm(dt.po.alphaFM); suppressMessages(gc());
  rm(dt.matrix.alphaFM); suppressMessages(gc());
  
  cat(" Making figures\n")
  
  p <- plot_posterior_intensities(dt.matrix, outdir=export.path, new_hh=TRUE)
  p <- plot_sliced_intensities(dt.matrix.alpha, outdir=export.path, new_hh=TRUE)
  p <- plot_sliced_intensities(dt.matrix, outdir=export.path, new_hh_intensity=TRUE)
  p <- plot_marginal_intensities(dt.margin, outdir=export.path, new_hh=TRUE, rate=TRUE)
  p <- plot_marginal_intensities(dt.margin.alpha, outdir=export.path, new_hh=TRUE)
  p <- plot_alpha(dt.matrix.alpha, outdir=export.fig.path)
  
  cat("\n DONE.\n")
}

cat("\n ########## Second Error Table ############")

group_var_0 <- c("part_idx", "age", "gender", "alter_gender", "alter_age_strata")
dt.matrix.alpha[, alter_age_strata := fcase(
  alter_age <= 4,  "0-4",
  alter_age <= 9,  "5-9",
  alter_age <= 14, "10-14",
  alter_age <= 19, "15-19",
  alter_age <= 24, "20-24",
  alter_age <= 34, "25-34",
  alter_age <= 44, "35-44",
  alter_age <= 54, "45-54",
  alter_age <= 64, "55-64",
  alter_age <= 69, "65-69",
  alter_age <= 74, "70-74",
  alter_age <= 79, "75-79",
  alter_age <= 84, "80-84",
  alter_age > 84, "85+",
  default = NA
)]


  dt.matrix.alpha[, cntct_intensity_predict := sum(M), by=c("part_idx", "age", "gender", "alter_gender", "alter_age_strata")]
  
  # ordering dt.cnt like in model
  covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
  # making sure order of factors in alter_age_strata is ascending instead of decreasing
  # note alter_age_strata_idx and age_strata_idx are the same, but they serve different purposes
  dt.cnt[, alter_age_strata_idx:=as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
  dt.cnt<- dt.cnt[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
  # check dt.cnt[n==0,]
  dt_intensity_MM <- dt.cnt[gender=="Male"& alter_gender=="Male"]
  dt_intensity_FF <- dt.cnt[gender=="Female"& alter_gender=="Female"]
  dt_intensity_MF <- dt.cnt[gender=="Male"& alter_gender=="Female"]
  dt_intensity_FM <- dt.cnt[gender=="Female"& alter_gender=="Male"]
  
  # adding a part_idx
  dt_intensity_MM[ , part_idx := .GRP, by = new_id]      
  dt_intensity_FF[ , part_idx := .GRP, by = new_id]      
  dt_intensity_MF[ , part_idx := .GRP, by = new_id]  
  dt_intensity_FM[ , part_idx := .GRP, by = new_id]      
  
  # check max(dt.matrix.alpha[gender=="Male" & alter_gender=="Male"]$part_idx) == max(dt_intensity_MM$part_idx)
  
  # bind everything
  dt_intensity <- rbind(dt_intensity_MM, dt_intensity_FF, dt_intensity_MF, dt_intensity_FM)
  setnames(dt_intensity, c("y"), c("cntct_intensity"))
  dt_intensity <- merge(dt_intensity, dt.matrix.alpha, by=c("part_idx", "age", "gender", "alter_age_strata", "alter_gender"), all.x = TRUE)
  dt_intensity <- unique(dt_intensity, by = group_var_0)
  


error_table_2 <- make_error_table(dt_intensity)

saveRDS(dt_intensity, file.path(outdir=export.path, "dt_intensity_2.rds"))
saveRDS(error_table_2, file.path(outdir=export.path, "error_table_2.rds"))

p5 <- ggplot(dt_intensity) +
  geom_tile(aes(x=age, y=alter_age_strata, fill = cntct_intensity)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  coord_equal() +
  viridis::scale_fill_viridis(na.value = "white", option="H") +
  labs(x="Participants' age", y="Contacts' age ", fill="Counts") +
  facet_grid( paste(alter_gender, "(Contacts)") ~ paste(gender, "(Participants)") ) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "bottom", text = element_text(size = 5),
        strip.background = element_rect(color=NA, fill = "transparent"))

p6 <- ggplot(dt_intensity) +
  geom_tile(aes(x=age, y=alter_age_strata, fill = cntct_intensity_predict)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  coord_equal() +
  viridis::scale_fill_viridis(na.value = "white", option="H") +
  labs(x="Participants' age", y="Contacts' age ", fill="Estimated contact intensity") +
  facet_grid( paste(alter_gender, "(Contacts)") ~ paste(gender, "(Participants)") ) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "bottom", text = element_text(size = 5),
        strip.background = element_rect(color=NA, fill = "transparent"))


ggsave(file.path(export.fig.path, "empirical_cntct_intensity_2.png"), plot = p5, height = 3, width = 7)
ggsave(file.path(export.fig.path, "estimated_cntct_intensity_2.png"), plot = p6, height = 3, width = 7)

cat("\n DONE.\n")
