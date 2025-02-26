# Run Stan models on COVIMOD data
cat("\n---------- Compiling and Running Stan Model ----------\n")

# Load libraries
library(optparse)
library(data.table)
library(stringr)
library(cmdstanr)
library(tidyr)
library(dplyr)

##### ---------- I/O ---------- #####
option_list <- list(
  optparse::make_option("--seed", type = "integer", default = 0721,
                        help = "Random number seed [default %default]",
                        dest = "seed"),
  optparse::make_option("--iter_warmup", type = "integer", default = 1,
                        help = "HMC warmup iterations [default %default]",
                        dest = 'iter.warmup'),
  optparse::make_option("--iter_sampling", type = "integer", default = 2,
                        help = "HMC of sampling iterations iterations [default %default]",
                        dest = 'iter.sampling'),
  optparse::make_option("--chains", type = "integer", default = 1,
                        help = "Number of MCMC chains",
                        dest = 'chains'),
  optparse::make_option("--model", type = "character", default = "hsgp-eq-rd-new-hh-dropping-all-zeros-symmetric-poisson",
                        help = "Name of Stan model",
                        dest = 'model.name'),
  optparse::make_option("--hsgp_c", type = "double", default = 1.5,
                        help = "The boundary inflation of the HSGP prior in any dimension [default \"%default\"]",
                        dest = "hsgp_c"),
  optparse::make_option("--hsgp_m", type = "integer", default = 20,
                        help = "The number of the HSGP basis functions in any dimension [default \"%default\"]",
                        dest = "hsgp_m"),
  # optparse::make_option("--repo_path", type = "character", default = "/rds/general/user/ssl219/home/bayes-rate-consistency-selena",
  #                       help = "Absolute file path to repository directory, used as long we don t build an R package [default]",
  #                       dest = 'repo.path'),
  # optparse::make_option("--data_path", type = "character", default = "/rds/general/user/ssl219/home",
  #                       help = "Absolute file path to data directory, used as long we don t build an R package [default]",
  #                       dest = 'data.path'),
  optparse::make_option("--repo_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate/bayes-rate-consistency-selena",
                        help = "Absolute file path to repository directory, used as long we don t build an R package [default]",
                        dest = 'repo.path'),
  optparse::make_option("--data_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate",
                        help = "Absolute file path to data directory, used as long we don t build an R package [default]",
                        dest = 'data.path'),
  optparse::make_option("--wave", type = "integer", default = 1,
                        help = "COVIMOD wave",
                        dest = "wave"),
  optparse::make_option("--mixing", type = "logical", default = TRUE,
                        help = "Whether to assess mixing",
                        dest = "mixing"),
  optparse::make_option("--ppc", type = "logical", default = TRUE,
                        help = "Whether to run posterior predictive checks",
                        dest = "ppc"),
  optparse::make_option("--plot", type = "logical", default = TRUE,
                        help = "Whether to plot posterior distributions",
                        dest = "plot")
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Load helpers
source(file.path(args$repo.path, "R/stan-utility.R"))

# Load data
covimod.single.new.hh <- readRDS(file.path(args$data.path, "data/COVIMOD/COVIMOD-single-new-hh.rds"))

dt.cnt <- covimod.single.new.hh$contacts[wave == args$wave]
dt.offsets <- covimod.single.new.hh$offsets[wave == args$wave]
dt.pop <- covimod.single.new.hh$pop

# Path to model
model.path <- paste0(file.path(args$repo.path, "stan_models", args$model.name), ".stan")

# Export path
export.path <- file.path(args$repo.path, "stan_fits")
if (!file.exists(export.path)) {
  cat(paste(" Making export directory:", export.path, "\n"))
  dir.create(export.path)
}

## Configure Stan data
cat(" Configuring Stan data ...")

# Initialize
stan_data <- init_stan_data()

# Add contact counts
stan_data <- add_contact_vector(stan_data, dt.cnt, survey="COVIMOD", single = TRUE, new_hh = TRUE)

# Add obs counts
stan_data <- add_N(stan_data, dt.cnt, survey = "COVIMOD", new_hh=TRUE)

# # Add missing u index
# dt.cnt[, u := fcase(wave == 1, 1)]

# Add flattened list of ages of contacts and corresponding cumulative list
stan_data <- add_ages_contacts(stan_data, dt.offsets)

# Add row major index
stan_data <- add_row_major_idx(stan_data, dt.cnt, survey="POLYMOD_2")

# Add household offsets
stan_data <- add_household_offsets(stan_data, dt.offsets)

# Map age to age strata
stan_data <- add_map_age_to_strata(stan_data)

# Map individual to age for each gender combination
stan_data <- add_map_indiv_to_age(stan_data, dt.cnt, dt.offsets)

# Add Non-nuisance index
stan_data <- add_nn_idx(stan_data)

# Add standardized age indexes
stan_data <- add_std_age_idx(stan_data)

# Add HSGP parameters
stan_data <- add_hsgp_parms(stan_data, args$hsgp_c, args$hsgp_m , args$hsgp_m)
cat(" DONE!\n")

cat(" Compiling Stan model ...")
# Compile stan program
model <- cmdstanr::cmdstan_model(model.path, force_recompile = TRUE)
cat(" DONE!\n")

cat(" Running Stan model ...\n")
fit <- model$sample(
  data = stan_data,
  chains = args$chains,
  seed = args$seed,
  refresh = 100,
  iter_warmup = args$iter.warmup,
  iter_sampling = args$iter.sampling,
  parallel_chains = args$chains,
  max_treedepth = 13
)
cat(" DONE!\n")

cat("\n Run Stan ALL DONE.\n")


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



# optparse::make_option("--repo_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate/bayes-rate-consistency-selena",
#                        help = "Absolute file path to repository directory, used as long we don t build an R package [default]",
#                        dest = "repo.path"),
# optparse::make_option("--data_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate",
#                        help = "Absolute file path to data directory, used as long we don t build an R package [default]",
#                        dest = 'data.path'),

data.path <- file.path(args$data.path, "data/COVIMOD/COVIMOD-single-new-hh.rds")

# Error handling
if(!file.exists(data.path)) {
  stop("The specified dataset does not exist!")
}

# Output directories
export.path <- file.path(args$repo.path, "results", args$model.name)
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
cat(paste("\n Data path:", data.path))

data <- readRDS(data.path)
dt.cnt <- data$contacts[wave == args$wave]
dt.offsets <- data$offsets[wave == args$wave]
dt.pop <- data$pop

source(file.path(args$repo.path, "R/covimod-utility.R"))
source(file.path(args$repo.path, "R/stan-utility.R"))
source(file.path(args$repo.path, "R/postprocess-diagnostic-single.R"))
source(file.path(args$repo.path, "R/postprocess-plotting-single.R"))

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
  make_ppd_check_covimod(po, dt.offsets, stan_data, outdir=export.path, fig.outdir=export.fig.path, new_hh=TRUE)
  
  cat("\n DONE.\n")
}

##### ---------- Plotting ---------- #####
if(args$plot){
  cat(" Extracting posterior contact intensities\n")
  po <- fit$draws(c("log_cnt_rate"), inc_warmup = FALSE, format="draws_matrix")
  dt.po <- extract_posterior_rates(po)
  dt.matrix <- posterior_contact_intensity(dt.po, dt.pop, type="matrix", outdir=export.path, new_hh=TRUE)
  # dt.margin <- posterior_contact_intensity(dt.po, dt.pop, type="marginal", outdir=export.path, new_hh=TRUE)
  
  rm(dt.po); suppressMessages(gc()); # Ease memory
  
  po.alphaMM <- fit$draws(c("alpha_age_MM"), inc_warmup = FALSE, format="draws_matrix")
  po.alphaFF <- fit$draws(c("alpha_age_FF"), inc_warmup = FALSE, format="draws_matrix")
  po.alphaMF <- fit$draws(c("alpha_age_MF"), inc_warmup = FALSE, format="draws_matrix")
  po.alphaFM <- fit$draws(c("alpha_age_FM"), inc_warmup = FALSE, format="draws_matrix")
  
  dt.po.alphaMM <- extract_posterior_alpha(po.alphaMM, gender_comb = "MM")
  dt.matrix.alphaMM <- posterior_alpha(fit, dt.po.alphaMM, stan_data, type="matrix", outdir=export.path, gender_comb="MM")
  
  
  dt.po.alphaFF <- extract_posterior_alpha(po.alphaFF, gender_comb = "FF")
  dt.matrix.alphaFF <- posterior_alpha(fit, dt.po.alphaFF, stan_data, type="matrix", outdir=export.path, gender_comb="FF")
  
  
  dt.po.alphaMF <- extract_posterior_alpha(po.alphaMF, gender_comb = "MF")
  dt.matrix.alphaMF <- posterior_alpha(fit, dt.po.alphaMF, stan_data, type="matrix", outdir=export.path, gender_comb="MF")
  
  
  dt.po.alphaFM <- extract_posterior_alpha(po.alphaFM, gender_comb = "FM")
  dt.matrix.alphaFM <- posterior_alpha(fit, dt.po.alphaFM, stan_data, type="matrix", outdir=export.path, gender_comb="FM")
  
  # combine datasets with all gender combinations
  dt.matrix.alpha <- rbind(dt.matrix.alphaMM, dt.matrix.alphaFF, dt.matrix.alphaMF, dt.matrix.alphaFM)
  
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
  # p <- plot_sliced_intensities(dt.matrix, outdir=export.path)
  # p <- plot_marginal_intensities(dt.margin, outdir=export.path)
  
  p <- plot_alpha(dt.matrix.alpha, outdir=export.fig.path)
}

# save at the end

cat(" Saving fitted model ...")
args$model.name <- paste(args$model.name, args$wave, sep="-")
fit$save_object(file = file.path(export.path, paste0(args$model.name, ".rds")))
cat(" DONE!\n")

