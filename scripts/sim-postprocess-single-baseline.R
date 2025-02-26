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
  optparse::make_option("--data_path", type = "character", default = "/rds/general/user/ssl219/home",
                        help = "Absolute file path to data directory, used as long we don t build an R package [default]",
                        dest = 'data.path'),
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
                        dest = "sim.no"),
  optparse::make_option("--hhsize", type = "integer", default = 4,
                        help = "Household size [default %default]",
                        dest = "hhsize"),
  optparse::make_option("--sample_size", type = "integer", default = 55,
                        help = "Number of participants with random age in the survey [default \"%default\"]",
                        dest = 'sample_size'),
  optparse::make_option("--scenario", type = "character", default = "flat",
                        help = "Scenario [default %default]",
                        dest = "scenario")
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
# args$repo.path <- "/Users/mac/Documents/M4R/code/bayes_consistency_rate/bayes-rate-consistency-selena"
# args$data.path <- "/Users/mac/Documents/M4R/code/bayes_consistency_rate"
# args$model.name <- "hsgp-eq-cd-1-nhh"

model.path <- file.path(args$repo.path, "stan_fits", paste0(args$model.name, ".rds"))
full.data.path <- file.path(args$data.path, "data/simulations/datasets/new-hh-both", paste0("hh", args$hhsize, "-", args$sample_size), paste0("dataset", args$sim.no), paste0("data-hh", args$hhsize, "-", args$scenario, "-", args$sample_size, "-amended-baseline.rds"))
# Error handling
if(!file.exists(model.path)) {
  cat("\n Model: ", model.path)
  stop("The specified model does not exist!")
}
if(!file.exists(full.data.path)) {
  stop("The specified dataset does not exists!")
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
cat(paste("\n Model path:", model.path))
cat(paste("\n Full data path:", full.data.path))

fit <- readRDS(model.path)
data <- readRDS(full.data.path)
dt.cnt <- data$contacts[wave == args$wave]
dt.offsets <- data$offsets[wave == args$wave]
dt.pop <- data$pop
setnames(dt.pop, "alter_age", "age")

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





##### ---------- Assess convergence and mixing ---------- #####
if(args$mixing){
  cat(" Assess convergence and mixing\n")
  
  # Make convergence diagnostic tables
  fit_summary <- make_convergence_diagnostic_stats(fit, outdir=export.path)
  
  # Make trace plots
  cat("\n Making trace plots")
  bayesplot::color_scheme_set(scheme = "mix-blue-pink")
  
  pars <- c('nu', 'gp_alpha', 'gp_rho_1', 'gp_rho_2')
  
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
  po <- fit$draws(c("yhat_strata", "log_cnt_rate"), inc_warmup = FALSE, format="draws_matrix")
  
  cat(" Making posterior predictive checks\n")
  dt.ppc <- make_ppd_check_covimod(po, dt.cnt, outdir=export.path, fig.outdir = export.fig.path)
  
  cat("\n DONE.\n")
}

##### ---------- Error Table ---------- #####

dt.ppc[, cntct_count_predict := M]
dt.ppc[, cntct_count := y]
dt.ppc[, cntct_intensity_predict := M/N]
dt.ppc[, cntct_intensity := y/N]
error_table <- make_error_table(dt.ppc, count=TRUE)
saveRDS(dt.ppc, file.path(outdir=export.path, "error_dt.rds"))
saveRDS(error_table, file.path(outdir=export.path, "error_table.rds"))
cat("\n DONE.\n")

##### ---------- Plotting ---------- #####
if(args$plot){
  cat(" Extracting posterior contact intensities\n")
  dt.po <- extract_posterior_rates(po)
  dt.matrix <- posterior_contact_intensity(dt.po, dt.pop, dt.sim.true.cntct=dt.sim.true.cntct, type="matrix", outdir=export.path, sim=TRUE)
  dt.margin <- posterior_contact_intensity(dt.po, dt.pop, type="marginal", outdir=export.path)
  
  rm(dt.po); suppressMessages(gc()); # Ease memory
  
  cat(" Making figures\n")
  
  p <- plot_posterior_intensities(dt.matrix, outdir=export.path)
  p <- plot_sliced_intensities(dt.matrix, outdir=export.path)
  p <- plot_marginal_intensities(dt.margin, outdir=export.path)
  
  cat("\n DONE.\n")
}

