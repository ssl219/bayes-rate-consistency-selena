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
  optparse::make_option("--model", type = "character", default = "hsgp-eq-rd-new-hh-dropping-all-zeros-symmetric-poisson-1-random",
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
                        dest = "plot")
)


# optparse::make_option("--repo_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate/bayes-rate-consistency-selena",
#                        help = "Absolute file path to repository directory, used as long we don t build an R package [default]",
#                        dest = "repo.path"),
# optparse::make_option("--data_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate",
#                        help = "Absolute file path to data directory, used as long we don t build an R package [default]",
#                        dest = 'data.path'),
args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

model.path <- file.path(args$repo.path, "stan_fits", paste0(args$model.name, ".rds"))
data.path <- file.path(args$data.path, "data/COVIMOD/COVIMOD-single-new-hh-random.rds")

# Error handling
if(!file.exists(model.path)) {
  cat("\n Model: ", model.path)
  stop("The specified model does not exist!")
}
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


#### ------------------- Stan data --------------------- #####

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
  dt.matrix <- posterior_contact_intensity(dt.po, dt.pop, type="matrix", outdir=export.path, new_hh=TRUE)
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