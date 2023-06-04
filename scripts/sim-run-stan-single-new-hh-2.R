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
  optparse::make_option("--iter_warmup", type = "integer", default = 500,
                        help = "HMC warmup iterations [default %default]",
                        dest = 'iter.warmup'),
  optparse::make_option("--iter_sampling", type = "integer", default = 1000,
                        help = "HMC of sampling iterations iterations [default %default]",
                        dest = 'iter.sampling'),
  optparse::make_option("--chains", type = "integer", default = 4,
                        help = "Number of MCMC chains",
                        dest = 'chains'),
  optparse::make_option("--hhsize", type = "integer", default = 4,
                        help = "Household size [default %default]",
                        dest = "hhsize"),
  optparse::make_option("--sample_size", type = "integer", default = 55,
                        help = "Number of participants with random age in the survey [default \"%default\"]",
                        dest = 'sample_size'),
  optparse::make_option("--scenario", type = "character", default = "flat",
                        help = "Scenario [default %default]",
                        dest = "scenario"),
  optparse::make_option("--model", type = "character", default = "hsgp-eq-rd-new-hh-dropping-all-zeros-symmetric-poisson",
                        help = "Name of Stan model",
                        dest = 'model.name'),
  optparse::make_option("--repo_path", type = "character", default = "/rds/general/user/ssl219/home/bayes-rate-consistency-selena",
                        help = "Absolute file path to repository directory, used as long we don t build an R package [default]",
                        dest = "repo.path"),
  optparse::make_option("--data_path", type = "character", default = "/rds/general/user/ssl219/home",
                        help = "Absolute file path to data directory, used as long we don t build an R package [default]",
                        dest = 'data.path'),
  optparse::make_option("--hsgp_c", type = "double", default = 1.5,
                        help = "The boundary inflation of the HSGP prior in any dimension [default \"%default\"]",
                        dest = "hsgp_c"),
  optparse::make_option("--hsgp_m", type = "integer", default = 30,
                        help = "The number of the HSGP basis functions in any dimension [default \"%default\"]",
                        dest = "hsgp_m"),
  optparse::make_option("--wave", type = "integer", default = 1,
                        help = "COVIMOD wave",
                        dest = "wave"),
  optparse::make_option("--sim.no", type = "integer", default = 1,
                        help = "Simulated Dataset Number [default %default]",
                        dest = "sim.no")
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
# args$repo.path = "/Users/mac/Documents/M4R/code/bayes_consistency_rate/bayes-rate-consistency-selena"
# args$data.path = "/Users/mac/Documents/M4R/code/bayes_consistency_rate"

# Load helpers
source(file.path(args$repo.path, "R/stan-utility.R"))

# Load data
data.path.test <- file.path(args$data.path, "data/simulations/datasets/new-hh-both", paste0("hh", args$hhsize, "-", args$sample_size), paste0("dataset", args$sim.no), paste0("data-hh", args$hhsize, "-", args$scenario, "-", args$sample_size, "-amended.rds"))
cat ("\n DATA PATH RUN", data.path.test)
covimod.single.new.hh <- readRDS(data.path.test)

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

cat(" Saving fitted model ...")
args$model.name <- paste(args$model.name, args$wave, sep="-")
fit$save_object(file = file.path(export.path, paste0(args$model.name, "-sim-hh", args$hhsize, "-", args$scenario, "-", args$sample_size, "-", args$sim.no, ".rds")))
cat(" DONE!\n")

cat("\n Run Stan ALL DONE.\n")
