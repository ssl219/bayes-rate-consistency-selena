# Creates simulated datasets from simulated intensity data
# For detailed documentation on how this works refer to the Simulating contact dataset notebook

cat("\n---------- Simulating contact dataset ----------\n")

# Load libraries
library(optparse)
library(data.table)
library(ggplot2)
library(ggpubr)

option_list <- list(
  optparse::make_option("--seed", type = "integer", default = 0721, dest = "seed"),
  optparse::make_option("--size", type = "integer", default = 85,
                        help = "Number of participants with random age in the survey [default \"%default\"]",
                        dest = 'size'),
  optparse::make_option("--nsim", type = "integer", default = 10,
                        help = "Number of simulated datasets [default \"%default\"]",
                        dest = "nsim"),
  optparse::make_option("--strata", type = "character", default = "COVIMOD",
                        help = "Age stratification scheme [default %default]",
                        dest = "strata"),
  optparse::make_option("--hhsize", type = "character", default = 4,
                        help = "Household size [default %default]",
                        dest = "hhsize"),
  optparse::make_option("--scenario", type = "character", default = "flat",
                        help = "Scenario [default %default]",
                        dest = "scenario"),
  optparse::make_option("--repo_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate/bayes-rate-consistency-selena",
                        help = "Absolute file path to repository directory",
                        dest = 'repo.path'),
  optparse::make_option("--data_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate",
                        help = "Absolute file path to data directory, used as long we don t build an R package [default]",
                        dest = 'data.path')
)
# "/rds/general/user/ssl219/home/bayes-rate-consistency-selena"
args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# args$repo.path <- "~/Imperial/covimod-gp"
# args$strata <- "5yr"

##### ---------- Error handling ---------- #####
if(is.na(args$repo.path)){
  stop("Please specify --repo_path")
}

if(is.na(args$scenario)){
  stop("Please specify --scenario")
}

###### ---------- Load data ---------- #####
source(file.path(args$repo.path, "R", "sim-dataset-utility.R"))


if (args$scenario == "flat"){
  dt <- as.data.table(readRDS( file.path(args$data.path, "flat-data.rds") ))
}

if (args$scenario == "boarding_school"){
  dt <- as.data.table(readRDS( file.path(args$data.path, "boarding-data.rds") ))
}


##### ---------- Data export ---------- #####
dir.name <- paste("new-hh", args$scenario, sep="-")
export.path <- file.path(args$data.path, "data/simulations/datasets", dir.name)
if(!dir.exists(export.path)){
  cat(paste("\n Making export directory:", export.path))
  dir.create(export.path)
}

##### ---------- Stratify age groups ---------- #####
cat("\n Stratifying age groups ...")
dt <- stratify_alter_age(dt, args$strata)

##### ---------- Simulating contact survey data ---------- ##########
cat("\n Generating contact dataset ...")
args$size = 85
args$strata = "COVIMOD-new-hh"
N_random <- args$size 
N = 85 + N_random
# Generate new_id for all participants (all distinct)
d.everything <- as.data.table(expand.grid(new_id = 1:N, alter_age = 0:84))
d.everything <- stratify_alter_age(d.everything, args$strata)
d.everything[, wave:=1L]
# hhsize = args$hhsize
d.everything[, household_size:=4L]

# Generate age of participants
set.seed(2002)
random_ages <- sample(0:84, N_random, replace=TRUE)
part_ages <- sort(c(random_ages, 0:84))
# create smaller dataset with new_id and part_ages
id_ages <- data.table(new_id=1:N, age=part_ages)
d.everything <- merge(d.everything, id_ages, by=c("new_id"), all=TRUE)

# Calculate average contact counts
dt <- merge(dt, dpop[, list(age, gender, part)], by=c("age", "gender"))
dt[, mu := round(cntct_intensity*part)]

# Simulate contact counts in survey
for (i in 1:args$nsim){
  set.seed(args$seed + i)

  dt[, y := rpois(nrow(dt), lambda=dt$mu)]

  # Stratify contact intensities and contact rates
  group_var <- c("age", "gender", "alter_age_strata", "alter_gender")
  dt[, y_strata := sum(y), by=group_var]
  dt[, cntct_intensity_strata := sum(cntct_intensity), by=group_var]
  dt[, pop_strata := sum(pop), by=group_var]
  dt[, cntct_rate_strata := sum(cntct_intensity)/sum(pop), by=group_var]

  # Save data
  saveRDS(dt, file = file.path(export.path, paste0("data","_",i,".rds")))
}

cat("\n DONE. \n")
