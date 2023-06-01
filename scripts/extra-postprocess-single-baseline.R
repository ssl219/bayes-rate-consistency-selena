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
                        dest = "plot")
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

args$model.name <- "hsgp-eq-rd-1-hh"
args$repo.path <- "/Users/mac/Documents/M4R/code/bayes_consistency_rate/bayes-rate-consistency-selena"
args$data.path <- "/Users/mac/Documents/M4R/code/bayes_consistency_rate"
intensity.path <- file.path("/Users/mac/Documents/M4R/hpc_results", paste0(args$model.name), "intensity_matrix.rds")
full.data.path <- file.path(args$data.path, "data/COVIMOD/COVIMOD-single-hh.rds")

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

intensity_matrix <- readRDS(intensity.path)
data <- readRDS(full.data.path)
dt.cnt <- data$contacts[wave == args$wave]
dt.offsets <- data$offsets[wave == args$wave]
dt.pop <- data$pop

source(file.path(args$repo.path, "R/covimod-utility.R"))
source(file.path(args$repo.path, "R/stan-utility.R"))
source(file.path(args$repo.path, "R/postprocess-diagnostic-single.R"))
source(file.path(args$repo.path, "R/postprocess-plotting-single.R"))

intensity_matrix[, alter_age_strata := fcase(
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

intensity_matrix <- intensity_matrix[,strata_cntct_intensity := sum(intensity_M), by = c("age", "gender", "alter_age_strata", "alter_gender") ]
intensity_matrix <- unique(intensity_matrix, by = c("age", "gender", "alter_age_strata", "alter_gender"))

dt_intensity <- merge(dt.cnt, intensity_matrix, by = c("age", "gender", "alter_age_strata", "alter_gender"), all.x = TRUE, all.y=FALSE)
dt_intensity <- unique(dt_intensity, by = c("age", "gender", "alter_age_strata", "alter_gender"))
dt_intensity[, cntct_intensity:=y/N]
dt_intensity[, cntct_intensity_predict:=strata_cntct_intensity]

error_table <- make_error_table(dt_intensity)
error_table
