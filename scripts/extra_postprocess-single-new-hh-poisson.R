# Preamble: Generates diagnostic statistics and result plots

cat("\n ---------- Begin Extra Post-processing ---------- \n")

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
library(data.table)   

##### ---------- I/O ---------- #####
option_list <- list(
  # optparse::make_option("--repo_path", type = "character", default = "/rds/general/user/ssl219/home/bayes-rate-consistency-selena",
  #                      help = "Absolute file path to repository directory, used as long we don t build an R package [default]",
  #                      dest = "repo.path"),
  # optparse::make_option("--data_path", type = "character", default = "/rds/general/user/ssl219/home",
  #                      help = "Absolute file path to data directory, used as long we don t build an R package [default]",
  #                      dest = 'data.path'),
  optparse::make_option("--repo_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate/bayes-rate-consistency-selena",
                        help = "Absolute file path to repository directory, used as long we don t build an R package [default]",
                        dest = "repo.path"),
  optparse::make_option("--data_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate",
                        help = "Absolute file path to data directory, used as long we don t build an R package [default]",
                        dest = 'data.path'),
  optparse::make_option("--wave", type="integer", default = 1,
                        help = "COVIMOD wave",
                        dest = "wave"),
  optparse::make_option("--model", type = "character", default = "hsgp-eq-cd-new-hh-dropping-all-zeros-symmetric-poisson-1-sim-flat-everyone-ppd",
                        #"hsgp-eq-rd-new-hh-dropping-all-zeros-symmetric-poisson",
                        help = "Name of the model",
                        dest = "model.name")
)

cat("\n before args")

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
cat("\n model name:", args$model.name)

cat("\n after args")

model.path <- file.path(args$repo.path, "stan_fits", paste0(args$model.name, ".rds"))
data.path <- file.path(args$data.path, "data/COVIMOD/COVIMOD-single-new-hh.rds")

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

##### ---------Error table--------- #####

# extract rate and mean predictions
group_var_0 <- c("part_idx", "age", "gender", "alter_gender", "alter_age_strata")
rate_matrix <- readRDS("~/Documents/M4R/hpc_results/hsgp-eq-rd-new-hh-dropping-all-zeros-symmetric-poisson-1/rate_matrix.rds")
# setnames(rate_matrix, c("M"), c("cntct_rate_predict"))
alpha_matrix <- readRDS("~/Documents/M4R/hpc_results/hsgp-eq-rd-new-hh-dropping-all-zeros-symmetric-poisson-1/alpha_matrix.rds")
rate_matrix[, alter_age_strata := fcase(
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
alpha_matrix[, alter_age_strata := fcase(
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

individual_level=TRUE
if (individual_level){
  alpha_matrix[, cntct_intensity_predict := sum(M), by=c("part_idx", "age", "gender", "alter_gender", "alter_age_strata")]
  
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
  
  # check max(alpha_matrix[gender=="Male" & alter_gender=="Male"]$part_idx) == max(dt_intensity_MM$part_idx)
  
  # bind everything
  dt_intensity <- rbind(dt_intensity_MM, dt_intensity_FF, dt_intensity_MF, dt_intensity_FM)
  setnames(dt_intensity, c("y"), c("cntct_intensity"))
  dt_intensity <- merge(dt_intensity, alpha_matrix, by=c("part_idx", "age", "gender", "alter_age_strata", "alter_gender"), all.x = TRUE)
  dt_intensity <- unique(dt_intensity, by = group_var_0)
  
}else{
  # summing to compare with baseline model
  alpha_matrix[, cntct_intensity_predict_indiv := sum(M), by=c("part_idx", "age", "gender", "alter_gender", "alter_age_strata")]
  alpha_matrix[, cntct_intensity_predict := mean(cntct_intensity_predict_indiv), by=c("age", "gender", "alter_gender", "alter_age_strata")]
  # alpha_matrix <- unique(alpha_matrix, by=c("age", "gender", "alter_gender", "alter_age_strata"))
  # ordering dt.cnt like in model
  covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
  # making sure order of factors in alter_age_strata is ascending instead of decreasing
  # note alter_age_strata_idx and age_strata_idx are the same, but they serve different purposes
  dt.cnt[, alter_age_strata_idx:=as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
  dt.cnt<- dt.cnt[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
  dt_intensity_MM <- dt.cnt[gender=="Male"& alter_gender=="Male"]
  dt_intensity_FF <- dt.cnt[gender=="Female"& alter_gender=="Female"]
  dt_intensity_MF <- dt.cnt[gender=="Male"& alter_gender=="Female"]
  dt_intensity_FM <- dt.cnt[gender=="Female"& alter_gender=="Male"]
  
  # adding a part_idx
  dt_intensity_MM[ , part_idx := .GRP, by = new_id]      
  dt_intensity_FF[ , part_idx := .GRP, by = new_id]      
  dt_intensity_MF[ , part_idx := .GRP, by = new_id]  
  dt_intensity_FM[ , part_idx := .GRP, by = new_id]      
  
  # check max(alpha_matrix[gender=="Male" & alter_gender=="Male"]$part_idx) == max(dt_intensity_MM$part_idx)
  
  # bind everything
  dt_intensity <- rbind(dt_intensity_MM, dt_intensity_FF, dt_intensity_MF, dt_intensity_FM)
  setnames(dt_intensity, c("y"), c("cntct_intensity_indiv"))
  dt_intensity <- merge(dt_intensity, alpha_matrix, by=c("part_idx", "age", "gender", "alter_age_strata", "alter_gender"), all.x = TRUE)
  dt_intensity[, cntct_intensity:=mean(cntct_intensity_indiv), by=c("age", "gender", "alter_age_strata", "alter_gender")]
  dt_intensity <- unique(dt_intensity, by = c("age", "gender", "alter_age_strata", "alter_gender"))
  
}


error_table <- make_error_table(dt_intensity)
error_table

p3 <- ggplot(dt_intensity) +
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

p4 <- ggplot(dt_intensity) +
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


# extract empirical rate and mean from original dataset


group_var_1 <- c("age", "alter_age_strata", "gender", "alter_gender")
rate_matrix[, cntct_rate_predict := sum(M), by=group_var_1]
rate_matrix <- unique(rate_matrix, by = group_var_1)
dt_rate <- subset(dt.cnt, select=c("new_id", "alter_gender", "alter_age_strata", "age", "gender", "y", "n", "household_size"))

# divide Hic_h by strata size
dt_rate[, divide_n := fcase(
  alter_age_strata == "0-4", as.numeric(n/5),
  alter_age_strata == "5-9", as.numeric(n/5),
  alter_age_strata == "10-14", as.numeric(n/5),
  alter_age_strata == "15-19", as.numeric(n/5),
  alter_age_strata == "20-24", as.numeric(n/5),
  alter_age_strata == "25-34", as.numeric(n/10),
  alter_age_strata == "35-44", as.numeric(n/10),
  alter_age_strata == "45-54", as.numeric(n/10),
  alter_age_strata == "55-64", as.numeric(n/10),
  alter_age_strata == "65-69", as.numeric(n/5),
  alter_age_strata == "70-74", as.numeric(n/5),
  alter_age_strata == "75-79", as.numeric(n/5),
  alter_age_strata == "80-84", as.numeric(n/5),
  alter_age_strata =="85+", as.numeric(n),
  default = NA
)]

# find empirical sum of rates (summing rates over strata)
dt_rate[, cntct_rate:=y/divide_n]
# just checking everything is unique
# unique(dt_rate, by=c("new_id", "alter_gender", "alter_age_strata", "age", "gender", "y", "n", "household_size"))
# unique(dt_rate, by=c("new_id", "alter_gender", "alter_age_strata", "age", "gender"))
# unique(dt_rate, by=c("alter_gender", "alter_age_strata", "age", "gender")) --> loss of uniqueness


# dt_rate <- dt_rate[, cntct_rate := y/n]

# remove values that are nan or infinity in case (should not remove anything as we initially omit n=0 in dt.cnt)
dt_rate <- dt_rate[!is.nan(cntct_rate) & !is.infinite(cntct_rate)]

dt_rate <- merge(dt_rate, rate_matrix, by=group_var_1, all.x = TRUE)
# dt_rate <- unique(dt_rate, by = group_var_1) --> don't do unique, keep at individual level

p1 <- ggplot(dt_rate) +
  geom_tile(aes(x=age, y=alter_age_strata, fill = cntct_rate)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  coord_equal() +
  viridis::scale_fill_viridis(na.value = "white", option="H") +
  labs(x="Participants' age", y="Contacts' age ", fill="Empirical contact rate") +
  facet_grid( paste(alter_gender, "(Contacts)") ~ paste(gender, "(Participants)") ) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "bottom", text = element_text(size = 5),
        strip.background = element_rect(color=NA, fill = "transparent"))

p2 <- ggplot(dt_rate) +
  geom_tile(aes(x=age, y=alter_age_strata, fill = cntct_rate_predict)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  coord_equal() +
  viridis::scale_fill_viridis(na.value = "white", option="H") +
  labs(x="Participants' age", y="Contacts' age ", fill="Estimated contact rate") +
  facet_grid( paste(alter_gender, "(Contacts)") ~ paste(gender, "(Participants)") ) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "bottom", text = element_text(size = 5),
              strip.background = element_rect(color=NA, fill = "transparent"))


make_error_table(dt_rate, rate=TRUE)
