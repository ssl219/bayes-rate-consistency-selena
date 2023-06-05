# Preamble ----
# the simulation script aims to simulate the data set containing the total number of reports, nb of reports with detailed information by age-gender clusters.

cat("\n---------- Generating contact intensities ----------\n")

# Load the required packages
library(optparse)
library(data.table)
library(ggplot2)

# Define input arguments that can be changed by users
option_list <- list(
  optparse::make_option("--seed", type = "integer", default = 721L,
                        help = "Random number seed [default %default]",
                        dest = "seed"),
  optparse::make_option("--repo_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate/bayes-rate-consistency-selena",
                        help = "Absolute file path to repository directory",
                        dest = 'repo.path'),
  optparse::make_option("--data_path", type = "character", default = "/Users/mac/Documents/M4R/code/bayes_consistency_rate",
                        help = "Absolute file path to data directory, used as long we don t build an R package [default]",
                        dest = 'data.path')
)

# hpc repo path: "/rds/general/user/ssl219/home/bayes-rate-consistency-selena"
args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

##### ---------- Error handling ---------- #####

if (is.na(args$repo.path)) {
  stop("Please specify --repo_path !")
}

##### ---------- Setup ---------- #####
# Setting up the export directory
export.path.part1 <- "data/simulations/intensity"
export.path.part2 <- "new-hh"
args$export.path <- file.path(args$data.path, export.path.part1, export.path.part2)

if(!file.exists(args$export.path)){
  cat(paste("\n Making export directory:", args$export.path))
  dir.create(args$export.path)
}

# Directory for saving figures
fig.path <- file.path(args$export.path, "figures")
if(!dir.exists(fig.path)){
  cat(paste("\n Making figures directory:", fig.path))
  dir.create(fig.path)
}

# Source functions
source( file.path(args$repo.path, "R", "sim-intensity-utility.R") )

##### ---------- Helper functions ---------- #####
#' Simulates contact intensities/rates by age
#'
#' @param dsim
#' @param dpop Population data
#'
#' @return Simulated contact intensity/rate data.table
cntct_sim_rates_by_age <- function(dsim)
{
  dsim[, DUMMY := 1L]
  tmp <- data.table( DUMMY = 1L, gender = c('Male','Female'))
  tmpp <- data.table(DUMMY = 1L, alter_gender = c('Male','Female'))
  dsim <- merge(dsim, tmp, by = 'DUMMY', allow.cartesian = TRUE)
  dsim <- merge(dsim, tmpp, by = 'DUMMY', allow.cartesian = TRUE)
  set(dsim, NULL, c('DUMMY'), NULL) # removes those columns

  # Generate gender-gender contact rate patterns
  #
  # The contact rate cntct_rate has the symmetry property
  # Male-Female contact rate pattern and Female-Male pattern are symmetric
  # Male-Male and Female-Female patterns are self-symmetric
  # first select the whole pattern of Male-Female and convert to Female-Male pattern
  # dMF <- dsim[gender == "Male" & alter_gender == "Female"]
  # dFM <- copy(dMF)
  # setnames(dFM, c("gender","age","alter_gender","alter_age"), c("alter_gender","alter_age","gender","age"))
  
  dMF <- dsim[gender == "Male" & alter_gender == "Female" & age >= alter_age]
  tmp <- copy(dMF)
  setnames(tmp, c("age","alter_age"), c("alter_age","age"))
  dMF <- rbind(tmp[age != alter_age], dMF)
  
  dFM <- dsim[gender == "Female" & alter_gender == "Male" & age >= alter_age]
  tmp <- copy(dFM)
  setnames(tmp, c("age","alter_age"), c("alter_age","age"))
  dFM <- rbind(tmp[age != alter_age], dFM)
  
  # next select the lower triangle pattern of Male-Male and Female-Female patterns, and fill out the other part by symmetry
  dMM <- dsim[gender == "Male" & alter_gender == "Male" & age >= alter_age]
  tmp <- copy(dMM)
  setnames(tmp, c("age","alter_age"), c("alter_age","age"))
  dMM <- rbind(tmp[age != alter_age], dMM)

  dFF <- dsim[gender == "Female" & alter_gender == "Female" & age >= alter_age]
  tmp <- copy(dFF)
  setnames(tmp, c("age","alter_age"), c("alter_age","age"))
  dFF <- rbind(tmp[age != alter_age], dFF)

  dsim <- rbind(dFM, dMF, dMM, dFF)

  return(dsim)
}


# ##### ---------- Define population sizes to determine contacts intensities and rates ---------- #####
# # Set seed----
# set.seed(args$seed)
# 
# cat("\n Define population sizes to determine contacts intensities and rates ... ")
# dpop <- as.data.table(expand.grid(gender = c("Male", "Female"), age = 6:54))
# tmp <- as.data.table(read.csv(file.path(args$repo.path, "data/germany-population-2011.csv")))
# dpop <- merge(dpop, tmp, by = c("gender", "age"))

##### ---------- Generate contact intensities by age and gender ---------- #####
set.seed(args$seed)

di_flat <- cntct_sim_rates_flat()
di_boarding <- cntct_sim_rates_boarding_school()

##### ---------- Generate contact rates ---------- #####
# Generate contact rates
cat("\n Generate gender- and age-specific contact rates ...")
di_flat <- cntct_sim_rates_by_age(di_flat)
di_boarding <- cntct_sim_rates_by_age(di_boarding)

# Save the dataset
saveRDS(di_flat, file = file.path(args$export.path, "flat-data.rds"))
saveRDS(di_boarding, file = file.path(args$export.path, "boarding-data.rds"))

##### ---------- Visualization ---------- #####
setkey(di_flat, alter_age, age, alter_gender, gender)
setkey(di_boarding, alter_age, age, alter_gender, gender)

# Contact rates
p_flat <- ggplot(di_flat) +
  geom_tile(aes(x = age, y = alter_age, fill = cntct_rate)) +
  labs(x = "Participants' age", y = "Contacts' age", fill = "Contact rate, flat" ) +
  coord_equal() +
  facet_grid(alter_gender ~ gender) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  viridis::scale_fill_viridis(na.value = "white", option="H", limits=c(NA,1.3)) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 10)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = NA, fill = "transparent")
  )

# Contact rates
p_flat <- ggplot(di_flat) +
  geom_tile(aes(x = age, y = alter_age, fill = cntct_rate)) +
  labs(x = "Participants' age", y = "Contacts' age", fill = "Contact rate, flat" ) +
  coord_equal() +
  facet_grid(alter_gender ~ gender) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  viridis::scale_fill_viridis(na.value = "white", option="H", limits=c(NA,1.3)) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 10)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = NA, fill = "transparent")
  )


p_flat_2 <- ggplot(di_flat) +
  geom_tile(aes(x = age, y = alter_age, fill = cntct_rate)) +
  labs(x = "Participants' age", y = "Contacts' age", fill = "Contact rate, flat" ) +
  coord_equal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  viridis::scale_fill_viridis(na.value = "white", option="H", limits=c(NA,1.3)) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 10)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = NA, fill = "transparent")
  )

ggsave(file.path(fig.path, "flat-rates.pdf"), plot = p_flat, width = 10, height = 6)
ggsave(file.path(fig.path, "flat-rates-2.pdf"), plot = p_flat_2, width = 10, height = 6)

# Contact rates
p_boarding <- ggplot(di_boarding) +
  geom_tile(aes(x = age, y = alter_age, fill = cntct_rate)) +
  labs(x = "Participants' age", y = "Contacts' age", fill = "Contact rate, boarding" ) +
  coord_equal() +
  facet_grid(alter_gender ~ gender) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  viridis::scale_fill_viridis(na.value = "white", option="H", limits=c(NA,1.3)) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 10)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = NA, fill = "transparent")
  )

ggsave(file.path(fig.path, "boarding-rates.pdf"), plot = p_boarding, width = 10, height = 6)

# Contact rates
p_boarding_2 <- ggplot(di_boarding) +
  geom_tile(aes(x = age, y = alter_age, fill = cntct_rate)) +
  labs(x = "Participants' age", y = "Contacts' age", fill = "Contact rate, boarding" ) +
  coord_equal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  viridis::scale_fill_viridis(na.value = "white", option="H", limits=c(NA,1.3)) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 10)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = NA, fill = "transparent")
  )

ggsave(file.path(fig.path, "boarding-rates-2.pdf"), plot = p_boarding_2, width = 10, height = 6)

# # Contact intensity aggregated by age
# tmp <- di[, list(cntct_intensity = sum(cntct_intensity)), by = c('gender','age')]
# p <- ggplot(tmp, aes(x = age, y = cntct_intensity) ) +
#   geom_step(aes(colour = gender)) +
#   scale_x_continuous(expand = c(0,0)) +3
#   scale_color_brewer(palette = "Set1") +
#   labs(x = "Participants' age", y = 'Contact rate', color = 'Gender') +
#   theme_bw() +
#   theme(
#     legend.position = "bottom",
#     strip.background = element_rect(color = NA, fill = "transparent")
#   )
# 
# ggsave(file.path(fig.path, "marginal-rates.pdf"), plot = p, width = 6, height = 5)

cat("\n DONE. \n")
