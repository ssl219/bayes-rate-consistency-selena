# Creates simulated datasets from simulated intensity data
# For detailed documentation on how this works refer to the Simulating contact dataset notebook

cat("\n---------- Simulating contact dataset ----------\n")

# Load libraries
library(optparse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(viridis)
library(pammtools)
library(dplyr)
library(stringr)
library(patchwork)

option_list <- list(
  optparse::make_option("--seed", type = "integer", default = 0721, dest = "seed"),
  optparse::make_option("--size", type = "integer", default = 85,
                        help = "Number of participants with random age in the survey [default \"%default\"]",
                        dest = 'size'),
  optparse::make_option("--nsim", type = "integer", default = 10,
                        help = "Number of simulated datasets [default \"%default\"]",
                        dest = "nsim"),
  optparse::make_option("--strata", type = "character", default = "COVIMOD-new-hh",
                        help = "Age stratification scheme [default %default]",
                        dest = "strata"),
  optparse::make_option("--hhsize", type = "integer", default = 4,
                        help = "Household size [default %default]",
                        dest = "hhsize"),
  optparse::make_option("--divide.Hicb", type = "logical", default = TRUE,
                        help = "Divide Hic_b [default %default]",
                        dest = "divide.Hicb"),
  optparse::make_option("--baseline", type = "logical", default = FALSE,
                        help = "Applied on baseline model [default %default]",
                        dest = "baseline"),
  optparse::make_option("--scenario", type = "character", default = "flat",
                        help = "Scenario [default %default]",
                        dest = "scenario"),
  optparse::make_option("--random", type = "logical", default = TRUE,
                        help = "Random Hic_b [default %default]",
                        dest = "random"),
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
args$scenario = "flat"
args$hhsize = 2
args$divide.Hicb = FALSE
args$baseline = FALSE
epsilon = 1e-13
args$size = 55
args$strata = "COVIMOD"
N_random <- args$size - 55
args$random = TRUE
args$drop_zero_Hicb = TRUE

##### ---------- Error handling ---------- #####
if(is.na(args$repo.path)){
  stop("Please specify --repo_path")
}

if(is.na(args$scenario)){
  stop("Please specify --scenario")
}

###### ---------- Load data ---------- #####
source(file.path(args$repo.path, "R", "sim-dataset-utility.R"))
source(file.path(args$repo.path, "R/covimod-utility.R"))
source(file.path(args$repo.path, "R/stan-utility.R"))


if (args$scenario == "flat"){
  dt <- as.data.table(readRDS( file.path(args$data.path, "data/simulations/intensity/new-hh/flat-data.rds") ))
}

if (args$scenario == "boarding_school"){
  dt <- as.data.table(readRDS( file.path(args$data.path, "data/simulations/intensity/new-hh/boarding-data.rds") ))
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


# order alter_age strata for plots
covimod_like_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54")

##### ---------- Simulating contact survey data ---------- ##########
cat("\n Generating contact dataset ...")

# Start with male participants
# upper bound age 54 for simplicity first and computational efficiency
N = (55 + N_random)
# Generate new_id for all participants (all distinct)
d.everything.M <- as.data.table(expand.grid(new_id = 1:N, alter_age = 0:54))
d.everything.M[, DUMMY:=1L]
tmp <- data.table( DUMMY = 1L, alter_gender = c('Male','Female'))
d.everything.M <- merge(d.everything.M, tmp, by = 'DUMMY', allow.cartesian = TRUE)
d.everything.M <- stratify_alter_age(d.everything.M, args$strata)
d.everything.M[, wave:=1L]
# hhsize = args$hhsize
# d.everything.M[, household_size:=4L]
# only male participants for now
d.everything.M[, gender:="Male"]

# Generate age of participants
# set.seed(12002)
random_ages <- sample(0:54, N_random, replace=TRUE)
# part_ages <- sort(c(0:54))
part_ages <- sort(c( 0:54, random_ages))
# create smaller dataset with new_id and part_ages
id_ages <- data.table(new_id=1:N, age=part_ages)
d.everything.M <- merge(d.everything.M, id_ages, by=c("new_id"), all=TRUE)
# check no NAs, should get an empty table
# d.everything.M[is.na(d.everything.M$age)==TRUE]

# Generate household structure

# subset by alter_gender type
# alter_gender = MALE
d.everything.MM <- d.everything.M[alter_gender=="Male"]
# alter_gender = FEMALE
d.everything.MF <- d.everything.M[alter_gender=="Female"]

# generate Hic_b for male contacts
d.everything.MM[, n_M:=0L]

# number of unique ids for MM contacts
N_MM <- max(d.everything.MM$new_id)
# taking alter_age in multiples of strata present in COVIMOD

# args$hhsize=0 is the case where all Hicb are 1
if (args$hhsize==0){
    d.everything.MM[, n_M:=1]
}

# have to use age categories of the strata given in COVIMOD
if (args$hhsize==4){
  for (i in 1:N_MM){
    d.everything.MM[, n_M:=fcase(
      new_id==i & age %in% 0:9 & alter_age %in% 0:9, as.numeric(sample(c(0, 1), size=1, prob=c(0.5, 0.5))),
      new_id==i & age %in% 0:9 & alter_age %in% 30:39, as.numeric(sample(c(0, 1, 2), size=1, prob=c(0.05, 0.9, 0.05))),
      
      new_id==i & age %in% 10:19 & alter_age %in% 10:19, as.numeric(sample(c(0, 1), size=1, prob=c(0.5, 0.5))),
      new_id==i & age %in% 10:19 & alter_age %in% 40:54, as.numeric(sample(c(0, 1, 2), size=1, prob=c(0.05, 0.9, 0.05))),
      
      new_id==i & age %in% 20:29 & alter_age %in% 20:29, as.numeric(sample(c(0, 1, 2, 3), size=1, prob=c(0.25, 0.25, 0.25, 0.25))),
      
      new_id==i & age %in% 30:39 & alter_age %in% 0:9, as.numeric(sample(c(0, 1, 2), size=1, prob=c(1/3, 1/3, 1/3))),
      new_id==i & age %in% 30:39 & alter_age %in% 30:39, as.numeric(sample(c(0, 1), size=1, prob=c(0.9, 0.1))),
      
      new_id==i & age %in% 40:54 & alter_age %in% 10:19, as.numeric(sample(c(0, 1, 2), size=1, prob=c(1/3, 1/3, 1/3))),
      new_id==i & age %in% 40:54 & alter_age %in% 40:54, as.numeric(sample(c(0, 1), size=1, prob=c(0.9, 0.1))),
      
      new_id < i, as.numeric(n_M),
      default = 0
    )]
  }  
}

if (args$hhsize==2){
  for (i in 1:N_MM){
    d.everything.MM[, n_M:=fcase(
      new_id==i & age %in% 20:29 & alter_age %in% 20:29, as.numeric(sample(c(0, 1), size=1, prob=c(0.5, 0.5))),
      new_id==i & age %in% 30:39 & alter_age %in% 30:39, as.numeric(sample(c(0, 1), size=1, prob=c(0.5, 0.5))),
      new_id==i & age %in% 40:54 & alter_age %in% 40:54, as.numeric(sample(c(0, 1), size=1, prob=c(0.5, 0.5))),
      new_id < i, as.numeric(n_M),
      default = 0
    )]
  }  
}

if (args$hhsize==8){
  for (i in 1:N_MM){
    d.everything.MM[, n_M:=fcase(
      new_id==i & age %in% 0:9 & alter_age %in% 0:9, as.numeric(sample(c(0, 1, 2, 3, 4, 5), size=1, prob=c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6))),
      new_id==i & age %in% 0:9 & alter_age %in% 30:39, as.numeric(sample(c(0, 1, 2), size=1, prob=c(0.05, 0.9, 0.05))),
      
      new_id==i & age %in% 10:19 & alter_age %in% 10:19, as.numeric(sample(c(0, 1, 2, 3, 4, 5), size=1, prob=c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6))),
      new_id==i & age %in% 10:19 & alter_age %in% 40:54, as.numeric(sample(c(0, 1, 2), size=1, prob=c(0.05, 0.9, 0.05))),
      
      new_id==i & age %in% 20:29 & alter_age %in% 20:29, as.numeric(sample(c(0, 1, 2, 3, 4, 5, 6, 7), size=1, prob=c(1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8))),
      
      new_id==i & age %in% 30:39 & alter_age %in% 0:9, as.numeric(sample(c(0, 1, 2, 3, 4, 5), size=1, prob=c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6))),
      new_id==i & age %in% 30:39 & alter_age %in% 30:39, as.numeric(sample(c(0, 1), size=1, prob=c(0.9, 0.1))),
      
      new_id==i & age %in% 40:54 & alter_age %in% 10:19, as.numeric(sample(c(0, 1, 2, 3, 4, 5), size=1, prob=c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6))),
      new_id==i & age %in% 40:54 & alter_age %in% 40:54, as.numeric(sample(c(0, 1), size=1, prob=c(0.9, 0.1))),
      
      new_id < i, as.numeric(n_M),
      default = 0
      )]
  }  
}

# generate Hic_b for female contacts
d.everything.MF[, n_F:=0L]

# get n_M from male contacts
d.everything.MF <- merge(d.everything.MF, d.everything.MM[, c("new_id", "n_M", "age", "alter_age")], by=c("new_id", "age", "alter_age"), all=TRUE)
# number of unique ids for MM contacts
N_MF <- max(d.everything.MF$new_id)

if (args$hhsize==0){
d.everything.MF[, n_F:=1]
    
  if (args$divide.Hicb){
  # Set Hic_b to 1
  d.everything.MM[, Hic_b := fcase(
    alter_age_strata == "0-4", as.numeric(n_M/5),
    alter_age_strata == "5-9", as.numeric(n_M/5),
    alter_age_strata == "10-14", as.numeric(n_M/5),
    alter_age_strata == "15-19", as.numeric(n_M/5),
    alter_age_strata == "20-24", as.numeric(n_M/5),
    alter_age_strata == "25-34", as.numeric(n_M/10),
    alter_age_strata == "35-44", as.numeric(n_M/10),
    alter_age_strata == "45-54", as.numeric(n_M/5),
    # alter_age_strata == "50-64", as.numeric(n_M/15),
    # alter_age_strata == "65-69", as.numeric(n_M/5),
    # alter_age_strata == "70-74", as.numeric(n_M/5),
    # alter_age_strata == "75-79", as.numeric(n_M/5),
    # alter_age_strata == "80-84", as.numeric(n_M/5),
    # alter_age_strata =="85+", as.numeric(n_M),
    default = NA
  )]
  
  d.everything.MF[, Hic_b := fcase(
    alter_age_strata == "0-4", as.numeric(n_M/5),
    alter_age_strata == "5-9", as.numeric(n_M/5),
    alter_age_strata == "10-14", as.numeric(n_M/5),
    alter_age_strata == "15-19", as.numeric(n_M/5),
    alter_age_strata == "20-24", as.numeric(n_M/5),
    alter_age_strata == "25-34", as.numeric(n_M/10),
    alter_age_strata == "35-44", as.numeric(n_M/10),
    alter_age_strata == "45-54", as.numeric(n_M/5),
    # alter_age_strata == "50-64", as.numeric(n_M/15),
    # alter_age_strata == "65-69", as.numeric(n_M/5),
    # alter_age_strata == "70-74", as.numeric(n_M/5),
    # alter_age_strata == "75-79", as.numeric(n_M/5),
    # alter_age_strata == "80-84", as.numeric(n_M/5),
    # alter_age_strata =="85+", as.numeric(n_M),
    default = NA
  )]
  }
else{
  d.everything.MM[, Hic_b := 1]
  d.everything.MF[, Hic_b := 1] 
}
}


if (args$hhsize==4){
for (i in 1:N_MF){
  d.everything.MF[, n_F:=fcase(
    new_id==i & age %in% 0:9 & alter_age %in% 0:9, as.numeric(1-n_M),
    new_id==i & age %in% 10:19 & alter_age %in% 10:19, as.numeric(1-n_M),
    
    new_id==i & age %in% 0:9 & alter_age %in% 30:39, as.numeric(2-n_M),
    new_id==i & age %in% 10:19 & alter_age %in% 40:54, as.numeric(2-n_M),
    
    new_id==i & age %in% 20:29 & alter_age %in% 20:29, as.numeric(3-n_M),
    
    new_id==i & age %in% 30:39 & alter_age %in% 0:9, as.numeric(2-n_M),
    new_id==i & age %in% 40:54 & alter_age %in% 10:19, as.numeric(2-n_M),
    
    new_id==i & age %in% 30:39 & alter_age %in% 30:39, as.numeric(1-n_M),
    new_id==i & age %in% 40:54 & alter_age %in% 40:54, as.numeric(1-n_M),
    new_id < i, as.numeric(n_F),
    default = 0
  )]
}  
}

if (args$hhsize==2){
  for (i in 1:N_MF){
    d.everything.MF[, n_F:=fcase(
      new_id==i & age %in% 20:29 & alter_age %in% 20:29, as.numeric(1-n_M),
      new_id==i & age %in% 30:39 & alter_age %in% 30:39, as.numeric(1-n_M),
      new_id==i & age %in% 40:54 & alter_age %in% 40:54, as.numeric(1-n_M),
      new_id < i, as.numeric(n_F),
      default = 0
    )]
  }  
}

if (args$hhsize==8){
  for (i in 1:N_MF){

    d.everything.MF[, n_F:=fcase(
      new_id==i & age %in% 0:9 & alter_age %in% 0:9, as.numeric(5-n_M),
      new_id==i & age %in% 10:19 & alter_age %in% 10:19, as.numeric(5-n_M),
      
      new_id==i & age %in% 0:9 & alter_age %in% 30:39, as.numeric(2-n_M),
      new_id==i & age %in% 10:19 & alter_age %in% 40:54, as.numeric(2-n_M),
      
      new_id==i & age %in% 20:29 & alter_age %in% 20:29, as.numeric(7-n_M),
      
      new_id==i & age %in% 30:39 & alter_age %in% 0:9, as.numeric(5-n_M),
      new_id==i & age %in% 40:54 & alter_age %in% 10:19, as.numeric(5-n_M),
      
      new_id==i & age %in% 30:39 & alter_age %in% 30:39, as.numeric(1-n_M),
      new_id==i & age %in% 40:54 & alter_age %in% 40:54, as.numeric(1-n_M),
      new_id < i, as.numeric(n_F),
      default = 0
    )]
  }  
}



  if (args$hhsize!=0){
    # divide by size of interval to get Hic_b
    d.everything.MM[, Hic_b:=fcase(
      alter_age < 40, as.numeric(n_M/10),
      alter_age >= 40, as.numeric(n_M/15),
      default = 0
    )]
    
    d.everything.MF[, Hic_b:=fcase(
      alter_age < 40, as.numeric(n_M/10),
      alter_age >= 40, as.numeric(n_M/15),
      default = 0
    )]
    
  }
  

# remove unwanted columns
set(d.everything.MM, NULL, c('DUMMY'), NULL)
set(d.everything.MF, NULL, c('n_M', 'DUMMY'), NULL)
# change names
setnames(d.everything.MM, c('n_M'), c('n') )
setnames(d.everything.MF, c('n_F'), c('n') )

d.everything.M.final <- rbind(d.everything.MM, d.everything.MF)
# DO NOT TAKE OUT 0 HIC_B
# d.everything.M.final <- d.everything.M.final[Hic_b!=0]
# if N_random = 0, should get a dataset of dimension 85*85*2 = 14450
#################################################################################################################################


# do the same for female participants
N = (55 + N_random)
# Generate new_id for all participants (all distinct)
d.everything.F <- as.data.table(expand.grid(new_id = (N+1):(2*N), alter_age = 0:54))
d.everything.F[, DUMMY:=1L]
tmp <- data.table( DUMMY = 1L, alter_gender = c('Male','Female'))
d.everything.F <- merge(d.everything.F, tmp, by = 'DUMMY', allow.cartesian = TRUE)
d.everything.F <- stratify_alter_age(d.everything.F, args$strata)
d.everything.F[, wave:=1L]
# hhsize = args$hhsize
# d.everything.F[, household_size:=4L]
# female participants now
d.everything.F[, gender:="Female"]

# Generate age of participants
set.seed(2002)
random_ages <- sample(0:54, N_random, replace=TRUE)
part_ages <- sort(c( 0:54, random_ages))
# part_ages <- sort(c(0:54))
# create smaller dataset with new_id and part_ages
id_ages <- data.table(new_id=(N+1):(2*N), age=part_ages)
d.everything.F <- merge(d.everything.F, id_ages, by=c("new_id"), all=TRUE)
# check no NAs, should get an empty table
# d.everything.F[is.na(d.everything.F$age)==TRUE]

# Generate household structure

# subset by alter_gender type
# alter_gender = MALE
d.everything.FM <- d.everything.F[alter_gender=="Male"]
# alter_gender = FEMALE
d.everything.FF <- d.everything.F[alter_gender=="Female"]

# generate Hic_b for male contacts
d.everything.FM[, n_M:=0L]

# number of unique ids for MM contacts
N_FM_U <- max(d.everything.FM$new_id)
N_FM_L <- min(d.everything.FM$new_id)

if (args$hhsize==4){
  # taking alter_age in multiples of strata present in COVIMOD
  for (i in N_FM_L:N_FM_U){
    d.everything.FM[, n_M:=fcase(
      new_id==i & age %in% 0:9 & alter_age %in% 0:9, as.numeric(sample(c(0, 1), size=1, prob=c(0.5, 0.5))),
      new_id==i & age %in% 0:9 & alter_age %in% 30:39, as.numeric(sample(c(0, 1, 2), size=1, prob=c(0.05, 0.9, 0.05))),
      
      new_id==i & age %in% 10:19 & alter_age %in% 10:19, as.numeric(sample(c(0, 1), size=1, prob=c(0.5, 0.5))),
      new_id==i & age %in% 10:19 & alter_age %in% 40:54, as.numeric(sample(c(0, 1, 2), size=1, prob=c(0.05, 0.9, 0.05))),
      
      new_id==i & age %in% 20:29 & alter_age %in% 20:29, as.numeric(sample(c(0, 1, 2, 3), size=1, prob=c(0.25, 0.25, 0.25, 0.25))),
      
      new_id==i & age %in% 30:39 & alter_age %in% 0:9, as.numeric(sample(c(0, 1, 2), size=1, prob=c(1/3, 1/3, 1/3))),
      new_id==i & age %in% 30:39 & alter_age %in% 30:39, as.numeric(sample(c(0, 1), size=1, prob=c(0.1, 0.9))),
      
      new_id==i & age %in% 40:54 & alter_age %in% 10:19, as.numeric(sample(c(0, 1, 2), size=1, prob=c(1/3, 1/3, 1/3))),
      new_id==i & age %in% 40:54 & alter_age %in% 40:54, as.numeric(sample(c(0, 1), size=1, prob=c(0.1, 0.9))),
      
      new_id < i, as.numeric(n_M),
      default = 0
    )]
  }  
}

if (args$hhsize==2){
  for (i in N_FM_L:N_FM_U){
    d.everything.FM[, n_M:=fcase(
      new_id==i & age %in% 20:29 & alter_age %in% 20:29, as.numeric(sample(c(0, 1), size=1, prob=c(0.5, 0.5))),
      new_id==i & age %in% 30:39 & alter_age %in% 30:39, as.numeric(sample(c(0, 1), size=1, prob=c(0.5, 0.5))),
      new_id==i & age %in% 40:54 & alter_age %in% 40:54, as.numeric(sample(c(0, 1), size=1, prob=c(0.5, 0.5))),
      new_id < i, as.numeric(n_M),
      default = 0
    )]
  }  
}

if (args$hhsize==8){
  for (i in N_FM_L:N_FM_U){
    d.everything.FM[, n_M:=fcase(
      new_id==i & age %in% 0:9 & alter_age %in% 0:9, as.numeric(sample(c(0, 1, 2, 3, 4, 5), size=1, prob=c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6))),
      new_id==i & age %in% 0:9 & alter_age %in% 30:39, as.numeric(sample(c(0, 1, 2), size=1, prob=c(0.05, 0.9, 0.05))),
      
      new_id==i & age %in% 10:19 & alter_age %in% 10:19, as.numeric(sample(c(0, 1, 2, 3, 4, 5), size=1, prob=c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6))),
      new_id==i & age %in% 10:19 & alter_age %in% 40:54, as.numeric(sample(c(0, 1, 2), size=1, prob=c(0.05, 0.9, 0.05))),
      
      new_id==i & age %in% 20:29 & alter_age %in% 20:29, as.numeric(sample(c(0, 1, 2, 3, 4, 5, 6, 7), size=1, prob=c(1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8))),
      
      new_id==i & age %in% 30:39 & alter_age %in% 0:9, as.numeric(sample(c(0, 1, 2, 3, 4, 5), size=1, prob=c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6))),
      new_id==i & age %in% 30:39 & alter_age %in% 30:39, as.numeric(sample(c(0, 1), size=1, prob=c(0.9, 0.1))),
      
      new_id==i & age %in% 40:54 & alter_age %in% 10:19, as.numeric(sample(c(0, 1, 2, 3, 4, 5), size=1, prob=c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6))),
      new_id==i & age %in% 40:54 & alter_age %in% 40:54, as.numeric(sample(c(0, 1), size=1, prob=c(0.9, 0.1))),
      
      new_id < i, as.numeric(n_M),
      default = 0
    )]
  }  
}

# generate Hic_b for female contacts
d.everything.FF[, n_F:=0L]

# get n_M from male contacts
d.everything.FF <- merge(d.everything.FF, d.everything.FM[, c("new_id", "n_M", "age", "alter_age")], by=c("new_id", "age", "alter_age"), all=TRUE)
# number of unique ids for MM contacts
N_FF_U <- max(d.everything.FF$new_id)
N_FF_L <- min(d.everything.FF$new_id)


if (args$hhsize==0){
  for (i in N_FM_L:N_FM_U){
    d.everything.FM[, n_M:=1]
  }  
  
  for (i in N_FF_L:N_FF_U){
    d.everything.FF[, n_F:=1]
  }  
  
  if (args$divide.Hicb){
    # SET EVERYTHING TO 1
    d.everything.FM[, Hic_b:=fcase(
      alter_age_strata == "0-4", as.numeric(n_M/5),
      alter_age_strata == "5-9", as.numeric(n_M/5),
      alter_age_strata == "10-14", as.numeric(n_M/5),
      alter_age_strata == "15-19", as.numeric(n_M/5),
      alter_age_strata == "20-24", as.numeric(n_M/5),
      alter_age_strata == "25-34", as.numeric(n_M/10),
      alter_age_strata == "35-44", as.numeric(n_M/10),
      alter_age_strata == "45-54", as.numeric(n_M/5),
      # alter_age_strata == "50-64", as.numeric(n_M/15),
      # alter_age_strata == "65-69", as.numeric(n_M/5),
      # alter_age_strata == "70-74", as.numeric(n_M/5),
      # alter_age_strata == "75-79", as.numeric(n_M/5),
      # alter_age_strata == "80-84", as.numeric(n_M/5),
      # alter_age_strata =="85+", as.numeric(n_M),
      default = NA
    )]
    
    d.everything.FF[, Hic_b:=fcase(
      alter_age_strata == "0-4", as.numeric(n_M/5),
      alter_age_strata == "5-9", as.numeric(n_M/5),
      alter_age_strata == "10-14", as.numeric(n_M/5),
      alter_age_strata == "15-19", as.numeric(n_M/5),
      alter_age_strata == "20-24", as.numeric(n_M/5),
      alter_age_strata == "25-34", as.numeric(n_M/10),
      alter_age_strata == "35-44", as.numeric(n_M/10),
      alter_age_strata == "45-54", as.numeric(n_M/5),
      # alter_age_strata == "50-64", as.numeric(n_M/15),
      # alter_age_strata == "65-69", as.numeric(n_M/5),
      # alter_age_strata == "70-74", as.numeric(n_M/5),
      # alter_age_strata == "75-79", as.numeric(n_M/5),
      # alter_age_strata == "80-84", as.numeric(n_M/5),
      # alter_age_strata =="85+", as.numeric(n_M),
      default = NA
    )]
  }
  
  else{
    d.everything.FM[, Hic_b:=1]
    d.everything.FF[, Hic_b:=1]
  }
}

if (args$hhsize==4){
for (i in N_FF_L:N_FF_U){
  d.everything.FF[, n_F:=fcase(
    new_id==i & age %in% 0:9 & alter_age %in% 0:9, as.numeric(1-n_M),
    new_id==i & age %in% 10:19 & alter_age %in% 10:19, as.numeric(1-n_M),
    
    new_id==i & age %in% 0:9 & alter_age %in% 30:39, as.numeric(2-n_M),
    new_id==i & age %in% 10:19 & alter_age %in% 40:54, as.numeric(2-n_M),
    
    new_id==i & age %in% 20:29 & alter_age %in% 20:29, as.numeric(3-n_M),
    
    new_id==i & age %in% 30:39 & alter_age %in% 0:9, as.numeric(2-n_M),
    new_id==i & age %in% 40:54 & alter_age %in% 10:19, as.numeric(2-n_M),
    
    new_id==i & age %in% 30:39 & alter_age %in% 30:39, as.numeric(1-n_M),
    new_id==i & age %in% 40:54 & alter_age %in% 40:54, as.numeric(1-n_M),
    new_id < i, as.numeric(n_F),
    default = 0
  )]
}  
}

if (args$hhsize==2){
  for (i in N_FF_L:N_FF_U){
    d.everything.FF[, n_F:=fcase(
      new_id==i & age %in% 20:29 & alter_age %in% 20:29, as.numeric(1-n_M),
      new_id==i & age %in% 30:39 & alter_age %in% 30:39, as.numeric(1-n_M),
      new_id==i & age %in% 40:54 & alter_age %in% 40:54, as.numeric(1-n_M),
      new_id < i, as.numeric(n_F),
      default = 0
    )]
  }  
}

if (args$hhsize==8){
  for (i in  N_FF_L:N_FF_U){
    d.everything.FF[, n_F:=fcase(
      new_id==i & age %in% 0:9 & alter_age %in% 0:9, as.numeric(5-n_M),
      new_id==i & age %in% 10:19 & alter_age %in% 10:19, as.numeric(5-n_M),
      
      new_id==i & age %in% 0:9 & alter_age %in% 30:39, as.numeric(2-n_M),
      new_id==i & age %in% 10:19 & alter_age %in% 40:54, as.numeric(2-n_M),
      
      new_id==i & age %in% 20:29 & alter_age %in% 20:29, as.numeric(7-n_M),
      
      new_id==i & age %in% 30:39 & alter_age %in% 0:9, as.numeric(5-n_M),
      new_id==i & age %in% 40:54 & alter_age %in% 10:19, as.numeric(5-n_M),
      
      new_id==i & age %in% 30:39 & alter_age %in% 30:39, as.numeric(1-n_M),
      new_id==i & age %in% 40:54 & alter_age %in% 40:54, as.numeric(1-n_M),
      new_id < i, as.numeric(n_F),
      default = 0
    )]
  }  
}

if (args$hhsize!=0){
# divide by size of interval to get Hic_b
d.everything.FM[, Hic_b:=fcase(
  alter_age < 40, as.numeric(n_M/10),
  alter_age >= 40, as.numeric(n_M/15),
  default = 0
)]

d.everything.FF[, Hic_b:=fcase(
  alter_age < 40, as.numeric(n_M/10),
  alter_age >= 40, as.numeric(n_M/15),
  default = 0
)]

}

# remove unwanted columns
set(d.everything.FM, NULL, c('DUMMY'), NULL)
set(d.everything.FF, NULL, c('n_M', 'DUMMY'), NULL)
# change names
setnames(d.everything.FM, c('n_M'), c('n') )
setnames(d.everything.FF, c('n_F'), c('n') )

d.everything.F.final <- rbind(d.everything.FM, d.everything.FF)
# DO NOT TAKE OUT ZERO HIC_B
# d.everything.F.final <- d.everything.F.final[Hic_b!=0]

# bind female and male participants
d.everything.final <- rbind(d.everything.M.final, d.everything.F.final)

# take random Hic_b instead of smoothing out

# note that for the case hhsize=0, already have individual Hic_b
if (args$random == TRUE & args$hhsize!=0){
  d.everything.final[, hh_structure_strata := fcase(
    alter_age %in% 0:9, "0-9",
    alter_age %in% 10:19, "10-19",
    alter_age %in% 20:29, "20-29",
    alter_age %in% 30:39, "30-39",
    alter_age %in% 40:54, "40-54",
    default = NA
  )]
  hh_structure_strata_comb <- d.everything.final %>% distinct(new_id, hh_structure_strata, alter_gender,
                                                              wave, .keep_all=TRUE)
  hh_structure_strata_comb_list <- hh_structure_strata_comb$hh_structure_strata
  extract_hh_structure_strata_comb <- str_extract_all(as.character(hh_structure_strata_comb_list), "[0-9]+")
  set.seed(1702652)
  hh_structure_strata_random = rep(0, length(hh_structure_strata_comb_list))
  for (i in 1:length(hh_structure_strata_comb_list)){
    hh_structure_strata_random[i] = sample(extract_hh_structure_strata_comb[[i]][1]:extract_hh_structure_strata_comb[[i]][2], size=1)
    # print(alter_age_strata_random[i])
  }
  hh_structure_strata_comb[,alter_age:=hh_structure_strata_random]
  # amend alter_age_strata (in case alter_age was chosen within hh_structure_strata but outside of alter_age_strata)
  stratify_alter_age(hh_structure_strata_comb, args$strata)
  
  setnames(hh_structure_strata_comb,c("Hic_b"),c("hh_structure_strata_Hic_b"))
  setnames(hh_structure_strata_comb, c("n"), c("Hic_b"))
  setnames(d.everything.final, c("Hic_b"), c("hh_structure_strata_Hic_b"))
  d.everything.final <- merge(d.everything.final, hh_structure_strata_comb, by = c("alter_age", "new_id", "alter_gender", "gender", "wave", "age", "alter_age_strata", "hh_structure_strata", "hh_structure_strata_Hic_b"), all.x=TRUE, all.y=TRUE)
  d.everything.final[is.na(Hic_b), Hic_b:=0]
  }

# calculate means alpha
d.everything.final <- merge(d.everything.final, dt, by=c("alter_age","age", "gender", "alter_gender", "alter_age_strata"), all.x=TRUE, all.y=FALSE)
d.everything.final[, alpha:=Hic_b*cntct_rate + epsilon]

# plot HH structure
d.everything.final[, alpha_hh_structure_plot:= mean(Hic_b*1 + epsilon), by=c("alter_age","age", "gender", "alter_gender")]
d.everything.final[, gender_comb := fcase(
  gender == "Male" & alter_gender == "Male", "Male to Male",
  gender == "Male" & alter_gender == "Female", "Male to Female",
  gender == "Female" & alter_gender == "Male", "Female to Male",
  gender == "Female" & alter_gender == "Female", "Female to Female",
  default = NA
)]

# order is FF, FM, MF, MM, do not omit for patchwork later
alpha_hh_structure_plot <- function(){
  ggplot(d.everything.final, aes(age, alter_age)) + 
  geom_tile(aes(fill = alpha_hh_structure_plot)) + 
  facet_wrap(~gender_comb, ncol = 4, nrow = 1) + 
  viridis::scale_fill_viridis(option='H') +
  theme_bw() + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Age of contacts", fill = "Average count") +
  guides(fill = guide_colourbar(barwidth = 0.8, barheight=4)) +
    theme( aspect.ratio = 1,
           axis.text.x = element_text(size = 7),
           axis.text.y = element_text(size = 7),
           axis.title.x = element_blank(),
           axis.title.y = element_text(size = 8),
           strip.background = element_blank(),
           # strip.text = element_blank(),
           strip.text = element_text(size = 8),
           legend.text = element_text(size = 8),
           legend.title = element_text(size = 8),
           legend.margin = margin(l = -6, unit = "cm"),
           plot.margin = margin())
}
# simulate contact counts
d.everything.final[, y := rpois(nrow(d.everything.final), lambda=d.everything.final$alpha)]
# take out rows with no hh structure, i.e. Hic_b = 0
# # fill in Hic_b
# for (i in 1:N_FM){
#   if (d.everything.FM$age[1, d.everything.FM$new_id==i] %in% 0:20){
#     d.everything.FM$Hic_b[d.everything.FM[, d.everything.FM$new_id==i & d.everything.FM$alter_age %in% 0:20]] <-  sample(c(0, 1), size=1, replace=FALSE, prob=c(0.5, 0.5)),
#   }
#   
# }

# Stratify contact intensities and contact rates
group_var_plot <- c("age", "gender", "alter_age_strata", "alter_gender")
group_var_indiv <- c("new_id", "age", "gender", "alter_age_strata", "alter_gender")
d_comb_no_dupl_plot <- copy(d.everything.final)
d_comb_no_dupl_plot[, y_plot_strata:= sum(y), by=group_var_plot]
d_comb_no_dupl_plot[, y_plot:= mean(y), by=group_var_plot]
d_comb_no_dupl_plot[, y_plot_indiv:= mean(y), by=group_var_indiv]

# add empirical contact rates for household model
# d_comb_no_dupl_plot[, emp_cntct_rates_indiv := y_plot / alpha_hh_structure_plot]
d_comb_no_dupl_plot[, emp_cntct_rates_indiv := y_plot_indiv / Hic_b]
d_comb_no_dupl_plot[is.nan(emp_cntct_rates_indiv), emp_cntct_rates_indiv := 0]
d_comb_no_dupl_plot[is.infinite(emp_cntct_rates_indiv), emp_cntct_rates_indiv := 0]
# d_comb_no_dupl_plot[emp_cntct_rates_indiv > 50, emp_cntct_rates_indiv := 0]
d_comb_no_dupl_plot[, mean_emp_cntct_rates_indiv := mean(emp_cntct_rates_indiv), by=c("age", "gender", "alter_age", "alter_gender")]
# d_comb_no_dupl_plot[is.na(mean_emp_cntct_rates_indiv), mean_emp_cntct_rates_indiv := 0]

# smooth out Hic_b to see effect 
d_comb_strata_plot <- copy(d_comb_no_dupl_plot)
d_comb_strata_plot[, strata_Hic_b_indiv:=mean(Hic_b), by=group_var_indiv]
# d_comb_strata_plot[, strata_Hic_b:=mean(Hic_b), by=group_var_plot]
d_comb_strata_plot[, emp_cntct_rates_indiv := y_plot_indiv / strata_Hic_b_indiv]
d_comb_strata_plot[is.nan(emp_cntct_rates_indiv), emp_cntct_rates_indiv := 0]
d_comb_strata_plot[is.infinite(emp_cntct_rates_indiv), emp_cntct_rates_indiv := 0]
d_comb_strata_plot[, mean_emp_cntct_rates_indiv := mean(emp_cntct_rates_indiv), by=c("age", "gender", "alter_age", "alter_gender")]
# d_comb_strata_plot[is.na(mean_emp_cntct_rates_indiv), mean_emp_cntct_rates_indiv := 0]

# add empirical contacts for baseline model
group_var_offset <- c("wave", "age", "gender")
d_comb_no_dupl_plot[, n_baseline := sum(n), by=group_var_offset]
d_comb_no_dupl_plot[, emp_cntct_int_baseline := y_plot / n_baseline]
d_comb_no_dupl_plot[is.nan(emp_cntct_int_baseline), emp_cntct_int_baseline := 0]

# merge with d_comb_no_dupl_plot
covimod <- load_covimod_data("/Users/mac/Documents/M4R/code")
dt.pop.plot <- covimod$pop
setnames(dt.pop.plot, "age", "alter_age")
d_comb_no_dupl_plot <- merge(d_comb_no_dupl_plot, dt.pop.plot, by = c("alter_age", "gender"), all.x=TRUE, all.y=FALSE)
d_comb_no_dupl_plot[, emp_cntct_rates_baseline := emp_cntct_int_baseline /pop]

saveRDS(d_comb_strata_plot, file=file.path(export.path, paste0("d_comb_strata_plot-hh", args$hhsize, "-", args$scenario, "-", args$size, ".rds")))
saveRDS(d_comb_no_dupl_plot, file=file.path(export.path, paste0("d_comb_no_dupl_plot-hh", args$hhsize, "-", args$scenario, "-", args$size, ".rds")))



simulated_counts_per_gender <- function(){
  ggplot(d_comb_no_dupl_plot, aes(age, factor(alter_age_strata, levels=covimod_like_strata_levels ))) +
    geom_tile(aes(fill = y_plot_strata)) +
    facet_wrap(~gender_comb, ncol = 4, nrow = 1) + 
    scale_fill_viridis(begin=0.3, option = "rocket") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    labs(x = "Age of participants", y = "Age of contacts", fill = "Count") +
    theme_bw() + 
    guides(fill = guide_colourbar(barwidth = 0.8, barheight=4)) +
    theme( aspect.ratio = 1,
           axis.text.x = element_text(size = 7),
           axis.text.y = element_text(size = 7),
           axis.title.x = element_blank(),
           axis.title.y = element_text(size = 8),
           strip.background = element_blank(),
           strip.text = element_blank(),
           # strip.text = element_text(size = 8),
           legend.text = element_text(size = 8),
           legend.title = element_text(size = 8),
           legend.margin = margin(l = -6.9, unit = "cm"),
           plot.margin = margin())
}

empirical_rates_new_hh <- function(){
  ggplot(d_comb_no_dupl_plot, aes(age, alter_age)) +
    geom_tile(aes(fill = mean_emp_cntct_rates_indiv)) +
    facet_wrap(~gender_comb, ncol = 4, nrow = 1) + 
    scale_fill_viridis(option = "H") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Age of participants", y = "Age of contacts", fill = "Intensity") +
    theme_bw() + 
    guides(fill = guide_colourbar(barwidth = 0.8, barheight=4)) +
    theme( aspect.ratio = 1,
           axis.text.x = element_text(size = 7),
           axis.text.y = element_text(size = 7),
           axis.title.x = element_blank(),
           axis.title.y = element_text(size = 8),
           strip.background = element_blank(),
           strip.text = element_blank(),
           # strip.text = element_text(size = 8),
           legend.text = element_text(size = 8),
           legend.title = element_text(size = 8),
           legend.margin = margin(l = -6.7, unit = "cm"),
           plot.margin = margin())
}

empirical_rates_new_hh_smooth_offset <- function(){
  ggplot(d_comb_strata_plot, aes(age, alter_age)) +
    geom_tile(aes(fill = mean_emp_cntct_rates_indiv)) +
    facet_wrap(~gender_comb, ncol = 4, nrow = 1) + 
    scale_fill_viridis(option = "H") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Age of participants", y = "Age of contacts", fill = "Intensity") +
    theme_bw() + 
    guides(fill = guide_colourbar(barwidth = 0.8, barheight=4)) +
    theme( aspect.ratio = 1,
           axis.text.x = element_text(size = 7),
           axis.text.y = element_text(size = 7),
           axis.title.x = element_blank(),
           axis.title.y = element_text(size = 8),
           strip.background = element_blank(),
           strip.text = element_blank(),
           # strip.text = element_text(size = 8),
           legend.text = element_text(size = 8),
           legend.title = element_text(size = 8),
           legend.margin = margin(l = -6.7, unit = "cm"),
           plot.margin = margin())
}


empirical_rates_baseline <- function(){
  ggplot(d_comb_no_dupl_plot, aes(age, alter_age)) +
    geom_tile(aes(fill = emp_cntct_rates_baseline)) +
    facet_wrap(~gender_comb, ncol = 4, nrow = 1) + 
    scale_fill_viridis(option = "H") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Age of participants", y = "Age of contacts", fill = "Intensity") +
    theme_bw() + 
    guides(fill = guide_colourbar(barwidth = 0.8, barheight=4)) +
    theme( aspect.ratio = 1,
           axis.text.x = element_text(size = 7),
           axis.text.y = element_text(size = 7),
           axis.title.x = element_text(size = 8),
           axis.title.y = element_text(size = 8),
           strip.background = element_blank(),
           strip.text = element_blank(),
           # strip.text = element_text(size = 8),
           legend.text = element_text(size = 8),
           legend.title = element_text(size = 8),
           legend.margin = margin(l = -6.5, unit = "cm"),
           plot.margin = margin())
  
}

if (!args$divide.Hicb){
  
  alpha_hh_structure_plot() / simulated_counts_per_gender() / empirical_rates_new_hh() / empirical_rates_new_hh_smooth_offset() / empirical_rates_baseline() + plot_layout(nrow = 5)
  ggsave(file.path(export.path, paste0("hh", args$hhsize, "-", args$scenario, "-", args$size, "-full.pdf")), width = 20, height = 16, units = "cm")
}else{
  
  alpha_hh_structure_plot() / simulated_counts_per_gender() / empirical_rates_new_hh() / empirical_rates_new_hh_smooth_offset() / empirical_rates_baseline() + plot_layout(nrow = 5)
  ggsave(file.path(export.path, paste0("hh", args$hhsize, "-", args$scenario, "-", args$size, "-full.pdf")), width = 20, height = 16, units = "cm")
  
}

# baseline model
group_var_baseline <- c("age", "gender", "alter_age_strata", "alter_gender")
d_comb_no_dupl_baseline <- copy(d.everything.final)
d_comb_no_dupl_baseline[, y_strata := sum(y), by=group_var_baseline]
d_comb_no_dupl_baseline[, cntct_rate_strata := sum(cntct_rate), by=group_var_baseline]
# add u for row major indexing
d_comb_no_dupl_baseline[, u:=1L]

d_comb_no_dupl_baseline <- d_comb_no_dupl_baseline %>% distinct(wave, age, alter_age_strata, gender, alter_gender, .keep_all=TRUE)
# if N_random = 0, should get dataset of dimension 85*2*13*2 = 4420, 
# 85* because 85 participants/ages, *2* because of 2 alter_genders, *13* because of strata, last *2 is for male/female participants, 
setnames(d_comb_no_dupl_baseline, c("y", "y_strata"), c("y_age", "y"))

# aggregate offsets by age and gender
group_var_offset_baseline <- c("wave", "age", "gender")
d.everything.final_baseline <- copy(d.everything.final)
d.everything.final_baseline[, n := sum(n), by=group_var_offset_baseline]
d.everything.final_baseline[, zeta := 1L]
d.everything.final_baseline <- d.everything.final_baseline %>% distinct(wave, age, gender, .keep_all=TRUE)
# N and zeta to the offsets dataset
setnames(d.everything.final_baseline, c("n"), c("N"))

if (args$drop_zero_Hicb){
  # omit Hic_b = 0 to be comparable with new-hh model
  d.everything.final <- d.everything.final[Hic_b !=0, ]
}

# new-hh model
  group_var <- c("new_id", "age", "gender", "alter_age_strata", "alter_gender")
  d_comb_no_dupl <- copy(d.everything.final)
  d_comb_no_dupl[, y_strata := sum(y), by=group_var]
  d_comb_no_dupl[, cntct_rate_strata := sum(cntct_rate), by=group_var]
  d_comb_no_dupl <- d_comb_no_dupl %>% distinct(new_id, wave, alter_age_strata, alter_gender, .keep_all=TRUE)
  # omit Hic_b = 0
  # if N_random = 0, should get dataset of dimension 85*2*13*2 = 4420, 
  # 85* because 85 participants/ages, *2* because of 2 alter_genders, *13* because of strata, last *2 is for male/female participants, 
  setnames(d_comb_no_dupl, c("y", "y_strata"), c("y_age", "y"))






# d_comb_no_dupl[, alter_age_strata_idx_simplot:=as.numeric(factor(alter_age_strata, levels=covimod_like_strata_levels ))]
# d_comb_no_dupl <- d_comb_no_dupl[order(age, new_id, alter_age_strata_idx_simplot, gender, alter_gender)]

# simulated_fine_counts_MM <- ggplot(d.everything.final[gender=="Male" & alter_gender=="Male"], aes(age, factor(alter_age))) + 
#   geom_tile(aes(fill = y_age)) + 
#   scale_fill_viridis(option = "F") + 
#   scale_x_continuous(expand = c(0,0)) + 
#   scale_y_discrete(expand = c(0,0)) + 
#   labs(x = "Age of participants", y = "Age of contacts", fill = "Counts MM") + 
#   theme(aspect.ratio = 1)


# d_comb_no_dupl_plot[, cntct_rate_strata_plot := sum(cntct_rate), by=group_var_plot]




# get population dt
dt.pop <- covimod$pop


covimod.single.new.hh.sim <- list(
  contacts = d_comb_no_dupl,
  offsets = d.everything.final,
  pop = dt.pop
)

covimod.single.new.hh.sim.baseline <- list(
  contacts = d_comb_no_dupl_baseline,
  offsets = d.everything.final_baseline,
  pop = dt.pop
)

  if (args$hhsize == 0){
    if (!args$divide.Hicb){
      if (!args$drop_zero_Hicb){
      saveRDS(covimod.single.new.hh.sim, file=file.path(export.path, paste0("nodivide-data-hh0-amended.rds")))
      saveRDS(covimod.single.new.hh.sim.baseline, file=file.path(export.path, paste0("nodivide-data-hh0-amended-baseline.rds")))
      }else{
        saveRDS(covimod.single.new.hh.sim, file=file.path(export.path, paste0("nodivide-data-hh0-amended-drop-zero-Hicb.rds")))
        saveRDS(covimod.single.new.hh.sim.baseline, file=file.path(export.path, paste0("nodivide-data-hh0-amended-baseline-drop-zero-Hicb.rds")))
      }
      }
    else{
      if (!args$drop_zero_Hicb){
      saveRDS(covimod.single.new.hh.sim, file=file.path(export.path, paste0("data-hh0-amended.rds")))
      saveRDS(covimod.single.new.hh.sim.baseline, file=file.path(export.path, paste0("data-hh0-amended-baseline-drop-zero-Hicb.rds")))
      }else{
        saveRDS(covimod.single.new.hh.sim, file=file.path(export.path, paste0("data-hh0-amended.rds")))
        saveRDS(covimod.single.new.hh.sim.baseline, file=file.path(export.path, paste0("data-hh0-amended-baseline-drop-zero-Hicb.rds")))
      }
    }
  }
  if (args$hhsize > 0){
    if (!args$drop_zero_Hicb){
    saveRDS(covimod.single.new.hh.sim, file=file.path(export.path, paste0("data-hh", args$hhsize, "-", args$scenario, "-", args$size, "-amended.rds")))
    saveRDS(covimod.single.new.hh.sim.baseline, file=file.path(export.path, paste0("data-hh", args$hhsize, "-", args$scenario, "-", args$size, "-amended-baseline.rds")))
    }else{
      saveRDS(covimod.single.new.hh.sim, file=file.path(export.path, paste0("data-hh", args$hhsize, "-", args$scenario, "-", args$size, "-amended-drop-zero-Hicb.rds")))
      saveRDS(covimod.single.new.hh.sim.baseline, file=file.path(export.path, paste0("data-hh", args$hhsize, "-", args$scenario, "-", args$size, "-amended-baseline-drop-zero-Hicb.rds")))
    }
  }

cat("\n DONE. \n")
