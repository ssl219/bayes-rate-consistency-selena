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
args$scenario = "flat"
args$hhsize = 0
args$divide.Hicb = FALSE

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

##### ---------- Simulating contact survey data ---------- ##########
cat("\n Generating contact dataset ...")

# Start with male participants
args$size = 0
args$strata = "COVIMOD-new-hh"
N_random <- args$size
N = 2*(85 + N_random)
# Generate new_id for all participants (all distinct)
d.everything.M <- as.data.table(expand.grid(new_id = 1:N, alter_age = 0:84))
d.everything.M[, DUMMY:=1L]
tmp <- data.table( DUMMY = 1L, alter_gender = c('Male','Female'))
d.everything.M <- merge(d.everything.M, tmp, by = 'DUMMY', allow.cartesian = TRUE)
d.everything.M <- stratify_alter_age(d.everything.M, args$strata)
d.everything.M[, wave:=1L]
# hhsize = args$hhsize
d.everything.M[, household_size:=4L]
# only male participants for now
d.everything.M[, gender:="Male"]

# Generate age of participants
# set.seed(12002)
random_ages <- sample(0:84, N_random, replace=TRUE)
part_ages <- sort(c(random_ages, 0:84))
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

if (args$hhsize==4){
  for (i in 1:N_MM){
    d.everything.MM[, n_M:=fcase(
      new_id==i & age %in% 0:20 & alter_age %in% 0:20, as.numeric(sample(c(0, 1), size=1, prob=c(0.5, 0.5))),
      new_id==i & age %in% 0:20 & alter_age %in% 35:44, as.numeric(sample(c(0, 1, 2), size=1, prob=c(0.05, 0.9, 0.05))),
      new_id==i & age %in% 20:34 & alter_age %in% 20:34, as.numeric(sample(c(0, 1, 2, 3), size=1, prob=c(0.5, 0.5, 0.5, 0.5))),
      new_id==i & age %in% 35:44 & alter_age %in% 0:20, as.numeric(sample(c(0, 1, 2), size=1, prob=c(0.5, 0.5, 0.5))),
      new_id==i & age %in% 35:44 & alter_age %in% 35:44, as.numeric(sample(c(0, 1), size=1, prob=c(0.9, 0.1))),
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
    alter_age_strata == "45-54", as.numeric(n_M/10),
    alter_age_strata == "55-64", as.numeric(n_M/10),
    alter_age_strata == "65-69", as.numeric(n_M/5),
    alter_age_strata == "70-74", as.numeric(n_M/5),
    alter_age_strata == "75-79", as.numeric(n_M/5),
    alter_age_strata == "80-84", as.numeric(n_M/5),
    alter_age_strata =="85+", as.numeric(n_M),
    default = NA
  )]
  
  d.everything.MF[, Hic_b := fcase(
    alter_age_strata == "0-4", as.numeric(n_F/5),
    alter_age_strata == "5-9", as.numeric(n_F/5),
    alter_age_strata == "10-14", as.numeric(n_F/5),
    alter_age_strata == "15-19", as.numeric(n_F/5),
    alter_age_strata == "20-24", as.numeric(n_F/5),
    alter_age_strata == "25-34", as.numeric(n_F/10),
    alter_age_strata == "35-44", as.numeric(n_F/10),
    alter_age_strata == "45-54", as.numeric(n_F/10),
    alter_age_strata == "55-64", as.numeric(n_F/10),
    alter_age_strata == "65-69", as.numeric(n_F/5),
    alter_age_strata == "70-74", as.numeric(n_F/5),
    alter_age_strata == "75-79", as.numeric(n_F/5),
    alter_age_strata == "80-84", as.numeric(n_F/5),
    alter_age_strata =="85+", as.numeric(n_F),
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
    new_id==i & age %in% 0:20 & alter_age %in% 0:20, as.numeric(1-n_M),
    new_id==i & age %in% 0:20 & alter_age %in% 35:44, as.numeric(2-n_M),
    new_id==i & age %in% 20:34 & alter_age %in% 20:34, as.numeric(3-n_M),
    new_id==i & age %in% 35:44 & alter_age %in% 0:20, as.numeric(2-n_M),
    new_id==i & age %in% 35:44 & alter_age %in% 35:44, as.numeric(1-n_M),
    new_id < i, as.numeric(n_F),
    default = 0
  )]
}  

# divide by size of interval to get Hic_b
d.everything.MM[, Hic_b:=fcase(
  age %in% 0:20 & alter_age %in% 0:20, as.numeric(n_M/21),
  age %in% 0:20 & alter_age %in% 35:44, as.numeric(n_M/10),
  age %in% 20:34 & alter_age %in% 20:34, as.numeric(n_M/15),
  age %in% 35:44 & alter_age %in% 0:20, as.numeric(n_M/21),
  age %in% 35:44 & alter_age %in% 35:44, as.numeric(n_M/10),
  default = 0
)]

d.everything.MF[, Hic_b:=fcase(
  age %in% 0:20 & alter_age %in% 0:20, as.numeric(n_F/21),
  age %in% 0:20 & alter_age %in% 35:44, as.numeric(n_F/10),
  age %in% 20:34 & alter_age %in% 20:34, as.numeric(n_F/15),
  age %in% 35:44 & alter_age %in% 0:20, as.numeric(n_F/21),
  age %in% 35:44 & alter_age %in% 35:44, as.numeric(n_F/10),
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
d.everything.M.final <- d.everything.M.final[Hic_b!=0]




# do the same for female participants
args$size = 85
args$strata = "COVIMOD-new-hh"
N_random <- args$size
N = 2*(85 + N_random)
# Generate new_id for all participants (all distinct)
d.everything.F <- as.data.table(expand.grid(new_id = (N+1):(2*N), alter_age = 0:84))
d.everything.F[, DUMMY:=1L]
tmp <- data.table( DUMMY = 1L, alter_gender = c('Male','Female'))
d.everything.F <- merge(d.everything.F, tmp, by = 'DUMMY', allow.cartesian = TRUE)
d.everything.F <- stratify_alter_age(d.everything.F, args$strata)
d.everything.F[, wave:=1L]
# hhsize = args$hhsize
d.everything.F[, household_size:=4L]
# female participants now
d.everything.F[, gender:="Female"]

# Generate age of participants
set.seed(2002)
random_ages <- sample(0:84, N_random, replace=TRUE)
part_ages <- sort(c(random_ages, 0:84))
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
      new_id==i & age %in% 0:20 & alter_age %in% 0:20, as.numeric(sample(c(0, 1), size=1, prob=c(0.5, 0.5))),
      new_id==i & age %in% 0:20 & alter_age %in% 35:44, as.numeric(sample(c(0, 1, 2), size=1, prob=c(0.05, 0.9, 0.05))),
      new_id==i & age %in% 20:34 & alter_age %in% 20:34, as.numeric(sample(c(0, 1, 2, 3), size=1, prob=c(0.5, 0.5, 0.5, 0.5))),
      new_id==i & age %in% 35:44 & alter_age %in% 0:20, as.numeric(sample(c(0, 1, 2), size=1, prob=c(0.5, 0.5, 0.5))),
      new_id==i & age %in% 35:44 & alter_age %in% 35:44, as.numeric(sample(c(0, 1), size=1, prob=c(0.9, 0.1))),
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
      alter_age_strata == "45-54", as.numeric(n_M/10),
      alter_age_strata == "55-64", as.numeric(n_M/10),
      alter_age_strata == "65-69", as.numeric(n_M/5),
      alter_age_strata == "70-74", as.numeric(n_M/5),
      alter_age_strata == "75-79", as.numeric(n_M/5),
      alter_age_strata == "80-84", as.numeric(n_M/5),
      alter_age_strata =="85+", as.numeric(n_M),
      default = NA
    )]
    
    d.everything.FF[, Hic_b:=fcase(
      alter_age_strata == "0-4", as.numeric(n_F/5),
      alter_age_strata == "5-9", as.numeric(n_F/5),
      alter_age_strata == "10-14", as.numeric(n_F/5),
      alter_age_strata == "15-19", as.numeric(n_F/5),
      alter_age_strata == "20-24", as.numeric(n_F/5),
      alter_age_strata == "25-34", as.numeric(n_F/10),
      alter_age_strata == "35-44", as.numeric(n_F/10),
      alter_age_strata == "45-54", as.numeric(n_F/10),
      alter_age_strata == "55-64", as.numeric(n_F/10),
      alter_age_strata == "65-69", as.numeric(n_F/5),
      alter_age_strata == "70-74", as.numeric(n_F/5),
      alter_age_strata == "75-79", as.numeric(n_F/5),
      alter_age_strata == "80-84", as.numeric(n_F/5),
      alter_age_strata =="85+", as.numeric(n_F),
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
    new_id==i & age %in% 0:20 & alter_age %in% 0:20, as.numeric(1-n_M),
    new_id==i & age %in% 0:20 & alter_age %in% 35:44, as.numeric(2-n_M),
    new_id==i & age %in% 20:34 & alter_age %in% 20:34, as.numeric(3-n_M),
    new_id==i & age %in% 35:44 & alter_age %in% 0:20, as.numeric(2-n_M),
    new_id==i & age %in% 35:44 & alter_age %in% 35:44, as.numeric(1-n_M),
    new_id < i, as.numeric(n_F),
    default = 0
  )]
}  

# divide by size of interval to get Hic_b
d.everything.FM[, Hic_b:=fcase(
  age %in% 0:20 & alter_age %in% 0:20, as.numeric(n_M/21),
  age %in% 0:20 & alter_age %in% 35:44, as.numeric(n_M/10),
  age %in% 20:34 & alter_age %in% 20:34, as.numeric(n_M/15),
  age %in% 35:44 & alter_age %in% 0:20, as.numeric(n_M/21),
  age %in% 35:44 & alter_age %in% 35:44, as.numeric(n_M/10),
  default = 0
)]

d.everything.FF[, Hic_b:=fcase(
  age %in% 0:20 & alter_age %in% 0:20, as.numeric(n_F/21),
  age %in% 0:20 & alter_age %in% 35:44, as.numeric(n_F/10),
  age %in% 20:34 & alter_age %in% 20:34, as.numeric(n_F/15),
  age %in% 35:44 & alter_age %in% 0:20, as.numeric(n_F/21),
  age %in% 35:44 & alter_age %in% 35:44, as.numeric(n_F/10),
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
d.everything.F.final <- d.everything.F.final[Hic_b!=0]


# bind female and male participants
d.everything.final <- rbind(d.everything.M.final, d.everything.F.final)


# calculate means alpha
d.everything.final <- merge(d.everything.final, dt, by=c("alter_age","age", "gender", "alter_gender", "alter_age_strata"), all.x=TRUE, all.y=FALSE)
d.everything.final[, alpha:=Hic_b*cntct_rate]

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
group_var <- c("age", "gender", "alter_age_strata", "alter_gender")
d_comb_no_dupl <- d.everything.final
d_comb_no_dupl[, y_strata := mean(y), by=group_var]
d_comb_no_dupl[, cntct_rate_strata := mean(cntct_rate), by=group_var]
setnames(d_comb_no_dupl, c("y", "y_strata"), c("y_age", "y"))
# order alter_age strate for plots
covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
d_comb_no_dupl[, alter_age_strata_idx_simplot:=as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
d_comb_no_dupl <- d_comb_no_dupl[order(age, new_id, alter_age_strata_idx_simplot, gender, alter_gender)]

simulated_counts <- ggplot(d_comb_no_dupl, aes(age, factor(alter_age_strata, levels=covimod_strata_levels))) + 
       geom_tile(aes(fill = y)) + 
       scale_fill_viridis(option = "F") + 
       scale_x_continuous(expand = c(0,0)) + 
       scale_y_discrete(expand = c(0,0)) + 
       labs(x = "Age of participants", y = "Age of contacts", fill = "Counts") + 
       theme(aspect.ratio = 1)

simulated_counts_MM <- ggplot(d_comb_no_dupl[gender=="Male" & alter_gender=="Female"], aes(age, factor(alter_age_strata, levels=covimod_strata_levels))) + 
  geom_tile(aes(fill = y)) + 
  scale_fill_viridis(option = "F") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  labs(x = "Age of participants", y = "Age of contacts", fill = "Counts MM") + 
  theme(aspect.ratio = 1)

simulated_counts_FF <- ggplot(d_comb_no_dupl[gender=="Female" & alter_gender=="Female"], aes(age, factor(alter_age_strata, levels=covimod_strata_levels))) + 
  geom_tile(aes(fill = y)) + 
  scale_fill_viridis(option = "F") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  labs(x = "Age of participants", y = "Age of contacts", fill = "Counts FF") + 
  theme(aspect.ratio = 1)

simulated_counts_MF <- ggplot(d_comb_no_dupl[gender=="Male" & alter_gender=="Female"], aes(age, factor(alter_age_strata, levels=covimod_strata_levels))) + 
  geom_tile(aes(fill = y)) + 
  scale_fill_viridis(option = "F") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  labs(x = "Age of participants", y = "Age of contacts", fill = "Counts MF") + 
  theme(aspect.ratio = 1)


simulated_counts_FM <- ggplot(d_comb_no_dupl[gender=="Female" & alter_gender=="Male"], aes(age, factor(alter_age_strata, levels=covimod_strata_levels))) + 
  geom_tile(aes(fill = y)) + 
  scale_fill_viridis(option = "F") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  labs(x = "Age of participants", y = "Age of contacts", fill = "Counts FM") + 
  theme(aspect.ratio = 1)

if (!args$divide.Hicb){
  ggsave(file.path(export.path, paste0("hh", args$hhsize, "-nodivide-simulated-counts.pdf")), plot = simulated_counts, width = 10, height = 6)
  ggsave(file.path(export.path, paste0("hh", args$hhsize,"-nodivide-simulated-counts-MM.pdf")), plot = simulated_counts_MM, width = 10, height = 6)
  ggsave(file.path(export.path, paste0("hh", args$hhsize,"-nodivide-simulated-counts-FF.pdf")), plot = simulated_counts_FF, width = 10, height = 6)
  ggsave(file.path(export.path, paste0("hh", args$hhsize,"-nodivide-simulated-counts-MF.pdf")), plot = simulated_counts_MF, width = 10, height = 6)
  ggsave(file.path(export.path, paste0("hh", args$hhsize,"-nodivide-simulated-counts-FM.pdf")), plot = simulated_counts_FM, width = 10, height = 6)
  }else{
  ggsave(file.path(export.path, paste0("hh", args$hhsize, "-simulated-counts.pdf")), plot = simulated_counts, width = 10, height = 6)
  ggsave(file.path(export.path, paste0("hh", args$hhsize,"-simulated-counts-MM.pdf")), plot = simulated_counts_MM, width = 10, height = 6)
  ggsave(file.path(export.path, paste0("hh", args$hhsize,"-simulated-counts-FF.pdf")), plot = simulated_counts_FF, width = 10, height = 6)
  ggsave(file.path(export.path, paste0("hh", args$hhsize,"-simulated-counts-MF.pdf")), plot = simulated_counts_MF, width = 10, height = 6)
  ggsave(file.path(export.path, paste0("hh", args$hhsize,"-simulated-counts-FM.pdf")), plot = simulated_counts_FM, width = 10, height = 6)
  }

covimod.single.new.hh.sim <- list(
  contacts = d_comb_no_dupl,
  offsets = d.everything.final,
  pop = d.everything.final
)

if (args$hhsize == 0){
  if (!args$divide.Hicb){
    saveRDS(covimod.single.new.hh.sim, file=file.path(export.path, paste0("nodivide-data-hh0.rds")))
  }
  else{
    saveRDS(covimod.single.new.hh.sim, file=file.path(export.path, paste0("data-hh0.rds")))
  }
  
}

if (args$hhsize == 4){
  if (!args$divide.Hicb){
    saveRDS(covimod.single.new.hh.sim, file=file.path(export.path, paste0("data-hh4.rds")))
  }
  else{
    saveRDS(covimod.single.new.hh.sim, file=file.path(export.path, paste0("-nodivide-data-hh4.rds")))
  }
}

cat("\n DONE. \n")
