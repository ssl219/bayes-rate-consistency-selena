#' Create a template stan data object
#'
#' @param A Number of participant age groups
#' @param C Number of contact age groups
#' @param T Number of waves (default NULL)
#'
#' @return list
#' @export
#'
#' @examples stan_data <- init_stan_data()
init_stan_data <- function(A = 85, C = 13, T = NULL){
  if(is.null(T)){ # Single waves
    list(A = A, C = C)
  } else { # Multiple waves
    list( T = T, U = sum(seq(1,T)), R = T, A = A, C = C )
  }
}

#' Makes a age x age strata grid
#'
#' @param A Number of participant age groups
#' @param U wave x repeat index
#' @param gender Whether to include gender columns
#'
#' @return a data.table (grid)
#'
#' @examples g <- make_grid(85, U = 3, gender = TRUE)
make_grid <- function(A, U = NULL, gender = FALSE){
  age_strata <- c("0-4","5-9","10-14","15-19","20-24","25-34","35-44",
                  "45-54","55-64","65-69","70-74","75-79","80-84")

  if(is.null(U)){ # Single wave
    if(gender){
      g <- as.data.table(expand.grid(age = seq(0,A-1),
                                     alter_age_strata = age_strata,
                                     gender = c("Male", "Female"),
                                     alter_gender = c("Male", "Female")))
    } else {

      g <- as.data.table(expand.grid(age = seq(0, A-1),
                                     alter_age_strata = age_strata))
    }
  } else { # Multiple waves
    if(gender){
      g <- as.data.table(expand.grid(u = seq(1,U),
                                     age = seq(0,A-1),
                                     alter_age_strata = age_strata,
                                     gender = c("Male", "Female"),
                                     alter_gender = c("Male", "Female")))
    } else {
      g <- as.data.table(expand.grid(u = seq(1,U),
                                     age = seq(0,A-1),
                                     alter_age_strata = age_strata))
    }
  }

  return(g)
}

add_contact_vector <- function(stan_data, contacts, single=FALSE, survey="COVIMOD", household_cnt=FALSE, new_hh=FALSE){
  if(survey == "COVIMOD"){
    
    if(!single){
      d <- contacts[order(u, age, alter_age_strata, gender, alter_gender)]

      stan_data$Y_MM <- d[gender == "Male" & alter_gender == "Male"]$y
      stan_data$Y_FF <- d[gender == "Female" & alter_gender == "Female"]$y
      stan_data$Y_MF <- d[gender == "Male" & alter_gender == "Female"]$y
      stan_data$Y_FM <- d[gender == "Female" & alter_gender == "Male"]$y

      return(stan_data)
      
      } 
    else { # Single wave
      
      # new household model that accounts for individual level households
      if(new_hh){
        d <- contacts
        covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
        # making sure order of factors in alter_age_strata is ascending instead of decreasing
        # note alter_age_strata_idx and age_strata_idx are the same, but they serve different purposes
        d[, alter_age_strata_idx:=as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
        d<- d[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
        
        # d <- contacts[order(age, new_id, alter_age_strata, gender, alter_gender)]
        
        stan_data$Y_MM <- d[gender == "Male" & alter_gender == "Male"]$y
        stan_data$Y_FF <- d[gender == "Female" & alter_gender == "Female"]$y
        stan_data$Y_MF <- d[gender == "Male" & alter_gender == "Female"]$y
        stan_data$Y_FM <- d[gender == "Female" & alter_gender == "Male"]$y
        
      }
      
      # household_cnt was the index for the old model, the model that doesn't account for individual level households
      if(household_cnt){
        d <- contacts
        covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
        # making sure order of factors in alter_age_strata is ascending instead of decreasing
        # note alter_age_strata_idx and age_strata_idx are the same, but they serve different purposes
        d[, alter_age_strata_idx:=as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
        d<- d[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
        # d <- d[order(u, age, alter_age_strata_idx, gender, alter_gender, o)]
        
        stan_data$Y_MM_1 <- d[gender == "Male" & alter_gender == "Male" & o == 1]$y
        stan_data$Y_FF_1 <- d[gender == "Female" & alter_gender == "Female" & o == 1]$y
        stan_data$Y_MF_1 <- d[gender == "Male" & alter_gender == "Female" & o == 1]$y
        stan_data$Y_FM_1 <- d[gender == "Female" & alter_gender == "Male" & o == 1]$y
        
        stan_data$Y_MM_2 <- d[gender == "Male" & alter_gender == "Male" & o == 2]$y
        stan_data$Y_FF_2 <- d[gender == "Female" & alter_gender == "Female" & o == 2]$y
        stan_data$Y_MF_2 <- d[gender == "Male" & alter_gender == "Female" & o == 2]$y
        stan_data$Y_FM_2 <- d[gender == "Female" & alter_gender == "Male" & o == 2]$y
        
        return(stan_data)
      }
      
    else{
      d <- contacts[order(age, alter_age_strata, gender, alter_gender)]
      
      stan_data$Y_MM <- d[gender == "Male" & alter_gender == "Male"]$y
      stan_data$Y_FF <- d[gender == "Female" & alter_gender == "Female"]$y
      stan_data$Y_MF <- d[gender == "Male" & alter_gender == "Female"]$y
      stan_data$Y_FM <- d[gender == "Female" & alter_gender == "Male"]$y
      
      return(stan_data)
    }  
      
    }
  }

  if(survey == "POLYMOD"){
    d <- contacts[order(age, alter_age, gender, alter_gender)]

    stan_data$Y_MM <- d[gender == "Male" & alter_gender == "Male"]$y
    stan_data$Y_FF <- d[gender == "Female" & alter_gender == "Female"]$y
    stan_data$Y_MF <- d[gender == "Male" & alter_gender == "Female"]$y
    stan_data$Y_FM <- d[gender == "Female" & alter_gender == "Male"]$y

    return(stan_data)
  }


}

# Add obs count
add_N <- function(stan_data, contacts, survey = "COVIMOD", household_cnt=FALSE, new_hh=FALSE){

  if(survey == "COVIMOD"){
    if (household_cnt){
      stan_data$N_MM_1 <- length(stan_data$Y_MM_1)
      stan_data$N_FF_1 <- length(stan_data$Y_FF_1)
      stan_data$N_MF_1 <- length(stan_data$Y_MF_1)
      stan_data$N_FM_1 <- length(stan_data$Y_FM_1)
      
      stan_data$N_MM_2 <- length(stan_data$Y_MM_2)
      stan_data$N_FF_2 <- length(stan_data$Y_FF_2)
      stan_data$N_MF_2 <- length(stan_data$Y_MF_2)
      stan_data$N_FM_2 <- length(stan_data$Y_FM_2)
    }
    
    if(new_hh){
      
      stan_data$N_MM <- length(stan_data$Y_MM)
      stan_data$N_FF <- length(stan_data$Y_FF)
      stan_data$N_MF <- length(stan_data$Y_MF)
      stan_data$N_FM <- length(stan_data$Y_FM)
      
      d <- contacts
      covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
      # making sure order of factors in alter_age_strata is ascending instead of decreasing
      # note alter_age_strata_idx and age_strata_idx are the same, but they serve different purposes
      d[, alter_age_strata_idx:=as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
      d<- d[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
      
      # d <- d[order(age, new_id, alter_age_strata, gender, alter_gender)]
      stan_data$P_MM <- length(unique(d[gender == "Male" & alter_gender == "Male"]$new_id))
      stan_data$P_FF <- length(unique(d[gender == "Female" & alter_gender == "Female"]$new_id))
      stan_data$P_MF <- length(unique(d[gender == "Male" & alter_gender == "Female"]$new_id))
      stan_data$P_FM <- length(unique(d[gender == "Female" & alter_gender == "Male"]$new_id))
      
    }
    else{
      stan_data$N_M <- length(stan_data$Y_MM)
      stan_data$N_F <- length(stan_data$Y_FF)
    }
    
  }

  if(survey == "POLYMOD"){
     stan_data$N_MM <- length(stan_data$Y_MM)
     stan_data$N_FF <- length(stan_data$Y_FF)
     stan_data$N_MF <- length(stan_data$Y_MF)
     stan_data$N_FM <- length(stan_data$Y_FM)

  }
  return(stan_data)
}

create_cumulative_list <- function(d_ordered_everything, no_participants, gender_comb="MM"){
  if (gender_comb == "MM"){
    d_ordered_everything <- d_ordered_everything[gender=="Male" & alter_gender=="Male"]
    d_ordered_everything <- d_ordered_everything[order(age, new_id, alter_age)]
  }

  if (gender_comb == "FF"){
    d_ordered_everything <- d_ordered_everything[gender=="Female" & alter_gender=="Female"]
    d_ordered_everything <- d_ordered_everything[order(age, new_id, alter_age)]
  }
  
  if (gender_comb == "MF"){
    d_ordered_everything <- d_ordered_everything[gender=="Male" & alter_gender=="Female"]
    d_ordered_everything <- d_ordered_everything[order(age, new_id, alter_age)]
  }
  
  if (gender_comb == "FM"){
    d_ordered_everything <- d_ordered_everything[gender=="Female" & alter_gender=="Male"]
    d_ordered_everything <- d_ordered_everything[order(age, new_id, alter_age)]
  }
  
  id_list = as.character(d_ordered_everything$new_id)
  duplicated_id_list = duplicated(id_list)
  cum_list = rep(0, no_participants + 1)
  j = 2 # R starts indexing at 1
  cum_list[1] = 0
  
  for (duplicated in duplicated_id_list){
    if (duplicated){
      cum_list[j-1] = cum_list[j-1] + 1
    }
    else{
      cum_list[j] = cum_list[j-1] + 1
      j = j + 1
    }
  }
  return(cum_list)
}


# for new household model and covimod survey only
# note we have to order the participants
# note we are using the dt.everything or aka the offsets dataframe for this
add_ages_contacts <- function(stan_data, everything){
  
  d <- everything[order(age, new_id, alter_age, gender, alter_gender)]
  stan_data$B_MM <- d[gender == "Male" & alter_gender == "Male"]$alter_age
  stan_data$B_FF <- d[gender == "Female" & alter_gender == "Female"]$alter_age
  stan_data$B_MF <- d[gender == "Male" & alter_gender == "Female"]$alter_age
  stan_data$B_FM <- d[gender == "Female" & alter_gender == "Male"]$alter_age
  
  stan_data$A_MM <-length(stan_data$B_MM)
  stan_data$A_FF <-length(stan_data$B_FF)
  stan_data$A_MF <-length(stan_data$B_MF)
  stan_data$A_FM <-length(stan_data$B_FM)
  
  stan_data$cum_MM <- create_cumulative_list(d, stan_data$P_MM, gender_comb="MM")
  stan_data$cum_FF <- create_cumulative_list(d, stan_data$P_FF, gender_comb="FF")
  stan_data$cum_MF <- create_cumulative_list(d, stan_data$P_MF, gender_comb="MF")
  stan_data$cum_FM <- create_cumulative_list(d, stan_data$P_FM, gender_comb="FM")
  return(stan_data)
}

add_row_major_idx <- function(stan_data, contacts, survey = "COVIMOD", household_cnt=FALSE){
  if (survey == "COVIMOD"){
    covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
    d <- contacts
    # making sure order of factors in alter_age_strata is ascending instead of decreasing
    # note alter_age_strata_idx and age_strata_idx are the same, but they serve different purposes
    d[, alter_age_strata_idx:=as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
    d <- d[order(u, age, alter_age_strata_idx, gender, alter_gender)]
    
    d[, age_idx := age + 1]
    d[, age_strata_idx := as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
    d[, row_major_idx := (age_idx-1)*stan_data$C + age_strata_idx]
    
    stan_data$ROW_MAJOR_IDX_M <- d[gender == "Male" & alter_gender == "Male"]$row_major_idx
    stan_data$ROW_MAJOR_IDX_F <- d[gender == "Female" & alter_gender == "Female"]$row_major_idx
    stan_data$ROW_MAJOR_IDX_MM <- d[gender == "Male" & alter_gender == "Male"]$row_major_idx
    stan_data$ROW_MAJOR_IDX_FF <- d[gender == "Female" & alter_gender == "Female"]$row_major_idx
    stan_data$ROW_MAJOR_IDX_MF <- d[gender == "Male" & alter_gender == "Female"]$row_major_idx
    stan_data$ROW_MAJOR_IDX_FM <- d[gender == "Female" & alter_gender == "Male"]$row_major_idx
    
    m = stan_data$ROW_MAJOR_IDX_M - 1
    f = stan_data$ROW_MAJOR_IDX_F - 1
    stan_data$ROW_MAJOR_IDX_M_AGE <- 1 + (m - m%%stan_data$C)/stan_data$C
    stan_data$ROW_MAJOR_IDX_F_AGE <- 1 + (f - f%%stan_data$C)/stan_data$C
    stan_data$ROW_MAJOR_IDX_M_STRATA <- 1 + m %% stan_data$C
    stan_data$ROW_MAJOR_IDX_F_STRATA <- 1 + f %% stan_data$C    
    return(stan_data)
  }

  if (survey == "POLYMOD"){
    d <- contacts[order(age, alter_age_strata, gender, alter_gender)]

    d[, age_idx := age + 1]
    d[, age_strata_idx := as.numeric(alter_age_strata) + 1]
    d[, row_major_idx := (age_idx-1)*85 + age_strata_idx]

    stan_data$ROW_MAJOR_IDX_MM <- d[gender == "Male" & alter_gender == "Male"]$row_major_idx
    stan_data$ROW_MAJOR_IDX_FF <- d[gender == "Female" & alter_gender == "Female"]$row_major_idx
    stan_data$ROW_MAJOR_IDX_MF <- d[gender == "Male" & alter_gender == "Female"]$row_major_idx
    stan_data$ROW_MAJOR_IDX_FM <- d[gender == "Female" & alter_gender == "Male"]$row_major_idx
  }
  
  if (survey == "simulation"){
    # this is assuming we have a single wave
    covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
    simulation_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54")
    C <- stan_data$C
    d <- contacts
    # making sure order of factors in alter_age_strata is ascending instead of decreasing
    # note alter_age_strata_idx and age_strata_idx are the same, but they serve different purposes
    d[, alter_age_strata_idx:=as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
    d<- d[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    # d <- contacts[order(age, new_id, alter_age_strata, gender, alter_gender)]
    
    d_MM <-  d[gender == "Male" & alter_gender == "Male"]
    d_FF <-  d[gender == "Female" & alter_gender == "Female"]
    d_MF <-  d[gender == "Male" & alter_gender == "Female"]
    d_FM <-  d[gender == "Female" & alter_gender == "Male"]
    
    d_MM <- d_MM[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    
    d_MM[, new_id_idx:=as.numeric(factor(new_id, levels=unique(new_id)))]
    d_MM[, age_strata_idx := as.numeric(factor(alter_age_strata, levels= covimod_strata_levels))]
    d_MM[, row_major_idx := (new_id_idx -1)*C + age_strata_idx]
    
    d_FF <- d_FF[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    d_FF[, new_id_idx:=as.numeric(factor(new_id, levels=unique(new_id)))]
    d_FF[, age_strata_idx := as.numeric(factor(alter_age_strata, levels= covimod_strata_levels))]
    d_FF[, row_major_idx := (new_id_idx -1)*C + age_strata_idx]
    
    d_MF <- d_MF[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    d_MF[, new_id_idx:=as.numeric(factor(new_id, levels=unique(new_id)))]
    d_MF[, age_strata_idx := as.numeric(factor(alter_age_strata, levels= covimod_strata_levels))]
    d_MF[, row_major_idx := (new_id_idx -1)*C + age_strata_idx]
    
    d_FM <- d_FM[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    d_FM[, new_id_idx:=as.numeric(factor(new_id, levels=unique(new_id)))]
    d_FM[, age_strata_idx := as.numeric(factor(alter_age_strata, levels= covimod_strata_levels))]
    d_FM[, row_major_idx := (new_id_idx -1)*C + age_strata_idx]
    
    stan_data$ROW_MAJOR_IDX_MM <- d_MM$row_major_idx
    stan_data$ROW_MAJOR_IDX_FF <- d_FF$row_major_idx
    stan_data$ROW_MAJOR_IDX_MF <- d_MF$row_major_idx
    stan_data$ROW_MAJOR_IDX_FM <- d_FM$row_major_idx
  }
    
  if (survey == "POLYMOD_2"){
    if (household_cnt){
      d <- contacts[order(u, age, alter_age_strata, gender, alter_gender)]
      
      d[, age_idx := age + 1]
      d[, age_strata_idx := as.numeric(alter_age_strata)]
      d[, row_major_idx := (age_idx-1)*13 + age_strata_idx]
      
      stan_data$ROW_MAJOR_IDX_MM_1 <- d[gender == "Male" & alter_gender == "Male", o == 1]$row_major_idx
      stan_data$ROW_MAJOR_IDX_FF_1 <- d[gender == "Female" & alter_gender == "Female", o == 1]$row_major_idx
      stan_data$ROW_MAJOR_IDX_MF_1 <- d[gender == "Male" & alter_gender == "Female", o == 1]$row_major_idx
      stan_data$ROW_MAJOR_IDX_FM_1 <- d[gender == "Female" & alter_gender == "Male", o == 1]$row_major_idx
      
      stan_data$ROW_MAJOR_IDX_MM_2 <- d[gender == "Male" & alter_gender == "Male", o == 2]$row_major_idx
      stan_data$ROW_MAJOR_IDX_FF_2 <- d[gender == "Female" & alter_gender == "Female", o == 2]$row_major_idx
      stan_data$ROW_MAJOR_IDX_MF_2 <- d[gender == "Male" & alter_gender == "Female", o == 2]$row_major_idx
      stan_data$ROW_MAJOR_IDX_FM_2 <- d[gender == "Female" & alter_gender == "Male", o == 2]$row_major_idx
      
    }
    else{
    # this is assuming we have a single wave
    covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
    d <- contacts
    # making sure order of factors in alter_age_strata is ascending instead of decreasing
    # note alter_age_strata_idx and age_strata_idx are the same, but they serve different purposes
    d[, alter_age_strata_idx:=as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
    d<- d[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    # d <- contacts[order(age, new_id, alter_age_strata, gender, alter_gender)]
    
    d_MM <-  d[gender == "Male" & alter_gender == "Male"]
    d_FF <-  d[gender == "Female" & alter_gender == "Female"]
    d_MF <-  d[gender == "Male" & alter_gender == "Female"]
    d_FM <-  d[gender == "Female" & alter_gender == "Male"]
    
    d_MM <- d_MM[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    
    d_MM[, new_id_idx:=as.numeric(factor(new_id, levels=unique(new_id)))]
    d_MM[, age_strata_idx := as.numeric(factor(alter_age_strata, levels= covimod_strata_levels))]
    d_MM[, row_major_idx := (new_id_idx -1)*stan_data$C + age_strata_idx]
    
    d_FF <- d_FF[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    d_FF[, new_id_idx:=as.numeric(factor(new_id, levels=unique(new_id)))]
    d_FF[, age_strata_idx := as.numeric(factor(alter_age_strata, levels= covimod_strata_levels))]
    d_FF[, row_major_idx := (new_id_idx -1)*stan_data$C + age_strata_idx]
    
    d_MF <- d_MF[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    d_MF[, new_id_idx:=as.numeric(factor(new_id, levels=unique(new_id)))]
    d_MF[, age_strata_idx := as.numeric(factor(alter_age_strata, levels= covimod_strata_levels))]
    d_MF[, row_major_idx := (new_id_idx -1)*stan_data$C + age_strata_idx]
    
    d_FM <- d_FM[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    d_FM[, new_id_idx:=as.numeric(factor(new_id, levels=unique(new_id)))]
    d_FM[, age_strata_idx := as.numeric(factor(alter_age_strata, levels= covimod_strata_levels))]
    d_FM[, row_major_idx := (new_id_idx -1)*stan_data$C + age_strata_idx]
    
    stan_data$ROW_MAJOR_IDX_MM <- d_MM$row_major_idx
    stan_data$ROW_MAJOR_IDX_FF <- d_FF$row_major_idx
    stan_data$ROW_MAJOR_IDX_MF <- d_MF$row_major_idx
    stan_data$ROW_MAJOR_IDX_FM <- d_FM$row_major_idx
    }
  }
  return(stan_data)
}

create_household_matrix <- function(d, P, A){
  H = matrix(0, nrow=P, ncol=A)
  for(i in 1:length(d$y)){
    H[d$new_id_idx[i], d$alter_age[i] + 1] = d$Hic_b[i]
  }
  return(H)
}

add_household_offsets <- function(stan_data, everything, no_log=FALSE){
  
  if(no_log){
    d <- everything
    covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
    # making sure order of factors in alter_age_strata is ascending instead of decreasing
    # note alter_age_strata_idx and age_strata_idx are the same, but they serve different purposes
    d[, alter_age_strata_idx:=as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
    d<- d[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    
    # d <- everything[order(age, new_id, alter_age_strata, gender, alter_gender)]
    
    d_MM <-  d[gender == "Male" & alter_gender == "Male"]
    d_FF <-  d[gender == "Female" & alter_gender == "Female"]
    d_MF <-  d[gender == "Male" & alter_gender == "Female"]
    d_FM <-  d[gender == "Female" & alter_gender == "Male"]
    
    covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
    
    d_MM <- d_MM[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    d_MM[, new_id_idx:=as.numeric(factor(new_id, levels=unique(new_id)))]
    
    d_FF <- d_FF[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    d_FF[, new_id_idx:=as.numeric(factor(new_id, levels=unique(new_id)))]
    
    d_MF <- d_MF[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    d_MF[, new_id_idx:=as.numeric(factor(new_id, levels=unique(new_id)))]

    d_FM <- d_FM[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
    d_FM[, new_id_idx:=as.numeric(factor(new_id, levels=unique(new_id)))]
    
    stan_data$H_MM <- create_household_matrix(d_MM, stan_data$P_MM, stan_data$A)
    stan_data$H_FF <- create_household_matrix(d_FF, stan_data$P_FF, stan_data$A)
    stan_data$H_MF <- create_household_matrix(d_MF, stan_data$P_MF, stan_data$A)
    stan_data$H_FM <- create_household_matrix(d_FM, stan_data$P_FM, stan_data$A)
    
    return(stan_data)
    
  }
  
  else{
    d_everything <- everything[order(age, new_id, alter_age, gender, alter_gender)]
    d_everything_MM <- d_everything[gender=="Male" & alter_gender=="Male"]
    d_everything_FF <- d_everything[gender=="Female" & alter_gender=="Female"]
    d_everything_MF <- d_everything[gender=="Male" & alter_gender=="Female"]
    d_everything_FM <- d_everything[gender=="Female" & alter_gender=="Male"]
    
    stan_data$log_H_MM <- log(d_everything_MM$Hic_b)
    stan_data$log_H_FF <- log(d_everything_FF$Hic_b)
    stan_data$log_H_MF <- log(d_everything_MF$Hic_b)
    stan_data$log_H_FM <- log(d_everything_FM$Hic_b)
    return(stan_data)
  }
  
}


# Add start end idx
make_start_end_idx <- function(contacts, .gender){
  d <- contacts[gender == .gender & alter_gender == .gender]

  g <- as.data.table(expand.grid(u = 1:max(contacts$u)))
  d[, idx := 1:.N]
  d <- d[, .(start_idx = min(idx), end_idx = max(idx)), by="u"]

  d <- merge(g, d, by = "u", all.x = TRUE)
  d[is.na(start_idx), start_idx := 0]
  d[is.na(end_idx), end_idx := 0]

  start_idx <- d$start_idx
  end_idx <- d$end_idx

  return(list(start_idx, end_idx))
}

add_start_end_idx <- function(stan_data, contacts){
  start_end_M <- make_start_end_idx(contacts, "Male")
  start_end_F <- make_start_end_idx(contacts, "Female")

  stan_data$START_IDX_M <- start_end_M[[1]]
  stan_data$START_IDX_F <- start_end_F[[1]]

  stan_data$END_IDX_M <- start_end_M[[2]]
  stan_data$END_IDX_F <- start_end_F[[2]]

  return(stan_data)
}

# Add participant size offsets
make_partsize_offsets <- function(offsets, .gender, U, A, survey = "COVIMOD"){
  if (survey == "COVIMOD"){
    d <- offsets[gender == .gender]
    g <- make_grid(A, U)

    d <- merge(g, d, by=c("u", "age"), all.x = TRUE)
    d[is.na(N), N := 1]
    d <- unique(d[,.(u, age, N)])

    d <- data.table::dcast(d, u ~ age, value.var = "N")[,-c("u")]

    return(log(as.matrix(d)))
  }

  if (survey == "POLYMOD"){
    d <- offsets[gender == .gender]
    g <- expand.grid(age = 0:84, alter_age = 0:84)

    d <- as.data.table( merge(g, d, by=c("age", "alter_age"), all.x = TRUE) )
    d[is.na(N), N := 1]
    d <- unique(d[, .(age, N)])

    return(log(d$N))
  }
}

add_partsize_offsets <- function(stan_data, offsets, survey = "COVIMOD"){
  U <- stan_data$U
  A <- stan_data$A

  if (survey == "COVIMOD"){
    stan_data$log_N_M <- make_partsize_offsets(offsets, "Male", U, A)
    stan_data$log_N_F <- make_partsize_offsets(offsets, "Female", U, A)
  }

  if (survey == "POLYMOD"){
    stan_data$log_N_M <- make_partsize_offsets(offsets, "Male", U, A, "POLYMOD")
    stan_data$log_N_F <- make_partsize_offsets(offsets, "Female", U, A, "POLYMOD")
  }

  return(stan_data)
}

add_part_offsets <- function(stan_data, contacts, offsets = NULL, survey = "COVIMOD", household_cnt=FALSE){
  d <- contacts

  if(survey == "COVIMOD"){
    stan_data$part_M <- d[gender == "Male" & alter_gender == "Male"]$N
    stan_data$part_F <- d[gender == "Female" & alter_gender == "Female"]$N

    stan_data$S_M <- d[gender == "Male" & alter_gender == "Male"]$zeta
    stan_data$S_F <- d[gender == "Female" & alter_gender == "Female"]$zeta
  }

  if(survey == "POLYMOD"){
    
    if(household_cnt){
      d1 <- complete(offsets[gender == "Male"], age = 0:84, fill = list(N = 1, zeta = 1))
      d2 <- complete(offsets[gender == "Female"], age = 0:84, fill = list(N = 1, zeta = 1))
      
      stan_data$log_N_M <- log( d1$N )
      stan_data$log_N_F <- log( d2$N )
      
      stan_data$log_S_M_1 <- log( d1[o == 1]$zeta )
      stan_data$log_S_F_1 <- log( d2[o == 1]$zeta )
      stan_data$log_S_M_2 <- log( d1[o == 2]$zeta )
      stan_data$log_S_F_2 <- log( d2[o == 2]$zeta )
    }
    else{
      # need the complete() function as we are missing age 83
      max_age <- stan_data$A - 1
      d1 <- complete(offsets[gender == "Male"], age = 0:max_age, fill = list(N = 1, zeta = 1))
      d2 <- complete(offsets[gender == "Female"], age = 0:max_age, fill = list(N = 1, zeta = 1))
      
      stan_data$log_N_M <- log( d1$N )
      stan_data$log_N_F <- log( d2$N )
      
      stan_data$log_S_M <- log( d1$zeta )
      stan_data$log_S_F <- log( d2$zeta )
    }
  }

  return(stan_data)
}

# Add zeta
make_zeta_offsets <- function(offsets, .gender, U, A){
  d <- offsets[gender == .gender]
  g <- make_grid(U, A)

  d <- merge(g, d, by=c("u", "age"), all.x=T)
  d[is.na(N), zeta := 1]
  d <- unique(d[,.(u, age, zeta)])
  d <- data.table::dcast(d, u ~ age, value.var = "zeta")[,-c("u")]

  return(log(as.matrix(d)))
}

add_zeta_offsets <- function(stan_data, offsets){
  U <- stan_data$U
  A <- stan_data$A

  stan_data$log_S_M <- make_zeta_offsets(offsets, "Male", U, A)
  stan_data$log_S_F <- make_zeta_offsets(offsets, "Female", U, A)

  return(stan_data)
}

# Add population offsets
add_pop_offsets <- function(stan_data, pop, survey = "COVIMOD"){
  A <- stan_data$A

  if(survey == "COVIMOD"){
    stan_data$pop_M <- dt.pop[gender == "Male" & age < A]$pop
    stan_data$pop_F <- dt.pop[gender == "Female" & age < A]$pop
  }

  if(survey == "POLYMOD"){
    stan_data$log_P_M <- log( dt.pop[gender == "Male" & age < A]$pop )
    stan_data$log_P_F <- log( dt.pop[gender == "Female" & age < A]$pop )
  }

  return(stan_data)
}

add_missing_u <- function(stan_data, u.missing){
  U <- stan_data$U

  missing_u <- rep(0, U)
  for(u in 1:U){
    if(u %in% u.missing) missing_u[u] <- 1
  }

  stan_data$missing_u <- missing_u # I'm missing U (haha)
  return(stan_data)
}

# Map U to T
add_map_u_to_t <- function(stan_data){
  U <- stan_data$U
  T <- stan_data$T

  map_u_to_t <- c()
  for (t in 1:T) map_u_to_t <- c(map_u_to_t, rep(t, t))

  stan_data$map_u_to_t <- map_u_to_t
  return(stan_data)
}

# Map U to R
add_map_u_to_r <- function(stan_data){
  T <- stan_data$T
  R <- stan_data$R

  map_u_to_r <- c()
  for (t in 1:T) map_u_to_r <- c(map_u_to_r, seq(1, t))

  stan_data$map_u_to_r <- map_u_to_r
  return(stan_data)
}

# Add indicator matrix for t
add_delta_t <- function(stan_data){
  T <- stan_data$T
  delta_t <- diag(1, T-1, T-1)
  delta_t <- rbind(rep(0,T-1), delta_t)

  stan_data$delta_t <- delta_t
  return(stan_data)
}

# Add indicator matrix for r
add_delta_r <- function(stan_data){
  R <- stan_data$R

  delta_r <- diag(1, R-1, R-1)
  delta_r <- rbind(rep(0,R-1), delta_r)

  stan_data$delta_r <- delta_r
  return(stan_data)
}

add_map_tr_to_u <- function(stan_data){
  U <- stan_data$U
  T <- stan_data$T
  R <- stan_data$R

  lower_tri_idx <- rep(0, (T*(T+1))%/%2)
  tmp_int = 1
  for (t in 1:T){
    for (r in 1:t){
      lower_tri_idx[tmp_int] = T*(t-1) + r;
      tmp_int = tmp_int + 1;
    }
  }

  tmp <- rep(0, T*R)
  tmp[lower_tri_idx] <- seq(1:U)

  map_tr_to_u <- matrix(tmp, nrow=T, ncol=R, byrow = T)

  stan_data$map_tr_to_u <- map_tr_to_u

  return(stan_data)
}

# Map array of ages of each participant of gender M with contacts of gender M

add_map_indiv_to_age <- function(stan_data, contact, everything, sim=FALSE){
  
  d <- everything[order(age, new_id, alter_age, gender, alter_gender)]
  d_MM <- d[gender=="Male" & alter_gender=="Male"]
  d_FF <- d[gender=="Female" & alter_gender=="Female"]
  d_MF <- d[gender=="Male" & alter_gender=="Female"]
  d_FM <- d[gender=="Female" & alter_gender=="Male"]
  
  d_MM_unique_ages <- d_MM %>% distinct(new_id, .keep_all=TRUE)
  d_FF_unique_ages <- d_FF %>% distinct(new_id, .keep_all=TRUE)
  d_MF_unique_ages <- d_MF %>% distinct(new_id, .keep_all=TRUE)
  d_FM_unique_ages <- d_FM %>% distinct(new_id, .keep_all=TRUE)
  
  stan_data$map_indiv_to_age_MM <- d_MM_unique_ages$age
  stan_data$map_indiv_to_age_FF <- d_FF_unique_ages$age
  stan_data$map_indiv_to_age_FM <- d_FM_unique_ages$age
  stan_data$map_indiv_to_age_MF <- d_MF_unique_ages$age
  # finding general map individual to age 
  tmp <- contact
  covimod_strata_levels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-34", "35-44", "45-54", "55-64", "65-69", "70-74", "75-79", "80-84")
  # making sure order of factors in alter_age_strata is ascending instead of decreasing
  # note alter_age_strata_idx and age_strata_idx are the same, but they serve different purposes
  tmp[, alter_age_strata_idx:=as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
  tmp<- tmp[order(age, new_id, alter_age_strata_idx, gender, alter_gender)]
  # tmp <- tmp[order(age, alter_age_strata, gender, alter_gender)]
  tmp[, age_idx := age + 1]
  tmp[, age_strata_idx := as.numeric(factor(alter_age_strata, levels=covimod_strata_levels))]
  tmp[, row_major_idx := (age_idx-1)*stan_data$C + age_strata_idx]
  tmp[, gender_comb_idx := fcase(
    gender == "Male" & alter_gender == "Male", 1,
    gender == "Female" & alter_gender == "Female", 2,
    gender == "Male" & alter_gender == "Female", 3,
    gender == "Female" & alter_gender == "Male", 4,
    default = NA
  )]
  tmp[, f:=row_major_idx*4 + gender_comb_idx]
  
  f_list = tmp$f
  unique_f_list = unique(f_list)
  nb_f = length(unique_f_list)
  N = length(tmp$y)
  
  map_indiv_to_age = matrix(0, nrow=N, ncol=nb_f)
  # print(map_indiv_to_age)
  # creating map individual to age matrix
  
  for (i in 1:N){
    for (w in which(unique_f_list==f_list[i])){
      map_indiv_to_age[i, w] <- 1
    }
  }
  
  stan_data$U <- nb_f
  # stan_data$map_indiv_to_age <- map_indiv_to_age
  
  return(stan_data)

}
# Map age to age strata
add_map_age_to_strata <- function(stan_data, survey = "COVIMOD"){
  if (survey == "COVIMOD"){
    age_strata <- c("0-4","5-9","10-14","15-19","20-24","25-34","35-44",
                    "45-54","55-64","65-69","70-74","75-79","80-84")

    alter_age_min <- as.numeric(str_extract(as.character(age_strata), "^[0-9]{1,2}"))

    strata_min <- unique(alter_age_min)
    strata_min_idx <- strata_min - min(strata_min) + 1

    A <- stan_data$A
    C <- stan_data$C

    map_age_to_strata <- matrix(0, nrow=A, ncol=C)
    for (c in 1:C) {
      if (c == C) {
        map_age_to_strata[strata_min_idx[c]:A, c] <- 1
      } else {
        map_age_to_strata[strata_min_idx[c]:(strata_min_idx[c+1]-1), c] <- 1
      }
    }

    stan_data$map_age_to_strata <- map_age_to_strata

    return(stan_data)
  }
  
  if (survey == "simulation"){
    age_strata <- c("0-4","5-9","10-14","15-19","20-24","25-34","35-44",
                    "45-54")
    
    alter_age_min <- as.numeric(str_extract(as.character(age_strata), "^[0-9]{1,2}"))
    
    strata_min <- unique(alter_age_min)
    strata_min_idx <- strata_min - min(strata_min) + 1
    
    A <- stan_data$A
    C <- stan_data$C
    
    map_age_to_strata <- matrix(0, nrow=A, ncol=C)
    for (c in 1:C) {
      if (c == C) {
        map_age_to_strata[strata_min_idx[c]:A, c] <- 1
      } else {
        map_age_to_strata[strata_min_idx[c]:(strata_min_idx[c+1]-1), c] <- 1
      }
    }
    
    stan_data$map_age_to_strata <- map_age_to_strata
    
    return(stan_data)
  }

  if (survey == "POLYMOD"){
    map_age_to_strata <- diag(rep(1, stan_data$A))

    stan_data$map_age_to_strata <- map_age_to_strata

    return(stan_data)
  }
}

reciprocal_matrix <- function(mat){
  for (i in 1:nrow(mat)){
    for (j in 1:ncol(mat)){
      if (mat[i, j] !=0){
        mat[i, j] = 1/mat[i, j]
      }
    }
  }
  return(mat)
}
# Add beta weights
add_beta_weights <- function(stan_data){
  
  # extract values needed
  A <- stan_data$A
  P_MM <- stan_data$P_MM
  P_FF <- stan_data$P_FF
  P_MF <- stan_data$P_MF
  P_FM <- stan_data$P_FM
  B_MM <- stan_data$B_MM
  B_FF <- stan_data$B_FF
  B_MF <- stan_data$B_MF
  B_FM <- stan_data$B_FM
  cum_MM <- stan_data$cum_MM
  cum_FF <- stan_data$cum_FF
  cum_MF <- stan_data$cum_MF
  cum_FM <- stan_data$cum_FM
  map_indiv_to_age_MM <- stan_data$map_indiv_to_age_MM
  map_indiv_to_age_FF <- stan_data$map_indiv_to_age_FF
  map_indiv_to_age_MF <- stan_data$map_indiv_to_age_MF
  map_indiv_to_age_FM <- stan_data$map_indiv_to_age_FM
  
  # initialise
  beta_weights_MM = matrix(0, A, A)
  beta_weights_FF = matrix(0, A, A)
  beta_weights_MF = matrix(0, A, A)
  beta_weights_FM = matrix(0, A, A)
  
  
  for (i in 1:P_MM){
    for (j in cum_MM[i]+1:cum_MM[i+1]){
      beta_weights_MM[map_indiv_to_age_MM[i]+1, B_MM[j]+1] = beta_weights_MM[map_indiv_to_age_MM[i]+1, B_MM[j]+1] + 1;
    }
  }
  
  for (i in 1:P_FF){
    for (j in cum_FF[i]+1:cum_FF[i+1]){
      beta_weights_FF[map_indiv_to_age_FF[i]+1, B_FF[j]+1] = beta_weights_FF[map_indiv_to_age_FF[i]+1, B_FF[j]+1] + 1;
    }
  }
  
  for (i in 1:P_MF){
    for (j in cum_MF[i]+1:cum_MF[i+1]){
      beta_weights_MF[map_indiv_to_age_MF[i]+1, B_MF[j]+1] =  beta_weights_MF[map_indiv_to_age_MF[i]+1, B_MF[j]+1] + 1;
    }
  }
  
  for (i in 1:P_FM){
    for (j in cum_FM[i]+1:cum_FM[i+1]){
      beta_weights_FM[map_indiv_to_age_FM[i]+1, B_FM[j]+1] = beta_weights_FM[map_indiv_to_age_FM[i]+1, B_FM[j]+1] + 1;
    }
  }
  
  stan_data$beta_weights_MM <- reciprocal_matrix(beta_weights_MM)
  stan_data$beta_weights_FF <- reciprocal_matrix(beta_weights_FF)
  stan_data$beta_weights_MF <- reciprocal_matrix(beta_weights_MF)
  stan_data$beta_weights_FM <- reciprocal_matrix(beta_weights_FM)
  
  return (stan_data)
}

# Add non-nuisance index
add_nn_idx <- function(stan_data){
  A <- stan_data$A

  NN_IDX <- rep(NA, A*A)
  for (i in 1:A){
    min_idx <- (2*A-1)*(i-1) + (A-i+1)
    max_idx <- (2*A-1)*(i-1) + (2*A-i)
    NN_IDX[A*(i-1)+1:A*1] <- seq.int(min_idx, max_idx, 1)
  }

  stan_data$NN_IDX <- NN_IDX
  return(stan_data)
}

# Add standardized age index
add_std_age_idx <- function(stan_data){
  A <- stan_data$A
  C <- stan_data$C
  
  age_idx <- seq.int(0,A-1,1)
  diff_idx <- seq.int(-A+1, A-1, 1)
  strata_idx <- seq.int(0, C-1, 1)
  
  age_idx_std <- (age_idx - mean(age_idx))/sd(age_idx)
  diff_idx_std <- (diff_idx - mean(diff_idx))/sd(diff_idx)
  strata_idx_std <- (strata_idx - mean(strata_idx))/sd(strata_idx)
  
  stan_data$age_idx_std <- age_idx_std
  stan_data$diff_idx_std <- diff_idx_std
  stan_data$strata_idx_std <- strata_idx_std
  return(stan_data)
}


# Add HSGP parameters
add_hsgp_parms <- function(stan_data, C, M1, M2){
  stan_data$C1 = C
  stan_data$C2 = C
  stan_data$M1 = M1
  stan_data$M2 = M2

  return(stan_data)
}

