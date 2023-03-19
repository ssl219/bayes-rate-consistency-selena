functions
{
#include gp-functions.stan
}

data
{
  int<lower=1> N_MM, N_FF, N_MF, N_FM; // Number of observations for each gender combination between the participant and the contact in strata space
  int<lower=1> A_MM, A_FF, A_MF, A_FM; // Number of observations for each gender combination between the participant and the contact in age space
  
  
  int<lower=1> P_MM, P_FF, P_FM, P_MF; // Unique number of participants with a given gender and with contacts of another gender
  int<lower=1> A;       // Number of age inputs
  int<lower=1> C;       // Number of age strata
  int<lower=1> U;       // Number of unique age and gender combinations in all observations


  array[N_MM] int Y_MM; // Contacts for all participants (ordered) to strata c
  array[N_FF] int Y_FF;
  array[N_MF] int Y_MF;
  array[N_FM] int Y_FM;

  array[A_MM] int B_MM; // flattened list of ages of contacts (ordered, ascending) in all gender combinations
  array[A_FF] int B_FF;
  array[A_MF] int B_MF;
  array[A_FM] int B_FM;

  array[P_MM + 1] int cum_MM; // list of cumulative indices for participants in B_MM
  array[P_FF + 1] int cum_FF;
  array[P_FM + 1] int cum_FM;
  array[P_MF + 1] int cum_MF;
  
  array[N_MM] int ROW_MAJOR_IDX_MM;
  array[N_FF] int ROW_MAJOR_IDX_FF;
  array[N_MF] int ROW_MAJOR_IDX_MF;
  array[N_FM] int ROW_MAJOR_IDX_FM;

  array[P_MM] int map_indiv_to_age_MM; // array of ages of each participant of gender M with contacts of gender M
  array[P_FF] int map_indiv_to_age_FF; 
  array[P_MF] int map_indiv_to_age_MF; 
  array[P_FM] int map_indiv_to_age_FM; 
  
  matrix[N_MM + N_FF + N_MF + N_FM, U] map_indiv_to_age; // map individual-age space to age-age space


  array[A_MM] real log_H_MM; // Log of household offsets (ordered in ascending order of participant age then contact age)
  array[A_FF] real log_H_FF;
  array[A_MF] real log_H_MF;
  array[A_FM] real log_H_FM;
 
  vector[A] age_idx_std;         // Standardized age index
  matrix[A,C] map_age_to_strata; // Indicator Matrix that maps age to age strata

  // HSGP parameters
  real<lower=0> C1; // Factor to determine the boundary value L (cohort age dimension)
  int<lower=1> M1;  // Number of basis functions (cohort age dimension)
  real<lower=0> C2; // Factor to determine the boundary value L for age of contacted individuals (age difference dimension)
  int<lower=1> M2;  // Number of basis functions (age difference dimension)
  



}

transformed data
{
  int N = N_MM + N_FF + N_MF + N_FM;  // Total number of observations
  int MM = 1, FF = 2, MF = 3, FM = 4;
  int G = 4; // gender combinations
  real gp_delta = 1e-9;               // GP nugget
  real epsilon = 1e-13;               // Prevent shape parameter to be 0

  real L1, L2;
  matrix[A,M1] PHI1;
  matrix[A,M2] PHI2;

  // Precompute HSGP basis functions
  L1 = C1 * max(age_idx_std);
  L2 = C2 * max(age_idx_std);
  PHI1 = PHI(A, M1, L1, age_idx_std);
  PHI2 = PHI(A, M2, L2, age_idx_std);

  // append data
  array[N] int Y = append_array( append_array( append_array(Y_MM, Y_FF), Y_MF), Y_FM);
}

parameters
{
  vector[G] beta_0;
  real<lower=0> nu;

  vector<lower=0>[G-1] gp_rho_1;
  vector<lower=0>[G-1] gp_rho_2;
  vector<lower=0, upper=pi()/2>[G-1] gp_alpha_unif;

  matrix[(G-1)*M1, M2] z; // HSGP basis function coefficients
}

transformed parameters
{
  vector<lower=0>[G-1] gp_alpha = tan(gp_alpha_unif); // Reparametrize Half-Cauchy for optimization in HMC


  matrix[A,A] f_MM, f_FF, f_MF;
  matrix<lower=0>[P_MM,A] alpha_MM;
  matrix<lower=0>[P_FF,A] alpha_FF;
  matrix<lower=0>[P_FM,A] alpha_FM;
  matrix<lower=0>[P_MF,A] alpha_MF;
  matrix<lower=0>[P_MM,C] alpha_strata_MM;
  matrix<lower=0>[P_FF,C] alpha_strata_FF;
  matrix<lower=0>[P_FM,C] alpha_strata_FM;
  matrix<lower=0>[P_MF,C] alpha_strata_MF;


  f_MM = hsgp(A, gp_alpha[MM], gp_rho_1[MM], gp_rho_2[MM],
              L1, L2, M1, M2, PHI1, PHI2, z[1:M1,]);
  f_FF = hsgp(A, gp_alpha[FF], gp_rho_1[FF], gp_rho_2[FF],
              L1, L2, M1, M2, PHI1, PHI2, z[(M1+1):2*M1,]);
  f_MF = hsgp(A, gp_alpha[MF], gp_rho_1[MF], gp_rho_2[MF],
              L1, L2, M1, M2, PHI1, PHI2, z[(2*M1+1):3*M1,]);

  // print("f_MM =", f_MM);
  // print("f_FF =", f_FF);
  // print("f_MF =", f_MF);

  alpha_MM = rep_matrix(0, P_MM, A); // initialize alpha_strata to zeros
  alpha_FF = rep_matrix(0, P_FF, A);
  alpha_MF = rep_matrix(0, P_MF, A);
  alpha_FM = rep_matrix(0, P_FM, A);

  // print("alpha_MM initialisation =", alpha_MM);
  // print("alpha_FF initialisation =", alpha_FF);
  // print("alpha_MF initialisation =", alpha_MF);
  // print("alpha_FM initialisation =", alpha_FM);


// have to add +1 to B_MM[j], cannot take index 0 for contact age 0, same thing for map_indiv_to_age_MM[i]
  for (i in 1:P_MM){
    for (j in cum_MM[i]+1:cum_MM[i+1]){
    alpha_MM[i, B_MM[j]+1]= exp(beta_0[MM] + symmetrize_from_lower_tri(f_MM)[map_indiv_to_age_MM[i]+1, B_MM[j]+1] + log_H_MM[j]);
    }
  }

  for (i in 1:P_MF){
    for (j in cum_MF[i]+1:cum_MF[i+1]){
    alpha_MF[i, B_MF[j]+1]=exp(beta_0[MF] + f_MF[map_indiv_to_age_MF[i]+1, B_MF[j]+1] + log_H_MF[j]) ;
    }
  }

  for (i in 1:P_FM){
    for (j in cum_FM[i]+1:cum_FM[i+1]){
    alpha_FM[i, B_FM[j]+1]=exp(beta_0[FM] + f_MF'[map_indiv_to_age_FM[i]+1, B_FM[j]+1] + log_H_FM[j]) ;
    }
  }

  for (i in 1:P_FF){
    for (j in cum_FF[i]+1:cum_FF[i+1]){
    alpha_FF[i, B_FF[j]+1]=exp(beta_0[FF] + symmetrize_from_lower_tri(f_FF)[map_indiv_to_age_FF[i]+1, B_FF[j]+1] + log_H_FF[j]) ;
    }
  }

  // print("alpha_MM after =", alpha_MM);
  // print("alpha_FF after =", alpha_FF);
  // print("alpha_MF after =", alpha_MF);
  // print("alpha_FM after =", alpha_FM);
  
  alpha_strata_MM = alpha_MM * map_age_to_strata / nu + epsilon;
  alpha_strata_MF = alpha_MF * map_age_to_strata / nu + epsilon;
  alpha_strata_FM = alpha_FM * map_age_to_strata / nu + epsilon;
  alpha_strata_FF = alpha_FF * map_age_to_strata / nu + epsilon;

  // print("alpha_strata_MM =", alpha_strata_MM);
  // print("alpha_strata_FF =", alpha_strata_FF);
  // print("alpha_strata_MF =", alpha_strata_MF);
  // print("alpha_strata_FM =", alpha_strata_FM);

}

model
{

  print("N_MM =", N_MM);
  print("N_FF =", N_FF);
  print("N_FM =", N_FM);
  print("N_MF =", N_MF);

  print("A_MM =", A_MM);
  print("A_FF =", A_FF);
  print("A_FM =", A_FM);
  print("A_MF =", A_MF);

  print("Y_MM =", Y_MM);
  print("Y_FF =", Y_FF);
  print("Y_FM =", Y_FM);
  print("Y_MF =", Y_MF);

  print("P_MM =", P_MM);
  print("P_FF =", P_FF);
  print("P_FM =", P_FM);
  print("P_MF =", P_MF);

  print("A =", A);
  print("C =", C);
  print("U =", U);

  print("B_MM =", B_MM);
  print("B_FF =", B_FF);
  print("B_FM =", B_FM);
  print("B_MF =", B_MF);

  print("cum_MM =", cum_MM);
  print("cum_FF =", cum_FF);
  print("cum_FM =", cum_FM);
  print("cum_MF =", cum_MF);

  print("ROW_MAJOR_IDX_MM =", ROW_MAJOR_IDX_MM);
  print("ROW_MAJOR_IDX_FF =", ROW_MAJOR_IDX_FF);
  print("ROW_MAJOR_IDX_FM =", ROW_MAJOR_IDX_FM);
  print("ROW_MAJOR_IDX_MF =", ROW_MAJOR_IDX_MF);

  print("map_indiv_to_age_MM =", map_indiv_to_age_MM);
  print("map_indiv_to_age_FF =", map_indiv_to_age_FF);
  print("map_indiv_to_age_FM =", map_indiv_to_age_FM);
  print("map_indiv_to_age_MF =", map_indiv_to_age_MF);

  print("map_indiv_to_age =", map_indiv_to_age);

  print("log_H_MM =", log_H_MM);
  print("log_H_FF =", log_H_FF);
  print("log_H_FM =", log_H_FM);
  print("log_H_MF =", log_H_MF);

  // print("age_idx_std =", age_idx_std);
  // print("map_age_to_strata =", map_age_to_strata);

  // GP priors
  target += inv_gamma_lpdf(gp_rho_1 | 5, 5);
  target += inv_gamma_lpdf(gp_rho_2 | 5, 5);
  target += cauchy_lpdf(gp_alpha | 0, 1);
  target += std_normal_lpdf( to_vector(z) );

  // Overdispersion
  target += exponential_lpdf(nu | 1);

  // baseline
  target += normal_lpdf(beta_0 | 0, 10);

  {
    vector[N] alpha_strata_flat_indiv =
    append_row(
        append_row(
          append_row(
            to_vector(alpha_strata_MM')[ROW_MAJOR_IDX_MM],
            to_vector(alpha_strata_FF')[ROW_MAJOR_IDX_FF]
          ),
          to_vector(alpha_strata_MF')[ROW_MAJOR_IDX_MF]
        ),
      to_vector(alpha_strata_FM')[ROW_MAJOR_IDX_FM]
    );
    target += neg_binomial_lpmf( Y | alpha_strata_flat_indiv, inv(nu));
  }
}

generated quantities
{
  // array[N] real log_lik;
  array[P_MM,C] int yhat_strata_MM;
  array[P_FF,C] int yhat_strata_FF;
  array[P_MF,C] int yhat_strata_MF;
  array[P_FM,C] int yhat_strata_FM;
  array[G] matrix[A,A] log_cnt_rate;
  
  log_cnt_rate[MM] = beta_0[MM] + symmetrize_from_lower_tri(f_MM);
  log_cnt_rate[FF] = beta_0[FF] + symmetrize_from_lower_tri(f_FF);
  log_cnt_rate[MF] = beta_0[MF] + f_MF;
  log_cnt_rate[FM] = beta_0[FM] + f_MF'; 


  {
    // row_vector[N] alpha_strata_flat_indiv_row = to_row_vector(alpha_strata_flat_indiv_row);
    
    // row_vector[U] alpha_strata_flat = alpha_strata_flat_indiv_row * map_indiv_to_age;
    // row_vector[N] Y_full = to_row_vector(Y);
    // row_vector[U] Y_age = Y_full * map_indiv_to_age;
    // array[U] int Y_age_array = to_array_1d(Y_age);
    // 
    // for(i in 1:U) {
    //   log_lik[i] = neg_binomial_lpmf( Y_age[i] | alpha_strata_flat[i], inv(nu) );
    // }
  }
  
  {
  //   alpha_strata_flat_full_indiv = 
  //   append_row(
  //       append_row(
  //         append_row(
  //           to_vector(alpha_strata_MM'),
  //           to_vector(alpha_strata_FF')
  //         ),
  //         to_vector(alpha_strata_MF')
  //       ),
  //     to_vector(alpha_strata_FM')
  //   );
  // 
  // alpha_strata_flat_full = alpha_strata_flat_full * map_indiv_to_age_full;
  
  // yhat_strata's in individual space
  for(i in 1:P_MM){
      yhat_strata_MM[i,:] = neg_binomial_rng( alpha_strata_MM[i,:], inv(nu) );
      }
  for(i in 1:P_FF){
      yhat_strata_FF[i,:] = neg_binomial_rng( alpha_strata_FF[i,:], inv(nu) );
      }
  for(i in 1:P_MF){
      yhat_strata_MF[i,:] = neg_binomial_rng( alpha_strata_MF[i,:], inv(nu) );
      }
  for(i in 1:P_FM){
      yhat_strata_FM[i,:] = neg_binomial_rng( alpha_strata_FM[i,:], inv(nu) );
      }
  
    }
  }

// In the words of Copilot: I have a model that I am trying to fit to some data. The model is a hierarchical model with a GP prior on the random effects. I am trying to fit the model using the NUTS sampler. The model is a bit complicated, but I have tried to simplify it as much as possible. 
  



 