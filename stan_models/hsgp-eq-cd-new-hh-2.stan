functions
{
#include gp-functions.stan
}

data
{
  int<lower=1> N_MM, N_FF, N_MF, N_FM; // Number of observations for each gender combination between the participant and the contact
  
  int<lower=1> P_MM, P_FF, P_FM, P_MF; // Unique number of participants with a given gender and with contacts of another gender
  int<lower=1> A;       // Number of age inputs
  int<lower=1> C;       // Number of age strata
  int<lower=1> U;       // Number of unique age and gender combinations in all observations

  array[N_MM] int Y_MM; // Contacts for all participants (ordered) to strata c
  array[N_FF] int Y_FF;
  array[N_MF] int Y_MF;
  array[N_FM] int Y_FM;
  
  array[P_MM] int map_indiv_to_age_MM; // array of ages of each participant of gender M with contacts of gender M
  array[P_FF] int map_indiv_to_age_FF; 
  array[P_FM] int map_indiv_to_age_FM; 
  array[P_MF] int map_indiv_to_age_MF; 
  
  matrix[N_MM + N_FF + N_MF + N_FM, U] map_indiv_to_age; // map individual-age space to age-age space

  matrix[P_MM, A] H_MM; // Household size offsets
  matrix[P_FF, A] H_FF;
  matrix[P_FM, A] H_FM;
  matrix[P_MF, A] H_MF;
  
  array[N_MM] int ROW_MAJOR_IDX_MM;
  array[N_FF] int ROW_MAJOR_IDX_FF;
  array[N_MF] int ROW_MAJOR_IDX_MF;
  array[N_FM] int ROW_MAJOR_IDX_FM;

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

  vector<lower=0>[G] gp_rho_1;
  vector<lower=0>[G] gp_rho_2;
  vector<lower=0, upper=pi()/2>[G] gp_alpha_unif;

  matrix[(G)*M1, M2] z; // HSGP basis function coefficients
}

transformed parameters
{
  vector<lower=0>[G] gp_alpha = tan(gp_alpha_unif); // Reparametrize Half-Cauchy for optimization in HMC

  matrix[A,A] f_MM, f_FF, f_MF, f_FM;
  matrix[P_MM,A] log_m_MM;
  matrix[P_FF,A] log_m_FF;
  matrix[P_FM,A] log_m_FM;
  matrix[P_MF,A] log_m_MF;
  matrix<lower=0>[P_MM,C] alpha_strata_MM;
  matrix<lower=0>[P_FF,C] alpha_strata_FF;
  matrix<lower=0>[P_FM,C] alpha_strata_FM;
  matrix<lower=0>[P_MF,C] alpha_strata_MF;
  matrix[P_MM, A] part_f_MM;
  matrix[P_MF, A] part_f_MF;
  matrix[P_FF, A] part_f_FF;
  matrix[P_FM, A] part_f_FM;

  f_MM = hsgp(A, gp_alpha[MM], gp_rho_1[MM], gp_rho_2[MM],
              L1, L2, M1, M2, PHI1, PHI2, z[1:M1,]);
  f_FF = hsgp(A, gp_alpha[FF], gp_rho_1[FF], gp_rho_2[FF],
              L1, L2, M1, M2, PHI1, PHI2, z[(M1+1):2*M1,]);
  f_MF = hsgp(A, gp_alpha[MF], gp_rho_1[MF], gp_rho_2[MF],
              L1, L2, M1, M2, PHI1, PHI2, z[(2*M1+1):3*M1,]);
  f_FM = hsgp(A, gp_alpha[FM], gp_rho_1[FM], gp_rho_2[FM],
              L1, L2, M1, M2, PHI1, PHI2, z[(3*M1+1):4*M1,]);
  
  // initialise          
  part_f_MM = rep_matrix(0, P_MM, A);
  part_f_FF = rep_matrix(0, P_FF, A);
  part_f_FM = rep_matrix(0, P_FM, A);
  part_f_MF = rep_matrix(0, P_MF, A);
              
  for (i in 1:P_MM){
    part_f_MM[i, :] = f_MM[map_indiv_to_age_MM[i]+1, :];
  }
  
  for (i in 1:P_FF){
    part_f_FF[i, :] = f_FF[map_indiv_to_age_FF[i]+1, :];
  }
  
  for (i in 1:P_FM){
    part_f_FM[i, :] = f_FM[map_indiv_to_age_FM[i]+1, :];
  }
  
  for (i in 1:P_MF){
    part_f_MF[i, :] = f_MF[map_indiv_to_age_MF[i]+1, :];
  }
  
  // print("part_f_MM =", part_f_MM)
  // print("part_f_FF =", part_f_FF)
  // print("part_f_FM =", part_f_FM)
  // print("part_f_MF =", part_f_MF)

  log_m_MM =  beta_0[MM] + part_f_MM;
  log_m_MF =  beta_0[MF] + part_f_MF;
  log_m_FM =  beta_0[FM] + part_f_FM;
  log_m_FF =  beta_0[FF] + part_f_FF;
  
  // will have to change this to a double for loop checking when Hib_c = 0, otherwise values of exp(log_m_MM) are too small!!
  alpha_strata_MM = (exp(log_m_MM).* H_MM) * map_age_to_strata / nu + epsilon;
  alpha_strata_MF = (exp(log_m_MF).* H_MF) * map_age_to_strata / nu + epsilon;
  alpha_strata_FM = (exp(log_m_FM).* H_FM) * map_age_to_strata / nu + epsilon;
  alpha_strata_FF = (exp(log_m_FF).* H_FF) * map_age_to_strata / nu + epsilon;
  
  // print("alpha_strata_MM =", alpha_strata_MM)
  // print("alpha_strata_FF =", alpha_strata_FF)
  // print("alpha_strata_FM =", alpha_strata_FM)
  // print("alpha_strata_MF =", alpha_strata_MF)
}

model
{
  // print("N_MM =", N_MM);
  // print("N_FF =", N_FF);
  // print("N_FM =", N_FM);
  // print("N_MF =", N_MF);
  // 
  // print("Y_MM =", Y_MM);
  // print("Y_FF =", Y_FF);
  // print("Y_FM =", Y_FM);
  // print("Y_MF =", Y_MF);
  // 
  // print("P_MM =", P_MM);
  // print("P_FF =", P_FF);
  // print("P_FM =", P_FM);
  // print("P_MF =", P_MF);
  // 
  // print("A =", A);
  // print("C =", C);
  // print("U =", U);
  // 
  // print("ROW_MAJOR_IDX_MM =", ROW_MAJOR_IDX_MM);
  // print("ROW_MAJOR_IDX_FF =", ROW_MAJOR_IDX_FF);
  // print("ROW_MAJOR_IDX_FM =", ROW_MAJOR_IDX_FM);
  // print("ROW_MAJOR_IDX_MF =", ROW_MAJOR_IDX_MF);
  // 
  // print("map_indiv_to_age_MM =", map_indiv_to_age_MM);
  // print("map_indiv_to_age_FF =", map_indiv_to_age_FF);
  // print("map_indiv_to_age_FM =", map_indiv_to_age_FM);
  // print("map_indiv_to_age_MF =", map_indiv_to_age_MF);
  // 
  // // print("map_indiv_to_age =", map_indiv_to_age);
  // 
  // print("H_MM =", H_MM);
  // print("H_FF =", H_FF);
  // print("H_FM =", H_FM);
  // print("H_MF =", H_MF);
  // 
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
  array[N] real log_lik;
  array[G,A,C] int yhat_strata;
  array[G] matrix[A,A] log_cnt_rate;
  
  log_cnt_rate[MM] = beta_0[MM] + f_MM;
  log_cnt_rate[FF] = beta_0[FF] + f_FF;
  log_cnt_rate[MF] = beta_0[MF] + f_MF;
  log_cnt_rate[FM] = beta_0[FM] + f_FM; 

  // {
  //   vector[N] alpha_strata_flat_indiv =
  //   append_row(
  //       append_row(
  //         append_row(
  //           to_vector(alpha_strata_MM')[ROW_MAJOR_IDX_MM],
  //           to_vector(alpha_strata_FF')[ROW_MAJOR_IDX_FF]
  //         ),
  //         to_vector(alpha_strata_MF')[ROW_MAJOR_IDX_MF]
  //       ),
  //     to_vector(alpha_strata_FM')[ROW_MAJOR_IDX_FM]
  //   );
  //   
  //   vector[U] alpha_strata_flat = alpha_strata_flat_indiv * map_indiv_to_age;
  //   vector[U] Y_age = Y * map_indiv_to_age;
  // 
  //   for(i in 1:U) {
  //     log_lik[i] = neg_binomial_lpmf( Y_age[i] | alpha_strata_flat[i], inv(nu) );
  //   }
  // }
}


