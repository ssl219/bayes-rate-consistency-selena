functions
{
#include gp-functions.stan
}
# 1 = household
# 2 = non-household
data
{
  int<lower=1> N_MM_1, N_FF_1, N_MF_1, N_FM_1, N_MM_2, N_FF_2, N_MF_2, N_FM_2; 
  // Number of observations for each gender and household combination

  int<lower=1> A;       // Number of age inputs
  int<lower=1> C;       // Number of age strata

  array[N_MM_1] int Y_MM_1; // Contacts for age i to ageband b
  array[N_FF_1] int Y_FF_1;
  array[N_MF_1] int Y_MF_1;
  array[N_FM_1] int Y_FM_1;

  array[N_MM_2] int Y_MM_2; 
  array[N_FF_2] int Y_FF_2;
  array[N_MF_2] int Y_MF_2;
  array[N_FM_2] int Y_FM_2;

  array[N_MM_1] int ROW_MAJOR_IDX_MM_1;
  array[N_FF_1] int ROW_MAJOR_IDX_FF_1;
  array[N_MF_1] int ROW_MAJOR_IDX_MF_1;
  array[N_FM_1] int ROW_MAJOR_IDX_FM_1;
  
  array[N_MM_2] int ROW_MAJOR_IDX_MM_2;
  array[N_FF_2] int ROW_MAJOR_IDX_FF_2;
  array[N_MF_2] int ROW_MAJOR_IDX_MF_2;
  array[N_FM_2] int ROW_MAJOR_IDX_FM_2;

  vector[A] log_N_M, log_N_F; // Participant size offsets
  vector[A] log_S_M_1, log_S_F_1, log_S_M_2, log_S_F_2; // Group contact offsets, replace log_S_M_1 with 1's? 
  row_vector[A] log_P_M, log_P_F; // Population size offsets

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
  int N = N_MM_1 + N_FF_1 + N_MF_1 + N_FM_1 + N_MM_2 + N_FF_2 + N_MF_2 + N_FM_2 // Total number of observations
  int MM_1 = 1, FF_1 = 2, MF_1 = 3, FM_1 = 4, MM_2 = 5, FF_2 = 6, MF_2 = 7, FM_2 = 8; // gender indexes
  int G = 8;                          // gender combinations
  real gp_delta = 1e-9;               // GP nugget
  real epsilon = 1e-13;               // Prevent shape parameter to be 0

  // Precompute offset terms
  array[G] matrix[A,A] log_offset;
  log_offset[MM_1] = rep_matrix(log_N_M + log_S_M_1, A) + rep_matrix(log_P_M, A);
  log_offset[FF_1] = rep_matrix(log_N_F + log_S_F_1, A) + rep_matrix(log_P_F, A);
  log_offset[MF_1] = rep_matrix(log_N_M + log_S_M_1, A) + rep_matrix(log_P_F, A);
  log_offset[FM_1] = rep_matrix(log_N_F + log_S_F_1, A) + rep_matrix(log_P_M, A);
  
  log_offset[MM_2] = rep_matrix(log_N_M + log_S_M_2, A) + rep_matrix(log_P_M, A);
  log_offset[FF_2] = rep_matrix(log_N_F + log_S_F_2, A) + rep_matrix(log_P_F, A);
  log_offset[MF_2] = rep_matrix(log_N_M + log_S_M_2, A) + rep_matrix(log_P_F, A);
  log_offset[FM_2] = rep_matrix(log_N_F + log_S_F_2, A) + rep_matrix(log_P_M, A);

  real L1, L2;
  matrix[A,M1] PHI1;
  matrix[A,M2] PHI2;

  // Precompute HSGP basis functions
  L1 = C1 * max(age_idx_std);
  L2 = C2 * max(age_idx_std);
  PHI1 = PHI(A, M1, L1, age_idx_std);
  PHI2 = PHI(A, M2, L2, age_idx_std);

  // append data
  array[N] int Y = append_array( 
    append_array( 
      append_array(
        append_array(
          append_array(
            append_array(
              append_array(
                append_array(Y_MM_1), Y_FF_1), Y_MF_1), Y_FM_1), Y_MM_2), Y_FF_2), Y_MF_2), Y_FM_2);
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

  matrix[A,A] f_MM_1, f_FF_1, f_MF_1, f_MM_2, f_FF_2, f_MF_2;
  array[G] matrix[A,A] log_cnt_rate;
  array[G] matrix<lower=0>[A,C] alpha_strata;

  f_MM_1 = hsgp(A, gp_alpha[MM_1], gp_rho_1[MM_1], gp_rho_2[MM_1],
              L1, L2, M1, M2, PHI1, PHI2, z[1:M1,]);
  f_FF_1 = hsgp(A, gp_alpha[FF_1], gp_rho_1[FF_1], gp_rho_2[FF_1],
              L1, L2, M1, M2, PHI1, PHI2, z[(M1+1):2*M1,]);
  f_MF_1 = hsgp(A, gp_alpha[MF_1], gp_rho_1[MF_1], gp_rho_2[MF_1],
              L1, L2, M1, M2, PHI1, PHI2, z[(2*M1+1):3*M1,]);
  
  f_MM_2 = hsgp(A, gp_alpha[MM_2], gp_rho_1[MM_2], gp_rho_2[MM_2],
              L1, L2, M1, M2, PHI1, PHI2, z[3*M1+1:4*M1,]);
  f_FF_2 = hsgp(A, gp_alpha[FF_2], gp_rho_1[FF_2], gp_rho_2[FF_2],
              L1, L2, M1, M2, PHI1, PHI2, z[(4*M1+1):5*M1,]);
  f_MF_2 = hsgp(A, gp_alpha[MF_2], gp_rho_1[MF_2], gp_rho_2[MF_2],
              L1, L2, M1, M2, PHI1, PHI2, z[(5*M1+1):6*M1,]);

  log_cnt_rate[MM_1] = beta_0[MM_1] + symmetrize_from_lower_tri(f_MM_1);
  log_cnt_rate[FF_1] = beta_0[FF_1] + symmetrize_from_lower_tri(f_FF_1);
  log_cnt_rate[MF_1] = beta_0[MF_1] + f_MF_1;
  log_cnt_rate[FM_1] = beta_0[FM_1] + f_MF_1'; # this is transposition?
  
  log_cnt_rate[MM_2] = beta_0[MM_2] + symmetrize_from_lower_tri(f_MM_2);
  log_cnt_rate[FF_2] = beta_0[FF_2] + symmetrize_from_lower_tri(f_FF_2);
  log_cnt_rate[MF_2] = beta_0[MF_2] + f_MF_2;
  log_cnt_rate[FM_2] = beta_0[FM_2] + f_MF_2'; 

  alpha_strata[MM_1] = exp(log_cnt_rate[MM_1] + log_offset[MM_1]) * map_age_to_strata / nu + epsilon;
  alpha_strata[FF_1] = exp(log_cnt_rate[FF_1] + log_offset[FF_1]) * map_age_to_strata / nu + epsilon;
  alpha_strata[MF_1] = exp(log_cnt_rate[MF_1] + log_offset[MF_1]) * map_age_to_strata / nu + epsilon;
  alpha_strata[FM_1] = exp(log_cnt_rate[FM_1] + log_offset[FM_1]) * map_age_to_strata / nu + epsilon;
  
  alpha_strata[MM_2] = exp(log_cnt_rate[MM_2] + log_offset[MM_2]) * map_age_to_strata / nu + epsilon;
  alpha_strata[FF_2] = exp(log_cnt_rate[FF_2] + log_offset[FF_2]) * map_age_to_strata / nu + epsilon;
  alpha_strata[MF_2] = exp(log_cnt_rate[MF_2] + log_offset[MF_2]) * map_age_to_strata / nu + epsilon;
  alpha_strata[FM_2] = exp(log_cnt_rate[FM_2] + log_offset[FM_2]) * map_age_to_strata / nu + epsilon;
}

model
{
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
    vector[N] alpha_strata_flat =
    append_row(
      append_row(
        append_row(
          append_row(
            append_row(
              append_row(
                append_row(
                  to_vector(alpha_strata[MM_1]')[ROW_MAJOR_IDX_MM_1],
                  to_vector(alpha_strata[FF_1]')[ROW_MAJOR_IDX_FF_1]
                ),
                to_vector(alpha_strata[MF_1]')[ROW_MAJOR_IDX_MF_1]
              ),
              to_vector(alpha_strata[FM_1]')[ROW_MAJOR_IDX_FM_1]
            ),
            to_vector(alpha_strata[MM_2]')[ROW_MAJOR_IDX_MM_2]
          ),
          to_vector(alpha_strata[FF_2]')[ROW_MAJOR_IDX_FF_2]
        ),
        to_vector(alpha_strata[MF_2]')[ROW_MAJOR_IDX_MF_2]
      ),
      to_vector(alpha_strata[FM_2]')[ROW_MAJOR_IDX_FM_2]
    );
    target += neg_binomial_lpmf( Y | alpha_strata_flat, inv(nu));
  }
}

generated quantities
{
  array[N] real log_lik;
  array[G,A,C] int yhat_strata;

  for(g in 1:G){
    for(i in 1:A){
      yhat_strata[g,i,:] = neg_binomial_rng( alpha_strata[g,i,:], inv(nu) );
    }
  }

  {
    vector[N] alpha_strata_flat =
    append_row(
      append_row(
        append_row(
          append_row(
            append_row(
              append_row(
                append_row(
                  to_vector(alpha_strata[MM_1]')[ROW_MAJOR_IDX_MM_1],
                  to_vector(alpha_strata[FF_1]')[ROW_MAJOR_IDX_FF_1]
                ),
                to_vector(alpha_strata[MF_1]')[ROW_MAJOR_IDX_MF_1]
              ),
              to_vector(alpha_strata[FM_1]')[ROW_MAJOR_IDX_FM_1]
            ),
            to_vector(alpha_strata[MM_2]')[ROW_MAJOR_IDX_MM_2]
          ),
          to_vector(alpha_strata[FF_2]')[ROW_MAJOR_IDX_FF_2]
        ),
        to_vector(alpha_strata[MF_2]')[ROW_MAJOR_IDX_MF_2]
      ),
      to_vector(alpha_strata[FM_2]')[ROW_MAJOR_IDX_FM_2]
    );

    for(i in 1:N) {
      log_lik[i] = neg_binomial_lpmf( Y[i] | alpha_strata_flat[i], inv(nu));
    }
  }
}

