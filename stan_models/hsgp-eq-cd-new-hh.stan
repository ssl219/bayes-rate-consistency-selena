functions
{
#include gp-functions.stan
}

data
{
  int<lower=1> N_M, N_F; // Number of observations for all participants for each gender of the contact
  
  int<lower=1> P;       // Number of participants
  int<lower=1> A;       // Number of age inputs
  int<lower=1> C;       // Number of age strata

  array[N_M] int Y_M; // Contacts for all participants (ordered) to strata c
  array[N_F] int Y_F;
  
  array[P] int map_indiv_to_age; // array of ages of each participant
  vector[P] int vec_gender_M, vec_gender_F; // vectors encoding gender of each participant

  matrix[P, A] H_M, H_F; // Household size offsets

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
  int N = N_M + N_F;  // Total number of observations
  int M = 1, F = 2; // gender indexes
  int MM = 1, FF = 2, MF = 3, FM = 4;
  int G = 2; // gender combinations
  int G_gp = 4
  real gp_delta = 1e-9;               // GP nugget
  real epsilon = 1e-13;               // Prevent shape parameter to be 0

  // Precompute offset terms
  array[G] matrix[P,A] offset;
  offset[M] = H_M;
  offset[F] = H_F;

  real L1, L2;
  matrix[A,M1] PHI1;
  matrix[A,M2] PHI2;

  // Precompute HSGP basis functions
  L1 = C1 * max(age_idx_std);
  L2 = C2 * max(age_idx_std);
  PHI1 = PHI(A, M1, L1, age_idx_std);
  PHI2 = PHI(A, M2, L2, age_idx_std);

  // append data
  array[N] int Y = append_array(Y_M, Y_F);
}

parameters
{
  vector[G] beta_0;
  real<lower=0> nu;

  vector<lower=0>[G_gp] gp_rho_1;
  vector<lower=0>[G_gp] gp_rho_2;
  vector<lower=0, upper=pi()/2>[G_gp] gp_alpha_unif;

  matrix[(G_gp)*M1, M2] z; // HSGP basis function coefficients
}

transformed parameters
{
  vector<lower=0>[G_gp] gp_alpha = tan(gp_alpha_unif); // Reparametrize Half-Cauchy for optimization in HMC

  matrix[A,A] f_MM, f_FF, f_MF, f_FM;
  array[G] matrix[P,C] log_m_h;
  array[G] matrix<lower=0>[P,C] alpha_strata_h;
  row_vector[A] vec_f_MM, vec_f_FF, vec_f_FM, vec_f_MF;
  

  f_MM = hsgp(A, gp_alpha[MM], gp_rho_1[MM], gp_rho_2[MM],
              L1, L2, M1, M2, PHI1, PHI2, z[1:M1,]);
  f_FF = hsgp(A, gp_alpha[FF], gp_rho_1[FF], gp_rho_2[FF],
              L1, L2, M1, M2, PHI1, PHI2, z[(M1+1):2*M1,]);
  f_MF = hsgp(A, gp_alpha[MF], gp_rho_1[MF], gp_rho_2[MF],
              L1, L2, M1, M2, PHI1, PHI2, z[(2*M1+1):3*M1,]);
              
  vec_f_MM[1, :] = f_MM[map_indiv_to_age, :]
  vec_f_FF[1, :] = f_FF[map_indiv_to_age, :]
  vec_f_FM[1, :] = f_FM[map_indiv_to_age, :]
  vec_f_MF[1, :] = f_MF[map_indiv_to_age, :]

  log_m_h[M] =  beta_0[M] + vec_gender_M * vec_f_MM + vec_gender_F * vec_f_FM
  log_m_h[F] =  beta_0[F] + vec_gender_F * vec_f_MF + vec_gender_F * vec_f_FF
  
  alpha_strata_h[M] = (exp(log_m_h[M]).* offset[M]) * map_age_to_strata / nu + epsilon;
  alpha_strata_h[F] = (exp(log_m_h[F]).* offset[F]) * map_age_to_strata / nu + epsilon;
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
            to_vector(alpha_strata[MM]')[ROW_MAJOR_IDX_MM],
            to_vector(alpha_strata[FF]')[ROW_MAJOR_IDX_FF]
          ),
          to_vector(alpha_strata[MF]')[ROW_MAJOR_IDX_MF]
        ),
      to_vector(alpha_strata[FM]')[ROW_MAJOR_IDX_FM]
    );
    target += neg_binomial_lpmf( Y | alpha_strata_flat, inv(nu));
  }
}

generated quantities
{
  array[N] real log_lik;
  array[G,A,C] int yhat_strata;
  array[G_gp] matrix[A,A] log_cnt_rate;
  
  log_cnt_rate[MM] = beta_0[M] + symmetrize_from_lower_tri(f_MM);
  log_cnt_rate[FF] = beta_0[F] + symmetrize_from_lower_tri(f_FF);
  log_cnt_rate[MF] = beta_0[F] + f_MF;
  log_cnt_rate[FM] = beta_0[M] + f_MF'; 
  

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
            to_vector(alpha_strata[MM]')[ROW_MAJOR_IDX_MM],
            to_vector(alpha_strata[FF]')[ROW_MAJOR_IDX_FF]
          ),
          to_vector(alpha_strata[MF]')[ROW_MAJOR_IDX_MF]
        ),
      to_vector(alpha_strata[FM]')[ROW_MAJOR_IDX_FM]
    );

    for(i in 1:N) {
      log_lik[i] = neg_binomial_lpmf( Y[i] | alpha_strata_flat[i], inv(nu));
    }
  }
}

