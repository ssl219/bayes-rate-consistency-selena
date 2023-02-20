functions
{
#include gp-functions.stan
}

data
{
  int<lower=1> Ni_M, Ni_F; // Number of observations for a given participant i for each gender of the contact
  

  int<lower=1> A;       // Number of age inputs
  int<lower=1> C;       // Number of age strata

  array[Ni_M] int Yi_M; // Contacts of individual i to strata c
  array[Ni_F] int Yi_F;
  
  array[2] int vec_gender i // vector encoding the gender of participant i 

  array[Ni_M] int ROW_MAJOR_IDX_iM;
  array[Ni_F] int ROW_MAJOR_IDX_iF;
  
  vector[A] log_Hi_M, log_Hi_F; // Household members for participant i offsets for precise ages b
                                // in practice, would have to change this to only Hi_h because take the log of Hib_h = 0

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
  int Ni = Ni_M + Ni_F;  // Total number of observations for participant i
  int M = 1, F = 2; // gender indexes
  int MM = 1, FF = 2, MF = 3, FM = 4;
  int G = 2; // gender combinations
  int G_gp = 4
  real gp_delta = 1e-9;               // GP nugget
  real epsilon = 1e-13;               // Prevent shape parameter to be 0

  // Precompute offset terms
  array[G] matrix[1,A] log_offset_i;
  log_offset_i[M] = log_Hi_M;
  log_offset_i[F] = log_Hi_F;

  real L1, L2;
  matrix[A,M1] PHI1;
  matrix[A,M2] PHI2;

  // Precompute HSGP basis functions
  L1 = C1 * max(age_idx_std);
  L2 = C2 * max(age_idx_std);
  PHI1 = PHI(A, M1, L1, age_idx_std);
  PHI2 = PHI(A, M2, L2, age_idx_std);

  // append data
  array[Ni] int Yi = append_array(Yi_M, Yi_F)
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
  array[G] matrix[1,A] log_cnt_rate_i;
  array[G] matrix<lower=0>[1,C] alpha_strata_i;

  f_MM = hsgp(A, gp_alpha[MM], gp_rho_1[MM], gp_rho_2[MM],
              L1, L2, M1, M2, PHI1, PHI2, z[1:M1,]);
  f_FF = hsgp(A, gp_alpha[FF], gp_rho_1[FF], gp_rho_2[FF],
              L1, L2, M1, M2, PHI1, PHI2, z[(M1+1):2*M1,]);
  f_MF = hsgp(A, gp_alpha[MF], gp_rho_1[MF], gp_rho_2[MF],
              L1, L2, M1, M2, PHI1, PHI2, z[(2*M1+1):3*M1,]);
  f_FM = hsgp(A, gp_alpha[FM], gp_rho_1[FM], gp_rho_2[FM],
              L1, L2, M1, M2, PHI1, PHI2, z[(3*M1+1):4*M1,]);

  log_cnt_rate_i[M] =  beta_0[M] + (vec_gender_i[M]*f_MM + vec_gender_i[F]*f_FM)[map_indiv_to_age, :]
  log_cnt_rate_i[F] =  beta_0[F] + (vec_gender_i[M]*f_MF + vec_gender_i[F]*f_FF)[map_indiv_to_age, :]

  alpha_strata_i[M] = exp(log_cnt_rate_i[M] + log_offset_i[M]) * map_age_to_strata / nu + epsilon;
  alpha_strata_i[F] = exp(log_cnt_rate_i[F] + log_offset_i[F]) * map_age_to_strata / nu + epsilon;
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
    vector[Ni] alpha_strata_i_flat =
          append_row(
            to_vector(alpha_strata_i[M]')[ROW_MAJOR_IDX_iM],
            to_vector(alpha_strata_i[F]')[ROW_MAJOR_IDX_iF]
          )
          ;
    target += neg_binomial_lpmf( Y | alpha_strata_i_flat, inv(nu));
  }
}

generated quantities
{
  array[Ni] real log_lik;
  array[G,1,C] int yhat_strata;

  for(g in 1:G){
      yhat_strata[g,1,:] = neg_binomial_rng( alpha_strata_i[g, 1,:], inv(nu) );
    }

  {
     vector[Ni] alpha_strata_i_flat =
          append_row(
            to_vector(alpha_strata_i[M]')[ROW_MAJOR_IDX_iM],
            to_vector(alpha_strata_i[F]')[ROW_MAJOR_IDX_iF]
          )
          ;

    for(i in 1:Ni) {
      log_lik[i] = neg_binomial_lpmf( Y[i] | alpha_strata_i_flat[i], inv(nu));
    }
  }
}

