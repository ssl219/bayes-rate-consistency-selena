functions
{
vector diagSPD_EQ(real alpha, real rho, real L, int M)
  {
    return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
  }

vector diagSPD_Matern32(real alpha, real rho, real L, int M)
{
  return 2*alpha * (sqrt(3)/rho)^1.5 * inv((sqrt(3)/rho)^2 + ((pi()/2/L) * linspaced_vector(M, 1, M))^2);
}

vector diagSPD_Matern52(real alpha, real rho, real L, int M)
{
  return 2*alpha * sqrt(4.0/3) * (sqrt(5)/rho)^2.5 * inv((sqrt(5)/rho)^2 + ((pi()/2/L) * linspaced_vector(M, 1, M))^2)^1.5;
}

matrix PHI(int N, int M, real L, vector x)
{
  return sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), linspaced_vector(M, 1, M)))/sqrt(L);
}

/** Kronecker multivariate product
  *
  * Enables efficient sampling from a 2D multivariate normal distribution.
  *
  * @param A
  * @param B
  * @param V A matrix of N(0,1) distributed random variables
  */
matrix kron_mvprod(matrix A, matrix B, matrix V)
{
  return (B*V) * transpose(A);
}

// Transform restructured matrix back to a squared matrix
matrix inverse_restruct(matrix B, array[] int nn_idx)
{
  return to_matrix(to_vector(B)[nn_idx], cols(B), cols(B), 0);
}

/** Kronecker Decomposed 2D Gaussian process
  *
  * @param x1, x2: participant and contact age
  * @param delta: GP nugget
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on the participants and contacts dimensions
  * @param z:
  * @return A two dimensional Gaussian process function
  */
matrix gp2d(array[] real x1, array[] real x2,
            real alpha, real rho1, real rho2, matrix z)
{
  int A = size(x1);
  int B = size(x2);

  // Compute the exponentiated quadratic covariance function:
  matrix[A,A] K1 = gp_exp_quad_cov(x1, alpha, rho1) + diag_matrix(rep_vector(1e-9, A));
  matrix[B,B] K2 = gp_exp_quad_cov(x2, alpha, rho2) + diag_matrix(rep_vector(1e-9, B));

  // Cholesky Decomposition K1 = L_K1 L_K1^T
  matrix[A,A] L_K1 = cholesky_decompose(K1);
  matrix[B,B] L_K2 = cholesky_decompose(K2);

  matrix[A,B] f = kron_mvprod(L_K2, L_K1, z);

  return(f);
}

/** Kronecker Decomposed 2D Gaussian process (Matern 3/2)
  *
  * @param x1, x2: participant and contact age
  * @param delta: GP nugget
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on the participants and contacts dimensions
  * @param z:
  * @return A two dimensional Gaussian process function
  */
matrix gp2d_matern32(array[] real x1, array[] real x2,
                     real alpha, real rho1, real rho2, matrix z)
{
  int A = size(x1);
  int B = size(x2);

  // Compute the exponentiated quadratic covariance function:
  matrix[A,A] K1 = gp_matern32_cov(x1, alpha, rho1) + diag_matrix(rep_vector(1e-9, A));
  matrix[B,B] K2 = gp_matern32_cov(x2, alpha, rho2) + diag_matrix(rep_vector(1e-9, B));

  // Cholesky Decomposition K1 = L_K1 L_K1^T
  matrix[A,A] L_K1 = cholesky_decompose(K1);
  matrix[B,B] L_K2 = cholesky_decompose(K2);

  matrix[A,B] f = kron_mvprod(L_K2, L_K1, z);

  return(f);
}

/** Kronecker Decomposed 2D Gaussian process (Matern 5/2)
  *
  * @param x1, x2: participant and contact age
  * @param delta: GP nugget
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on the participants and contacts dimensions
  * @param z:
  * @return A two dimensional Gaussian process function
  */
matrix gp2d_matern52(array[] real x1, array[] real x2,
                     real alpha, real rho1, real rho2, matrix z)
{
  int A = size(x1);
  int B = size(x2);

  // Compute the exponentiated quadratic covariance function:
  matrix[A,A] K1 = gp_matern52_cov(x1, alpha, rho1) + diag_matrix(rep_vector(1e-9, A));
  matrix[B,B] K2 = gp_matern52_cov(x2, alpha, rho2) + diag_matrix(rep_vector(1e-9, B));

  // Cholesky Decomposition K1 = L_K1 L_K1^T
  matrix[A,A] L_K1 = cholesky_decompose(K1);
  matrix[B,B] L_K2 = cholesky_decompose(K2);

  matrix[A,B] f = kron_mvprod(L_K2, L_K1, z);

  return(f);
}

/** Restructured Kronecker Decomposed 2D Gaussian process
  *
  * @param x1, x2: Age difference and cohort age
  * @param delta GP nugget
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on age difference and cohort age dimensions
  * @param z:
  * @param nn_idx: Index of non-nuiance parameters
  * @return A two dimensional Gaussian process fuction
  */
matrix gp2d_restruct(array[] real x1, array[] real x2,
                     real alpha, real rho1, real rho2, matrix z,
                     array[] int nn_idx)
{
  int A = size(x1);
  int B = size(x2);
  matrix[A,B] f_restruct = gp2d(x1, x2, alpha, rho1, rho2, z);

  return(inverse_restruct(f_restruct, nn_idx));
}

/** Restructured Kronecker Decomposed 2D Gaussian process (Matern 3/2)
  *
  * @param x1, x2: Age difference and cohort age
  * @param delta GP nugget
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on age difference and cohort age dimensions
  * @param z:
  * @param nn_idx: Index of non-nuiance parameters
  * @return A two dimensional Gaussian process fuction
  */
matrix gp2d_matern32_restruct(array[] real x1, array[] real x2,
                              real alpha, real rho1, real rho2, matrix z,
                              array[] int nn_idx)
{
  int A = size(x1);
  int B = size(x2);
  matrix[A,B] f_restruct = gp2d_matern32(x1, x2, alpha, rho1, rho2, z);

  return(inverse_restruct(f_restruct, nn_idx));
}

/** Restructured Kronecker Decomposed 2D Gaussian process (Matern 5/2)
  *
  * @param x1, x2: Age difference and cohort age
  * @param delta GP nugget
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on age difference and cohort age dimensions
  * @param z:
  * @param nn_idx: Index of non-nuiance parameters
  * @return A two dimensional Gaussian process fuction
  */
matrix gp2d_matern52_restruct(array[] real x1, array[] real x2,
                              real alpha, real rho1, real rho2, matrix z,
                              array[] int nn_idx)
{
  int A = size(x1);
  int B = size(x2);
  matrix[A,B] f_restruct = gp2d_matern52(x1, x2, alpha, rho1, rho2, z);

  return(inverse_restruct(f_restruct, nn_idx));
}

/** Hilbert Space approximate 2D Gaussian process
  *
  * @param A: Number of ages for participants
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on participants and contacts age dimensions
  * @param L1, L2: HSGP parameters
  * @param M1, M2: HSGP parameters
  * @param z:
  * @return A two dimensional Gaussian process function
  */
matrix hsgp(int A, real alpha, real rho1, real rho2, real L1, real L2, int M1, int M2,
            matrix PHI1, matrix PHI2, matrix z)
{
  vector[M1] sqrt_spd_1 = diagSPD_EQ(alpha, rho1, L1, M1);
  vector[M2] sqrt_spd_2 = diagSPD_EQ(alpha, rho2, L2, M2);

  matrix[A,A] f = kron_mvprod(
    diag_post_multiply( PHI1, sqrt_spd_1 ),
    diag_post_multiply( PHI2, sqrt_spd_2 ),
    z
  );

  return(f);
}

/** Hilbert Space approximate 2D Gaussian process (Matern 3/2)
  *
  * @param A: Number of ages for particiapnts
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on participants and contacts age dimensions
  * @param L1, L2: HSGP parameters
  * @param M1, M2: HSGP parameters
  * @param z:
  * @return A two dimensional Gaussian process fuction
  */
matrix hsgp_matern32(int A, real alpha, real rho1, real rho2, real L1, real L2, int M1, int M2,
            matrix PHI1, matrix PHI2, matrix z)
{
  vector[M1] sqrt_spd_1 = diagSPD_Matern32(alpha, rho1, L1, M1);
  vector[M2] sqrt_spd_2 = diagSPD_Matern32(alpha, rho2, L2, M2);

  matrix[A,A] f = kron_mvprod(
    diag_post_multiply( PHI1, sqrt_spd_1 ),
    diag_post_multiply( PHI2, sqrt_spd_2 ),
    z
  );

  return(f);
}

/** Hilbert Space approximate 2D Gaussian process (Matern 5/2)
  *
  * @param A: Number of ages for particiapnts
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on participants and contacts age dimensions
  * @param L1, L2: HSGP parameters
  * @param M1, M2: HSGP parameters
  * @param z:
  * @return A two dimensional Gaussian process fuction
  */
matrix hsgp_matern52(int A, real alpha, real rho1, real rho2, real L1, real L2, int M1, int M2,
            matrix PHI1, matrix PHI2, matrix z)
{
  vector[M1] sqrt_spd_1 = diagSPD_Matern52(alpha, rho1, L1, M1);
  vector[M2] sqrt_spd_2 = diagSPD_Matern52(alpha, rho2, L2, M2);

  matrix[A,A] f = kron_mvprod(
    diag_post_multiply( PHI1, sqrt_spd_1 ),
    diag_post_multiply( PHI2, sqrt_spd_2 ),
    z
  );

  return(f);
}

/** Restructured Hilbert Space approximate 2D Gaussian process
  *
  * @param A: Number of ages for particiapnts
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on age difference and cohort age dimensions
  * @param L1, L2: HSGP parameters
  * @param M1, M2: HSGP parameters
  * @param z:
  * @param nn_idx: Index of non-nuiance parameters
  * @return A two dimensional Gaussian process fuction
  */
matrix hsgp_restruct(int A, real alpha, real rho1, real rho2, real L1, real L2, int M1, int M2,
                     matrix PHI1, matrix PHI2, matrix z, array[] int nn_idx)
{
  vector[M1] sqrt_spd_1 = diagSPD_EQ(alpha, rho1, L1, M1);
  vector[M2] sqrt_spd_2 = diagSPD_EQ(alpha, rho2, L2, M2);

  matrix[2*A-1, A] f_restruct = kron_mvprod(
    diag_post_multiply( PHI2, sqrt_spd_2 ),
    diag_post_multiply( PHI1, sqrt_spd_1 ),
    z
  );

  return inverse_restruct(f_restruct, nn_idx);
}

/** Restructured Hilbert Space approximate 2D Gaussian process (Matern 3/2)
  *
  * @param A: Number of ages for particiapnts
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on age difference and cohort age dimensions
  * @param L1, L2: HSGP parameters
  * @param M1, M2: HSGP parameters
  * @param z:
  * @param nn_idx: Index of non-nuiance parameters
  * @return A two dimensional Gaussian process fuction
  */
matrix hsgp_matern32_restruct(int A, real alpha, real rho1, real rho2, real L1, real L2, int M1, int M2,
                     matrix PHI1, matrix PHI2, matrix z, array[] int nn_idx)
{
  vector[M1] sqrt_spd_1 = diagSPD_Matern32(alpha, rho1, L1, M1);
  vector[M2] sqrt_spd_2 = diagSPD_Matern32(alpha, rho2, L2, M2);

  matrix[2*A-1, A] f_restruct = kron_mvprod(
    diag_post_multiply( PHI2, sqrt_spd_2 ),
    diag_post_multiply( PHI1, sqrt_spd_1 ),
    z
  );

  return inverse_restruct(f_restruct, nn_idx);
}

/** Restructured Hilbert Space approximate 2D Gaussian process (Matern 5/2)
  *
  * @param A: Number of ages for particiapnts
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on age difference and cohort age dimensions
  * @param L1, L2: HSGP parameters
  * @param M1, M2: HSGP parameters
  * @param z:
  * @param nn_idx: Index of non-nuiance parameters
  * @return A two dimensional Gaussian process fuction
  */
matrix hsgp_matern52_restruct(int A, real alpha, real rho1, real rho2, real L1, real L2, int M1, int M2,
                              matrix PHI1, matrix PHI2, matrix z, array[] int nn_idx)
{
  vector[M1] sqrt_spd_1 = diagSPD_Matern52(alpha, rho1, L1, M1);
  vector[M2] sqrt_spd_2 = diagSPD_Matern52(alpha, rho2, L2, M2);

  matrix[2*A-1, A] f_restruct = kron_mvprod(
    diag_post_multiply( PHI2, sqrt_spd_2 ),
    diag_post_multiply( PHI1, sqrt_spd_1 ),
    z
  );

  return inverse_restruct(f_restruct, nn_idx);
}

}

data
{
  int<lower=1> N_MM, N_FF, N_MF, N_FM; // Number of observations for each gender pair

  int<lower=1> A;       // Number of age inputs
  int<lower=1> C;       // Number of age strata

  array[N_MM] int Y_MM; // Contacts for age i to ageband b
  array[N_FF] int Y_FF;
  array[N_MF] int Y_MF;
  array[N_FM] int Y_FM;

  array[N_MM] int ROW_MAJOR_IDX_MM;
  array[N_FF] int ROW_MAJOR_IDX_FF;
  array[N_MF] int ROW_MAJOR_IDX_MF;
  array[N_FM] int ROW_MAJOR_IDX_FM;

  vector[A] log_N_M, log_N_F; // Participant size offsets
  vector[A] log_S_M, log_S_F; // Group contact offsets
  row_vector[A] log_P_M, log_P_F; // Population size offsets

  vector[A] age_idx_std;          // Standardized age index
  vector[2*A-1] diff_idx_std;     // Standardized age difference index
  matrix[A,C] map_age_to_strata;  // Indicator Matrix that maps age to age strata
  array[A*A] int NN_IDX;          // Index indicating the locations of the non-nuisance parameters in the resturctured HSGP matrix

  // HSGP parameters
  real<lower=0> C1; // Factor to determine the boundary value L (cohort age dimension)
  int<lower=1> M1;  // Number of basis functions (cohort age dimension)
  real<lower=0> C2; // Factor to determine the boundary value L for age of contacted individuals (age difference dimension)
  int<lower=1> M2;  // Number of basis functions (age difference dimension)
}

transformed data
{
  int N = N_MM + N_FF + N_MF + N_FM;  // Total number of observations
  int MM = 1, FF = 2, MF = 3, FM = 4; // gender indexes
  int G = 4;                          // gender combinations
  real gp_delta = 1e-9;               // GP nugget
  real epsilon = 1e-13;               // Prevent shape parameter to be 0

  // Precompute offset terms
  array[G] matrix[A,A] log_offset;
  log_offset[MM] = rep_matrix(log_N_M + log_S_M, A) + rep_matrix(log_P_M, A);
  log_offset[FF] = rep_matrix(log_N_F + log_S_F, A) + rep_matrix(log_P_F, A);
  log_offset[MF] = rep_matrix(log_N_M + log_S_M, A) + rep_matrix(log_P_F, A);
  log_offset[FM] = rep_matrix(log_N_F + log_S_F, A) + rep_matrix(log_P_M, A);

  real L1, L2;
  matrix[2*A-1, M1] PHI1;
  matrix[A, M2] PHI2;

  // Precompute HSGP basis functions
  L1 = C1 * max(diff_idx_std);
  L2 = C2 * max(age_idx_std);
  PHI1 = PHI(2*A-1, M1, L1, diff_idx_std);
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
  vector<lower=0, upper=pi()/2 >[G-1] gp_alpha_unif;

  matrix[(G-1)*M1, M2] z; // HSGP basis function coefficients
}

transformed parameters
{
  vector<lower=0>[G-1] gp_alpha = tan(gp_alpha_unif); // Reparametrize Half-Cauchy for optimization in HMC

  matrix[A, A] f_MM, f_FF, f_MF;
  array[G] matrix[A, A] log_cnt_rate;
  array[G] matrix<lower=0>[A, C] alpha_strata;

  f_MM = hsgp_restruct(A, gp_alpha[MM], gp_rho_1[MM], gp_rho_2[MM],
                       L1, L2, M1, M2, PHI1, PHI2,
                       z[1:M1,], NN_IDX);
  f_FF = hsgp_restruct(A, gp_alpha[FF], gp_rho_1[FF], gp_rho_2[FF],
                       L1, L2, M1, M2, PHI1, PHI2,
                       z[(M1+1):2*M1,], NN_IDX);
  f_MF = hsgp_restruct(A, gp_alpha[MF], gp_rho_1[MF], gp_rho_2[MF],
                       L1, L2, M1, M2, PHI1, PHI2,
                       z[(2*M1+1):3*M1,], NN_IDX);

  log_cnt_rate[MM] = beta_0[MM] + symmetrize_from_lower_tri(f_MM);
  log_cnt_rate[FF] = beta_0[FF] + symmetrize_from_lower_tri(f_FF);
  log_cnt_rate[MF] = beta_0[MF] + f_MF;
  log_cnt_rate[FM] = beta_0[FM] + f_MF';

  alpha_strata[MM] = exp(log_cnt_rate[MM] + log_offset[MM]) * map_age_to_strata / nu + epsilon;
  alpha_strata[FF] = exp(log_cnt_rate[FF] + log_offset[FF]) * map_age_to_strata / nu + epsilon;
  alpha_strata[MF] = exp(log_cnt_rate[MF] + log_offset[MF]) * map_age_to_strata / nu + epsilon;
  alpha_strata[FM] = exp(log_cnt_rate[FM] + log_offset[FM]) * map_age_to_strata / nu + epsilon;
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

