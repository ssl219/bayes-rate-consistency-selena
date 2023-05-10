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
  * @param A: Number of ages for particiapnts
  * @param alpha: GP scaling parameter
  * @param rho1, rho2: GP length-scale parameter on participants and contacts age dimensions
  * @param L1, L2: HSGP parameters
  * @param M1, M2: HSGP parameters
  * @param z:
  * @return A two dimensional Gaussian process fuction
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
  
  // matrix[N_MM + N_FF + N_MF + N_FM, U] map_indiv_to_age; // map individual-age space to age-age space


  array[A_MM] real log_H_MM; // Log of household offsets (ordered in ascending order of participant age then contact age)
  array[A_FF] real log_H_FF;
  array[A_MF] real log_H_MF;
  array[A_FM] real log_H_FM;
 
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
  int MM = 1, FF = 2, MF = 3, FM = 4;
  int G = 4; // gender combinations
  real gp_delta = 1e-9;               // GP nugget
  real epsilon = 1e-13;               // Prevent shape parameter to be 0
  
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
  vector [G] beta_0;

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


 f_MM = hsgp_restruct(A, gp_alpha[MM], gp_rho_1[MM], gp_rho_2[MM],
                       L1, L2, M1, M2, PHI1, PHI2,
                       z[1:M1,], NN_IDX);
  f_FF = hsgp_restruct(A, gp_alpha[FF], gp_rho_1[FF], gp_rho_2[FF],
                       L1, L2, M1, M2, PHI1, PHI2,
                       z[(M1+1):2*M1,], NN_IDX);
  f_MF = hsgp_restruct(A, gp_alpha[MF], gp_rho_1[MF], gp_rho_2[MF],
                       L1, L2, M1, M2, PHI1, PHI2,
                       z[(2*M1+1):3*M1,], NN_IDX);

  alpha_MM = rep_matrix(epsilon, P_MM, A); // initialize alpha_strata to zeros
  alpha_FF = rep_matrix(epsilon, P_FF, A);
  alpha_MF = rep_matrix(epsilon, P_MF, A);
  alpha_FM = rep_matrix(epsilon, P_FM, A);


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
  
  alpha_strata_MM = alpha_MM * map_age_to_strata + epsilon;
  alpha_strata_MF = alpha_MF * map_age_to_strata + epsilon;
  alpha_strata_FM = alpha_FM * map_age_to_strata + epsilon;
  alpha_strata_FF = alpha_FF * map_age_to_strata + epsilon;

}

model
{
  // GP priors
  target += inv_gamma_lpdf(gp_rho_1 | 5, 5);
  target += inv_gamma_lpdf(gp_rho_2 | 5, 5);
  target += cauchy_lpdf(gp_alpha | 0, 1);
  target += std_normal_lpdf( to_vector(z) );

  // Overdispersion
  // target += exponential_lpdf(nu | 1);

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
    target += poisson_lpmf( Y | alpha_strata_flat_indiv);
  }
}

generated quantities
{
  array[G] matrix[A,A] log_cnt_rate;
  matrix[P_MM, A] alpha_age_MM;
  matrix[P_FF, A] alpha_age_FF;
  matrix[P_FM, A] alpha_age_FM;
  matrix[P_MF, A] alpha_age_MF;
  array[P_MM] int part_age_MM; 
  array[P_FF] int part_age_FF; 
  array[P_FM] int part_age_FM; 
  array[P_MF] int part_age_MF; 
  array[A_MM] int contact_age_MM; 
  array[A_FF] int contact_age_FF;
  array[A_MF] int contact_age_MF;
  array[A_FM] int contact_age_FM;

  array[P_MM + 1] int cumulative_MM; 
  array[P_FF + 1] int cumulative_FF;
  array[P_FM + 1] int cumulative_FM;
  array[P_MF + 1] int cumulative_MF;
  array[N] real log_lik;
  
  /*array[N_MM] int yhat_strata_MM;
  array[N_FF] int yhat_strata_FF;
  array[N_MF] int yhat_strata_MF;
  array[N_FM] int yhat_strata_FM;*/
  
  array[N_MM] int row_major_idx_MM;
  array[N_FF] int row_major_idx_FF;
  array[N_MF] int row_major_idx_MF;
  array[N_FM] int row_major_idx_FM;
  
 
  {
  log_cnt_rate[MM] = beta_0[MM] + symmetrize_from_lower_tri(f_MM);
  log_cnt_rate[FF] = beta_0[FF] + symmetrize_from_lower_tri(f_FF);
  log_cnt_rate[MF] = beta_0[MF] + f_MF;
  log_cnt_rate[FM] = beta_0[FM] + f_MF'; 
  
  alpha_age_MM = alpha_MM;
  alpha_age_MF = alpha_MF;
  alpha_age_FM = alpha_FM;
  alpha_age_FF = alpha_FF;

  part_age_MM = map_indiv_to_age_MM;
  part_age_FF = map_indiv_to_age_FF;
  part_age_FM = map_indiv_to_age_FM;
  part_age_MF = map_indiv_to_age_MF;
  
  row_major_idx_MM = ROW_MAJOR_IDX_MM;
  row_major_idx_FF = ROW_MAJOR_IDX_FF;
  row_major_idx_MF = ROW_MAJOR_IDX_MF;
  row_major_idx_FM = ROW_MAJOR_IDX_FM;
  
  /*// yhat_strata, not extracting values for all strata ages to improve computational efficiency
  vector[N_MM] alpha_strata_flat_MM = to_vector(alpha_strata_MM')[ROW_MAJOR_IDX_MM];
  vector[N_FF] alpha_strata_flat_FF = to_vector(alpha_strata_FF')[ROW_MAJOR_IDX_FF];
  vector[N_FM] alpha_strata_flat_FM = to_vector(alpha_strata_FM')[ROW_MAJOR_IDX_FM];
  vector[N_MF] alpha_strata_flat_MF = to_vector(alpha_strata_MF')[ROW_MAJOR_IDX_MF];
  
  yhat_strata_MM = poisson_rng(alpha_strata_flat_MM);
  yhat_strata_FF = poisson_rng(alpha_strata_flat_FF);
  yhat_strata_MF = poisson_rng(alpha_strata_flat_MF);
  yhat_strata_FM = poisson_rng(alpha_strata_flat_FM);
  }*/
  
  //try log-lik for individual-strata space, not age-age space as done below
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
    for(i in 1:N) {
    log_lik[i] = poisson_lpmf( Y[i] | alpha_strata_flat_indiv[i]);
    }
  }
  
}

  
  

// In the words of Copilot: I have a model that I am trying to fit to some data. The model is a hierarchical model with a GP prior on the random effects. I am trying to fit the model using the NUTS sampler. The model is a bit complicated, but I have tried to simplify it as much as possible. 
  


 