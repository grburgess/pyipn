functions {
  vector log_f_predict(int N, real alpha, real rho, real[] x_predict, vector f_tilde) {

    matrix[N, N] cov;
    matrix[N, N] L_cov;
    cov =   cov_exp_quad(x_predict, alpha, rho) + diag_matrix(rep_vector(1e-10, N));
    L_cov = cholesky_decompose(cov);    
    return  L_cov * f_tilde;


  }

}



data {
  int<lower=1> N;
  real time[N];
  int counts[N];

  int<lower=1> N2;
  real time2[N2];
  int counts2[N2];
  

}

transformed data {
}

parameters {

  real<lower=0> rho;
  real<lower=0> alpha;

  real<lower=0> m;


  real<lower=0> m2;
  real<lower=0> dT;
  
  
  vector[N] f_tilde;
  vector[N2] f_tilde2;
}

transformed parameters {  
  real x_shift[N2];
  for(n in 1:N2) {
    x_shift[n] = time2[n] - dT;
  }
  
}

model {
  
  rho ~ inv_gamma(6.57075, 41.8423);
  alpha ~ normal(0, 2);

  m ~ normal(40,10);  
  m2 ~ normal(40,10);

  dT ~ normal(10,50);
  
  
  f_tilde ~ std_normal();
  f_tilde2 ~ std_normal();


  counts ~ poisson(exp(log_f_predict(N, alpha, rho, time, f_tilde))  + m);

  counts2 ~ poisson(exp(log_f_predict(N2, alpha, rho, x_shift, f_tilde2)) + m2);

  
}

generated quantities {
  vector[N] f_predict = exp(log_f_predict(N, alpha, rho, time, f_tilde));
  //vector[N] y_predict;

  vector[N2] f_predict2 = exp(log_f_predict(N2, alpha, rho, time2, f_tilde2));
  //vector[N2] y_predict2;

  
/*   for (n in 1:N) */
/*     y_predict[n] = poisson_log_rng(log_f_predict[n]); */


  
/*   for (n in 1:N2) */
/*     y_predict2[n] = poisson_log_rng(log_f_predict2[n]); */


  
}
