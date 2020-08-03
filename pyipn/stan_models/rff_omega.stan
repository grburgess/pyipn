functions {
#include functions.stan
}

data {


  int N_detectors; // number of detectors used
  int N_time_bins[N_detectors]; // the number of time bins in each detector
  int max_N_time_bins; // the max number of time bins

  int counts[N_detectors, max_N_time_bins]; // matrix of counts per detector

  vector[max_N_time_bins] time[N_detectors]; // matrix of time bin midpoints per detector
  vector[max_N_time_bins] exposure[N_detectors]; // matrix of exposures per detector


  vector[3] sc_pos[N_detectors]; // the 3D vector for each space craft

  int<lower=1> k; // number of FFs

  int grainsize[N_detectors];


}
transformed data {


}
parameters {

  vector[k] beta1; // the amplitude along the cos basis
  vector[k] beta2; // the amplitude along the sin basis

  row_vector[k] omega_var[2]; // this weird MC integration thing.

  real  log_bkg;
  
  vector[2] log_scale;
  
  //positive_ordered[2] bw;
  ordered[2] log_bw;


}

transformed parameters {

  real bkg = exp(log_bkg);



  vector[2] scale = exp(log_scale) * inv_sqrt(k);
  vector[2] bw = exp(log_bw);
  

  row_vector[k] omega[2]; // this weird MC integration thing.
  
  // non-center
  omega[1] = omega_var[1] * bw[1];
  omega[2] = omega_var[2] * bw[2];
  
  
}

model {

  // priors

  beta1 ~ std_normal();
  beta2 ~ std_normal();

  log_scale ~ std_normal();

  log_bkg ~ normal(log(500), log(100));
  log_bw ~ normal(0, .5);

  //bw ~ std_normal();
  
  omega_var[1] ~ std_normal();
  omega_var[2] ~ std_normal();

  target += reduce_sum(partial_log_like_bw_multi_scale, counts[1], grainsize[1],
                       time[1], exposure[1],
                       omega[1], omega[2], beta1, beta2,
                       0., bkg, scale[1], scale[2], 1., k);

}

generated quantities {

}
