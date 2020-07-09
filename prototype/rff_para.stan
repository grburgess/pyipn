functions {
#include functions.stan  
}

data {
  int<lower=1> N1; // number of time bins for LC 1
  int<lower=1> N2; // number of time bins for LC 2
  vector[N1] time1; // mid-points of LC 1
  vector[N2] time2; // mid-points of LC 2
  vector[N1] exposure1; // mid-points of LC 1
  vector[N2] exposure2; // mid-points of LC 2
  int counts1[N1]; // counts in LC 1
  int counts2[N2]; // counts in LC 2

  int grainsize;
  vector[3] sc_pos1;
  vector[3] sc_pos2;

  
  // real bw; // the band width of the FFs; lower => smoother
  int<lower=1> k; // number of FFs

  

  // predcition for plotting
  //int N_model;
  //vector[N_model] predict_time;

}
transformed data {

  real bw = 1.;
  

  row_vector[k] omega1; // this weird MC integration thing. I suppose I could do this in stan
  row_vector[k] omega2; // this weird MC integration thing. I suppose I could do this in stan
  

  
  for (i in 1:k) {

    omega1[i] = normal_rng(0, 1);
    omega2[i] = normal_rng(0, 1);

  }
  
}
parameters {

  vector[k] beta1;
  vector[k] beta2;

  unit_vector[3] grb_xyz;

  real log_scale;
  
  vector[2]  log_bkg;    
  vector[2] log_amplitude; // independent amplitude1 of LC 1; probably do not need right now...

}
transformed parameters {

  vector[2] bkg = exp(log_bkg);
  vector[2] amplitude = exp(log_amplitude);
  
  real scale = exp(log_scale) * inv_sqrt(k);

  
  real dt = time_delay(grb_xyz, sc_pos1, sc_pos2);


  
}

model {

  // priors

  beta1 ~ std_normal();
  beta2 ~ std_normal();
  
  log_scale ~ std_normal();
  
  log_amplitude ~ std_normal();
  log_bkg ~ normal(log(50), 1);  

  

  target += reduce_sum(partial_log_like, counts1, grainsize, time1, exposure1, omega1, omega2, beta1, beta2, 0., bkg[1], scale, amplitude[1]);

  target += reduce_sum(partial_log_like, counts2, grainsize, time2, exposure2, omega1, omega2, beta1, beta2, dt, bkg[2], scale, amplitude[2]);

}

generated quantities {

  real grb_phi = atan2(grb_xyz[2], grb_xyz[1]);
  real grb_theta = -( acos(grb_xyz[3]) - 0.5*pi());

  
  /* int ppc1[N1]; */
  /* int ppc2[N2]; */

  /* // PPCs */

  /* for (n in 1:N1) { */

  /*   ppc1[n] = poisson_rng( expected_count1[n]); */

  /* } */

  /* for (n in 1:N2) { */

  /*   ppc2[n] = poisson_rng( expected_count2[n]); */

  /* } */

}
