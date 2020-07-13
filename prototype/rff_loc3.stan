functions {
#include functions.stan
}

data {
  int<lower=1> N1; // number of time bins for LC 1
  int<lower=1> N2; // number of time bins for LC 2
  int<lower=1> N3; // number of time bins for LC 3

  vector[N1] time1; // mid-points of LC 1
  vector[N2] time2; // mid-points of LC 2
  vector[N3] time3; // mid-points of LC 3
  vector[N1] exposure1; // mid-points of LC 1
  vector[N2] exposure2; // mid-points of LC 2
  vector[N3] exposure3; // mid-points of LC 3
  int counts1[N1]; // counts in LC 1
  int counts2[N2]; // counts in LC 2
  int counts3[N3]; // counts in LC 3

  vector[3] sc_pos1;
  vector[3] sc_pos2;
  vector[3] sc_pos3;


  /* real bw; // the band width of the FFs; lower => smoother */
  /* real bw2; // the band width of the FFs; lower => smoother; */
  int<lower=1> k; // number of FFs

  int grainsize;



}
transformed data {

  real bw = 1.;

  row_vector[k] omega1; // this weird MC integration thing. I suppose I could do this in stan
  row_vector[k] omega2; // this weird MC integration thing. I suppose I could do this in stan



  // for the non-delayed LC, let's go ahead and compute the fucking matrices

  for (i in 1:k) {

    omega1[i] = normal_rng(0, 1);
    omega2[i] = normal_rng(0, 1);
  }
  
}
parameters {

  vector[k] beta1; // the amplitude along the cos basis
  vector[k] beta2; // the amplitude along the sin basis

  vector[3]  log_bkg;


  //  real log_bw;
  real log_scale;

 

  unit_vector[3] grb_xyz;

  vector[3] log_amplitude; // independent amplitude1 of LC 1; probably do not need right now...
}

transformed parameters {
 
  vector[3] bkg = exp(log_bkg);
  vector[3] amplitude = exp(log_amplitude);

  //  real bw = exp(log_bw);
  real scale = exp(log_scale) * inv_sqrt(k);


  real dt_1_2 = time_delay(grb_xyz, sc_pos1, sc_pos2);
  real dt_1_3 = time_delay(grb_xyz, sc_pos1, sc_pos3);
  real dt_2_3 = time_delay(grb_xyz, sc_pos2, sc_pos3);

  

}

model {

  // priors

  beta1 ~ std_normal();
  beta2 ~ std_normal();

  log_scale ~ std_normal();

  log_bkg ~ normal(log(50), 1);

  //  log_bw ~ normal(-1, 2);

  log_amplitude ~ std_normal();



  //  log_duration ~ normal(1,.2);
  //tstart ~ normal(1,5);


  target += reduce_sum(partial_log_like, counts1, grainsize, time1, exposure1, omega1, omega2, beta1, beta2, bw, 0., bkg[1], scale, amplitude[1]);
  
  target += reduce_sum(partial_log_like, counts2, grainsize, time2, exposure2, omega1, omega2, beta1, beta2, bw, dt_1_2, bkg[2], scale, amplitude[2]);
  
  target += reduce_sum(partial_log_like, counts3, grainsize, time3, exposure3, omega1, omega2, beta1, beta2, bw, dt_1_3, bkg[3], scale, amplitude[3]);

  

}

generated quantities {

  real grb_phi = atan2(grb_xyz[2], grb_xyz[1]);
  real grb_theta = -( acos(grb_xyz[3]) - 0.5*pi());


  /* // PPCs */

  /* for (n in 1:N1) { */

  /*   ppc1[n] = poisson_rng( expected_count1[n]); */


  /* } */

  /* for (n in 1:N2) { */

  /*   ppc2[n] = poisson_rng( expected_count2[n]); */

  /* } */


  /* for (n in 1:N3) { */

  /*   ppc3[n] = poisson_rng( expected_count3[n]); */

  /* } */


}
