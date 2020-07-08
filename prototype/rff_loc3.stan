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




}
transformed data {

  real bw = 1.;
  
  matrix [N1, k] time_omega1;
  matrix [N1, k] time_omega2;
  vector[k] omega1; // this weird MC integration thing. I suppose I could do this in stan
  vector[k] omega2; // this weird MC integration thing. I suppose I could do this in stan


  //  real dt = 29.64;
  real tstart = -5;
  real tstop = 20;
  real strength = 10.;
  vector[N1] window1 = filter(time1, tstart, tstop  , strength);

  // for the non-delayed LC, let's go ahead and compute the fucking matrices

  for (i in 1:k) {

    omega1[i] = normal_rng(0, 1);
    omega2[i] = normal_rng(0, 1);
  }

  time_omega1 = time1 * omega1';
  time_omega2 = time1 * omega2';

  
}
parameters {

  vector[k] beta1; // the amplitude along the cos basis
  vector[k] beta2; // the amplitude along the sin basis

  vector[3]  log_bkg;


  real log_scale;

 

  unit_vector[3] grb_xyz;

  vector[3] log_amplitude; // independent amplitude1 of LC 1; probably do not need right now...
}

transformed parameters {
  vector[N1] fhat1; // raw GP for LC 1
  vector[N2] fhat2; // raw GP for LC 2
  vector[N3] fhat3; // raw GP for LC 3

  vector[3] bkg = exp(log_bkg);


  //  real bw = exp(log_bw);
  vector[N1] expected_count1;
  vector[N2] expected_count2;
  vector[N3] expected_count3;

  real scale = exp(log_scale) * inv_sqrt(k);


  real dt_1_2 = time_delay(grb_xyz, sc_pos1, sc_pos2);
  real dt_1_3 = time_delay(grb_xyz, sc_pos1, sc_pos3);
  real dt_2_3 = time_delay(grb_xyz, sc_pos2, sc_pos3);

  {

      // do not multiply a matrix by a constant
    matrix [N1, k] features_one1 = bw * time_omega1;
    matrix [N1, k] features_two1 = bw * time_omega2;

    matrix[N1,k] tmp[2] = cos_sin_features_nonstationary(N1, k, features_one1, features_two1);

    // mulitply by the filter... maybe remove
    fhat1 = window1 .* exp( scale * (tmp[1,:,:] * beta1 + tmp[2,:, :] * beta2) + log_amplitude[1]);



  }

  // have to compute the matrices on the fly for the delayed LC
  {

    matrix [N2, k] features_one2 = (bw * (time2 - dt_1_2)) * omega1';
    matrix [N2, k] features_two2 = (bw * (time2 - dt_1_2)) * omega2';
    


    matrix[N2,k] tmp[2] = cos_sin_features_nonstationary(N2, k, features_one2, features_two2);

    fhat2 = filter(time2 - dt_1_2, tstart, tstop ,  strength) .* exp( scale * (tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2) + log_amplitude[2] );
    //fhat2 =  exp(tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2 + log_amplitude2);

  }


  // have to compute the matrices on the fly for the delayed LC
  {

    matrix [N3, k] features_one3 = (bw * (time3 - dt_1_3)) * omega1';
    matrix [N3, k] features_two3 = (bw * (time3 - dt_1_3)) * omega2';
    
    matrix[N3,k] tmp[2] = cos_sin_features_nonstationary(N3, k, features_one3, features_two3);

    fhat3 = filter(time3 - dt_1_3, tstart, tstop ,  strength) .* exp( scale * (tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2) + log_amplitude[2] );


  }


  expected_count1 = exposure1 .* (fhat1 + bkg[1]);
  expected_count2 = exposure2 .* (fhat2 + bkg[2]);
  expected_count3 = exposure3 .* (fhat3 + bkg[3]);





}

model {

  // priors

  beta1 ~ std_normal();
  beta2 ~ std_normal();

  log_scale ~ std_normal();

  log_bkg ~ normal(log(50), 1);

  // log_bw ~ std_normal();

  log_amplitude ~ std_normal();



  //  log_duration ~ normal(1,.2);
  //tstart ~ normal(1,5);

  counts1 ~ poisson( expected_count1 );
  counts2 ~ poisson( expected_count2 );
  counts3 ~ poisson( expected_count3 );


}

generated quantities {

  real grb_phi = atan2(grb_xyz[2], grb_xyz[1]);
  real grb_theta = -( acos(grb_xyz[3]) - 0.5*pi());


  int ppc1[N1];
  int ppc2[N2];
  int ppc3[N3];

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
