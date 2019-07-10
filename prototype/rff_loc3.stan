functions {

  matrix[] cos_sin_features_nonstationary(int N, int k, vector time, vector omega, real bw_one, real bw_two) {
    /*
      random fourier features based off of
      https://bitbucket.org/flaxter/random-fourier-features-in-stan/src/master/
    */

    // store as an array of matrices so that I can return
    // both at the same time
    matrix[N,k] cos_sin_features[2];
    real scale;
    matrix[N,k] features_one;
    matrix[N,k] features_two;

    features_one = time * omega' * bw_one;
    features_two = time * omega' * bw_two;

    scale = sqrt(2.0/k);
    cos_sin_features[1,:,:] = (cos(features_one)+cos(features_two)) * scale;
    cos_sin_features[2,:,:] = (sin(features_one)+sin(features_two)) * scale;
    return cos_sin_features;

  }

  
  matrix[] cos_sin_features(int N, int k, vector time, vector omega, real bw) {
    /*
      random fourier features based off of
      https://bitbucket.org/flaxter/random-fourier-features-in-stan/src/master/

    */

    // store as an array of matrices so that I can return
    // both at the same time
    matrix[N,k] cos_sin_features[2];
    real scale;
    matrix[N,k] features;

    features = time * omega' * bw;

    scale = sqrt(2.0/k);
    cos_sin_features[1,:,:] = cos(features) * scale;
    cos_sin_features[2,:,:] = sin(features) * scale;
    return cos_sin_features;

  }

  /* vector filter(vector time, real tstart, real tstop,  real strength) { */
  /*   /\* */

  /*     A two-sided filter that forces the light curve prediction to zero */
  /*     outside of a start and stop time. This may not be needed in the future */
  /*   *\/ */
  /*   return inv_logit(strength * (time - tstart) ) .* (1 - inv_logit(strength * (time- tstop) )); */

  /* } */



  vector filter(vector time, real tstart, real tstop,  real strength) {
    /*

      A two-sided filter that forces the light curve prediction to zero
      outside of a start and stop time. This may not be needed in the future
    */
    return 0.25 * (1+tanh(strength * (time - tstart) )) .* (1 - tanh(strength * (time- tstop) ));

  }


  

  real time_delay( vector grb_xyz, vector sc_pos1, vector sc_pos2) {


    real c = 299792.46; // km/s

    real t1 = dot_product(grb_xyz, sc_pos1);
    real t2 = dot_product(grb_xyz, sc_pos2);


    return (t1 - t2)/c;
    
    



  }
  
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

  
  real bw; // the band width of the FFs; lower => smoother
  real bw2; // the band width of the FFs; lower => smoother
  int<lower=1> k; // number of FFs

  vector[k] omega; // this weird MC integration thing. I suppose I could do this in stan

  // predcition for plotting
  int N_model;
  vector[N_model] predict_time;

}
transformed data {
  matrix[N1,k] cosfeatures1;
  matrix[N1,k] sinfeatures1;

  matrix[N_model,k] predict_cosfeatures;
  matrix[N_model,k] predict_sinfeatures;

  //  real dt = 29.64;
  real tstart = -5; 
  real tstop = 20;
  real strength = 20.;

  // for the non-delayed LC, let's go ahead and compute the fucking matrices

  {

    matrix[N1,k] tmp[2] = cos_sin_features_nonstationary(N1, k, time1, omega, bw, bw2);

    cosfeatures1 = tmp[1,:,:];
    sinfeatures1 = tmp[2,:,:];

  }

  // sample for the plotting stuffs

  {

    matrix[N_model,k] tmp[2] = cos_sin_features_nonstationary(N_model, k, predict_time, omega, bw, bw2);

    predict_cosfeatures = tmp[1,:,:];
    predict_sinfeatures = tmp[2,:,:];

  }

}
parameters {

  vector[k] beta1; // the amplitude along the cos basis
  vector[k] beta2; // the amplitude along the sin basis

  real<lower=0> bkg1; // the bkg for LC 1;  right now this is a constant
  real<lower=0> bkg2; // the bkg for LC 2;  right now this is a constant
  real<lower=0> bkg3; // the bkg for LC 2;  right now this is a constant

  unit_vector[3] grb_xyz;

  //  real <lower=-10> tstart;
  //  real<lower=-3, upper=2> log_duration;
  //real<lower=0, upper=max(time2)> dt; // the time delay

  
  real log_amplitude1; // independent amplitude1 of LC 1; probably do not need right now...
  real log_amplitude2; // independent amplitude1 of LC 2; probably do not need right now...
  real log_amplitude3; // independent amplitude1 of LC 2; probably do not need right now...

}
transformed parameters {
  vector[N1] fhat1; // raw GP for LC 1
  vector[N2] fhat2; // raw GP for LC 2
  vector[N3] fhat3; // raw GP for LC 3

  real<lower = -pi(), upper = pi()> grb_phi = atan2(grb_xyz[2], grb_xyz[1]);
  real<lower=0, upper = pi()> grb_theta = acos(grb_xyz[3]);

  
  //  real duration = 10^log_duration;

  //  real tstop = tstart + duration;

  
  real dt_1_2 = time_delay(grb_xyz, sc_pos1, sc_pos2);
  real dt_1_3 = time_delay(grb_xyz, sc_pos1, sc_pos3);
  real dt_2_3 = time_delay(grb_xyz, sc_pos2, sc_pos3);

  // mulitply by the filter... maybe remove
  fhat1 = filter(time1, tstart, tstop  , strength) .* exp(cosfeatures1 * beta1 + sinfeatures1*beta2 + log_amplitude1);
  //fhat1 =  exp(cosfeatures1 * beta1 + sinfeatures1*beta2 + log_amplitude1);

  // have to compute the matrices on the fly for the delayed LC
  {

    matrix[N2,k] tmp[2] = cos_sin_features_nonstationary(N2, k, time2 - dt_1_2, omega, bw, bw2);

    fhat2 = filter(time2 - dt_1_2, tstart, tstop ,  strength) .* exp(tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2 + log_amplitude2);
    //fhat2 = exp(tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2 + log_amplitude2);

  }

    // have to compute the matrices on the fly for the delayed LC
  {

    matrix[N3,k] tmp[2] = cos_sin_features_nonstationary(N3, k, time3 - dt_1_3, omega, bw, bw2);

    fhat3 = filter(time3 - dt_1_3, tstart, tstop ,  strength) .* exp(tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2 + log_amplitude3);
    //    fhat3 =  exp(tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2 + log_amplitude3);

  }


  
}

model {

  // priors

  beta1 ~ std_normal();
  beta2 ~ std_normal();
  bkg1 ~ normal(50,10);
  bkg2 ~ normal(50,10);
  bkg3 ~ normal(50,10);

  log_amplitude1 ~ normal(0,1);
  log_amplitude2 ~ normal(0,1);
  log_amplitude3 ~ normal(0,1);


  //  log_duration ~ normal(1,.2);
  //tstart ~ normal(1,5);

  counts1 ~ poisson( exposure1 .* (fhat1 + bkg1));
  counts2 ~ poisson( exposure2 .* (fhat2 + bkg2));
  counts3 ~ poisson( exposure3 .* (fhat3 + bkg3));
}

generated quantities {

  /* //  vector[N_model] predict = filter(predict_time, tstart, tstop , strength) .* exp(predict_cosfeatures * beta1 + predict_sinfeatures * beta2); */
  /* vector[N_model] predict =  exp(predict_cosfeatures * beta1 + predict_sinfeatures * beta2); */

  int ppc1[N1];
  int ppc2[N2];
  int ppc3[N3];

  // PPCs

  for (n in 1:N1) {

    ppc1[n] = poisson_rng( exposure1[n] *(fhat1[n] + bkg1) );

  }

  for (n in 1:N2) {

    ppc2[n] = poisson_rng( exposure2[n] *( fhat2[n] + bkg2) );

  }


  for (n in 1:N3) {

    ppc3[n] = poisson_rng( exposure3[n] *( fhat3[n] + bkg3) );

  }

  
}
