functions {

  matrix[] cos_sin_features_nonstationary(int N, int k, vector time, vector omega, real bw) {
    /*
      random fourier features based off of
      https://bitbucket.org/flaxter/random-fourier-features-in-stan/src/master/
    */

    // store as an array of matrices so that I can return
    // both at the same time
    matrix[N,k] cos_sin_features[2];

    matrix[N,k] features_one;
    matrix[N,k] features_two;

    features_one = (bw * time ) * omega';
    features_two = 0.5 * features_one;

    cos_sin_features[1,:,:] = (cos(features_one)+cos(features_two));
    cos_sin_features[2,:,:] = (sin(features_one)+sin(features_two));
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

    matrix[N,k] features;

    features = time * omega' * bw;


    cos_sin_features[1,:,:] = cos(features);
    cos_sin_features[2,:,:] = sin(features);
    return cos_sin_features;

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
  vector[N1] time1; // mid-points of LC 1
  vector[N2] time2; // mid-points of LC 2
  vector[N1] exposure1; // mid-points of LC 1
  vector[N2] exposure2; // mid-points of LC 2
  int counts1[N1]; // counts in LC 1
  int counts2[N2]; // counts in LC 2

  vector[3] sc_pos1;
  vector[3] sc_pos2;

  
  real bw; // the band width of the FFs; lower => smoother
  int<lower=1> k; // number of FFs

  int<lower=1> k2;

  // predcition for plotting
  //int N_model;
  //  vector[N_model] predict_time;

}
transformed data {
  matrix[N1,k] cosfeatures1;
  matrix[N1,k] sinfeatures1;

  matrix[N1,k2] cosfeatures_onoff1;
  matrix[N1,k2] sinfeatures_onoff1;

  real scale = 2.0 * inv_sqrt(k);
  real scale2 = 2.0 * inv_sqrt(k2);
  vector[k] omega; // this weird MC integration thing. I suppose I could do this in stan
  vector[k2] omega2; // this weird MC integration thing. I suppose I could do this in stan

  for (i in 1:k) {

    omega[i] =  normal_rng(0,1);

  }

   for (i in 1:k2) {

    omega2[i] =  normal_rng(0,1);

  }


  {
    
    matrix[N1,k] tmp[2] = cos_sin_features_nonstationary(N1, k, time1, omega, bw);
    
    cosfeatures1 = tmp[1,:,:];
    sinfeatures1 = tmp[2,:,:];

    
  }

  

  // for the non-delayed LC, let's go ahead and compute the fucking matrices

  {

    matrix[N1,k2] tmp[2] = cos_sin_features(N1, k2, time1, omega2, 1.);

    cosfeatures_onoff1 = tmp[1,:,:];
    sinfeatures_onoff1 = tmp[2,:,:];

  }


}
parameters {

  vector[k] beta1; // the amplitude along the cos basis
  vector[k] beta2; // the amplitude along the sin basis

  vector[k2] beta1_onoff; // the amplitude along the cos basis
  vector[k2] beta2_onoff; // the amplitude along the sin basis

  
  real<lower=0> bkg1; // the bkg for LC 1;  right now this is a constant
  real<lower=0> bkg2; // the bkg for LC 2;  right now this is a constant

  unit_vector[3] grb_xyz;

  //  real <lower=-10> tstart;
  //  real<lower=-3, upper=2> log_duration;
  //real<lower=0, upper=max(time2)> dt; // the time delay

  
  real log_amplitude1; // independent amplitude1 of LC 1; probably do not need right now...
  real log_amplitude2; // independent amplitude1 of LC 2; probably do not need right now...


  
}
transformed parameters {
  vector[N1] fhat1; // raw GP for LC 1
  vector[N2] fhat2; // raw GP for LC 2
  vector[N1] grb_window1;
  vector[N2] grb_window2;
  
  real dt = time_delay(grb_xyz, sc_pos1, sc_pos2);

  fhat1 = exp((cosfeatures1*beta1 + sinfeatures1*beta2) * scale + log_amplitude1);  
  // mulitply by the filter... maybe remove

  grb_window1 = inv_logit(-3 + (cosfeatures_onoff1 * beta1_onoff + sinfeatures_onoff1 * beta2_onoff) * scale2 );

  {
    
    matrix[N2, k2] tmp[2] = cos_sin_features(N2, k2, time2 - dt, omega2, 1.);
    
    grb_window2 = inv_logit(-3 + (tmp[1,:,:] * beta1_onoff + tmp[2,:,:] * beta2_onoff) * scale2 ); 
    
  }
  
  // have to compute the matrices on the fly for the delayed LC
  {
    
    matrix[N2,k] tmp[2] = cos_sin_features_nonstationary(N2, k, time2 - dt, omega, bw);
    
    fhat2 = exp((tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2) * scale  + log_amplitude2);
    
  }

}

model {

  // priors

  vector[N1] expected_counts1;
  vector[N2] expected_counts2;
  
  beta1 ~ std_normal();
  beta2 ~ std_normal();
  beta1_onoff ~ std_normal();
  beta2_onoff ~ std_normal();
  bkg1 ~ normal(50,10);
  bkg2 ~ normal(50,10);

  //log_bw ~normal(-1,.5);
  log_amplitude1 ~ normal(0,1);
  log_amplitude2 ~ normal(0,1);


  //  log_duration ~ normal(1,.2);
  //tstart ~ normal(1,5);

  expected_counts1 = exposure1 .* (fhat1 + bkg1);
  expected_counts2 = exposure2 .* (fhat2 + bkg2);

  for (n in 1:N1) {
    target += log_mix(grb_window1[n],
		      poisson_lpmf(counts1[n]| expected_counts1[n]),
		      poisson_lpmf(counts1[n]| bkg1) );

  }

  for (n in 1:N2) {
    target += log_mix(grb_window2[n],
		      poisson_lpmf(counts2[n]| expected_counts2[n]),
		      poisson_lpmf(counts2[n]| bkg2) );

  }

}

generated quantities {

  /* // vector[N_model] predict = filter(predict_time, tstart, tstop , strength) .* exp(predict_cosfeatures * beta1 + predict_sinfeatures * beta2); */
  /* vector[N_model] predict =  exp(predict_cosfeatures * beta1 + predict_sinfeatures * beta2); */

  real grb_phi = atan2(grb_xyz[2], grb_xyz[1]);
  real grb_theta = -(acos(grb_xyz[3]) -0.5 * pi());

  
  /* int ppc1[N1]; */
  /* int ppc2[N2]; */

  /* // PPCs */

  /* for (n in 1:N1) { */

  /*   ppc1[n] = poisson_rng( exposure1[n] *(fhat1[n] + bkg1) ); */

  /* } */

  /* for (n in 1:N2) { */

  /*   ppc2[n] = poisson_rng( exposure2[n] *( fhat2[n] + bkg2) ); */

  /* } */

}
