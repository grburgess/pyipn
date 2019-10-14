functions {


  matrix[] cos_sin_features_nonstationary(int N, int k, matrix features_one, matrix features_two) {
    /*
      random fourier features based off of
      https://bitbucket.org/flaxter/random-fourier-features-in-stan/src/master/
    */

    // store as an array of matrices so that I can return
    // both at the same time
    matrix[N,k] cos_sin_features[2];

    cos_sin_features[1,:,:] = (cos(features_one)+cos(features_two));
    cos_sin_features[2,:,:] = (sin(features_one)+sin(features_two));
    return cos_sin_features;

  }


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
  vector[N1] time1; // mid-points of LC 1
  vector[N2] time2; // mid-points of LC 2
  vector[N1] exposure1; // mid-points of LC 1
  vector[N2] exposure2; // mid-points of LC 2
  int counts1[N1]; // counts in LC 1
  int counts2[N2]; // counts in LC 2

  vector[3] sc_pos1;
  vector[3] sc_pos2;

  
  // real bw; // the band width of the FFs; lower => smoother
  int<lower=1> k; // number of FFs

  

  // predcition for plotting
  //int N_model;
  //vector[N_model] predict_time;

}
transformed data {

  matrix [N1, k] time_omega1;
  vector[k] omega; // this weird MC integration thing. I suppose I could do this in stan
  


  //  real dt = 29.64;
  real tstart = -5; 
  real tstop = 20;
  real strength = 10.;
  vector[N1] window1 = filter(time1, tstart, tstop  , strength);

  
  for (i in 1:k) {

    omega[i] = normal_rng(0, 1);

  }
  
  time_omega1 = time1 * omega';
  
}
parameters {

  vector[k] beta1;
  vector[k] beta2;



  unit_vector[3] grb_xyz;

  real log_scale;
  real log_bw;
  
  vector[2]  log_bkg;    
  vector[2] log_amplitude; // independent amplitude1 of LC 1; probably do not need right now...

}
transformed parameters {
  vector[N1] fhat1; // raw GP for LC 1
  vector[N2] fhat2; // raw GP for LC 2
  vector[2] bkg = exp(log_bkg);

  real bw = exp(log_bw);
  vector[N1] expected_count1;
  vector[N2] expected_count2;
  
  real scale = exp(log_scale) * inv_sqrt(k);

  
  //  real duration = 10^log_duration;

  //  real tstop = tstart + duration;

  
  real dt = time_delay(grb_xyz, sc_pos1, sc_pos2);


  {

      // do not multiply a matrix by a constant
    matrix [N1, k] features_one1 = bw * time_omega1;
    matrix [N1, k] features_two1 = 0.5 * features_one1;

    matrix[N1,k] tmp[2] = cos_sin_features_nonstationary(N1, k, features_one1, features_two1);
  
      // mulitply by the filter... maybe remove
    fhat1 = window1 .* exp( scale * (tmp[1,:,:] * beta1 + tmp[2,:, :] * beta2) + log_amplitude[1]);
    // fhat1 = exp( scale * (tmp[1,:,:] * beta1 + tmp[2,:, :] * beta2) + log_amplitude[1]);
    


  }
  
  // have to compute the matrices on the fly for the delayed LC
  {

    matrix [N2, k] features_one2 = (bw * (time2 -dt)) * omega';
    matrix [N2, k] features_two2 = 0.5 * features_one2;


    matrix[N2,k] tmp[2] = cos_sin_features_nonstationary(N2, k, features_one2, features_two2);

    fhat2 = filter(time2 - dt, tstart, tstop ,  strength) .* exp( scale * (tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2) + log_amplitude[2] );
    //fhat2 = exp( scale * (tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2) + log_amplitude[2] );
    

  }


  // pre compute for PPC
  
  expected_count1 = exposure1 .* (fhat1 + bkg[1]);
  expected_count2 = exposure2 .* (fhat2 + bkg[2]);


  
}

model {

  // priors

  beta1 ~ std_normal();
  beta2 ~ std_normal();
  
  log_scale ~ std_normal();
  log_bw ~ std_normal();
  
  log_amplitude ~ std_normal();
  log_bkg ~ normal(log(50), 1);  


  //  log_duration ~ normal(1,.2);
  //tstart ~ normal(1,5);

  counts1 ~ poisson( expected_count1 );
  counts2 ~ poisson( expected_count2 );
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
