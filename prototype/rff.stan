functions {

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
    /* for(i in 1:N) */
    /*   for(j in 1:k) { */
    /*  cos_sin_features[1, i,j] = cos(features[i,j]); */
    /*  cos_sin_features[2, i,j] = sin(features[i,j]); */
    /*   } */
    cos_sin_features[1,:,:] = cos(features) * scale;
    cos_sin_features[2,:,:] = sin(features) * scale;
    return cos_sin_features;

  }

  vector filter(vector time, real tstart, real tstop,  real strength) {
    /*

      A two-sided filter that forces the light curve prediction to zero
      outside of a start and stop time. This may not be needed in the future
    */
    return inv_logit(strength * (time - tstart) ) .* (1 - inv_logit(strength * (time- tstop) ));

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

  real bw; // the band width of the FFs; lower => smoother
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
  /* real tstart = -2; */
  /* real tstop = 20; */
  real strength = 20.;

  // for the non-delayed LC, let's go ahead and compute the fucking matrices

  {

    matrix[N1,k] tmp[2] = cos_sin_features(N1, k, time1, omega, bw);

    cosfeatures1 = tmp[1,:,:];
    sinfeatures1 = tmp[2,:,:];

  }

  // sample for the plotting stuffs

  {

    matrix[N_model,k] tmp[2] = cos_sin_features(N_model, k, predict_time, omega, bw);

    predict_cosfeatures = tmp[1,:,:];
    predict_sinfeatures = tmp[2,:,:];

  }

}
parameters {

  vector[k] beta1; // the amplitude along the cos basis
  vector[k] beta2; // the amplitude along the sin basis

  real<lower=0> bkg1; // the bkg for LC 1;  right now this is a constant
  real<lower=0> bkg2; // the bkg for LC 2;  right now this is a constant

  real <lower=-10> tstart;
  real log_duration;
  //real<lower=0, upper=max(time2)> dt; // the time delay

  real<upper=log10(max(time2))> log_dt;

  
  real log_amplitude1; // independent amplitude1 of LC 1; probably do not need right now...
  real log_amplitude2; // independent amplitude1 of LC 2; probably do not need right now...

}
transformed parameters {
  vector[N1] fhat1; // raw GP for LC 1
  vector[N2] fhat2; // raw GP for LC 2

  real dt = 10^log_dt;
  real duration = 10^log_duration;

  real tstop = tstart +duration;

  
  //  real<lower=0> dt = 10^log_dt;

  // mulitply by the filter... maybe remove
  fhat1 = filter(time1, tstart, tstop  , strength) .* exp(cosfeatures1 * beta1 + sinfeatures1*beta2 + log_amplitude1);

  // have to compute the matrices on the fly for the delayed LC
  {

    matrix[N2,k] tmp[2] = cos_sin_features(N2, k, time2 - dt, omega, bw);

    fhat2 = filter(time2 - dt, tstart, tstop ,  strength) .* exp(tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2 + log_amplitude2);

  }

}

model {

  // priors

  beta1 ~ std_normal();
  beta2 ~ std_normal();
  bkg1 ~ normal(50,10);
  bkg2 ~ normal(50,10);

  log_amplitude1 ~ normal(0,1);
  log_amplitude2 ~ normal(0,1);

  log_dt ~ normal(-1,.5);
  log_duration ~ normal(1,1.);
  tstart ~ cauchy(1,5);

  counts1 ~ poisson( exposure1 .* (fhat1 + bkg1));
  counts2 ~ poisson( exposure2 .* (fhat2 + bkg2));
}

generated quantities {

  vector[N_model] predict = filter(predict_time, tstart, tstop , strength) .* exp(predict_cosfeatures * beta1 + predict_sinfeatures * beta2);

  int ppc1[N1];
  int ppc2[N2];

  // PPCs

  for (n in 1:N1) {

    ppc1[n] = poisson_rng( exposure1[n] *(fhat1[n] + bkg1) );

  }

  for (n in 1:N2) {

    ppc2[n] = poisson_rng( exposure2[n] *( fhat2[n] + bkg2) );

  }

}
