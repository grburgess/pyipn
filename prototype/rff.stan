functions {

  matrix[] cos_sin_features(int N, int k, vector time, vector omega, real bw) {
    matrix[N,k] cos_sin_features[2];
    real scale;
    matrix[N,k] features;
  
    features = time * omega' * bw;
  
    scale = sqrt(2.0/k);
    for(i in 1:N)
      for(j in 1:k) {
	cos_sin_features[1, i,j] = cos(features[i,j]);
	cos_sin_features[2, i,j] = sin(features[i,j]);
      }
    cos_sin_features[1,:,:] = cos_sin_features[1,:,:] * scale;
    cos_sin_features[2,:,:] = cos_sin_features[2,:,:] * scale;
    return cos_sin_features;




  }

  vector filter(vector time, real tstart, real tstop,  real strength) {
  
    return inv_logit(strength * (time - tstart) ) .* (1 - inv_logit(strength * (time- tstop) ));
  
}
  


}


data {
  int<lower=1> N1;
  int<lower=1> N2;
  vector[N1] time1;
  vector[N2] time2;
  int counts1[N1];
  int counts2[N2];

  real bw;
  int<lower=1> k;
  
  vector[k] omega;

  int N_model;
  vector[N_model] predict_time;

  

  
}
transformed data {
  matrix[N1,k] cosfeatures1;
  matrix[N1,k] sinfeatures1;

  matrix[N_model,k] predict_cosfeatures;
  matrix[N_model,k] predict_sinfeatures;
  
  //  real dt = 29.64;
  real tstart = -1;
  real tstop = 10;

  {

    matrix[N1,k] tmp[2] = cos_sin_features(N1, k, time1, omega, bw);

    cosfeatures1 = tmp[1,:,:];
    sinfeatures1 = tmp[2,:,:];


  }


    {

    matrix[N_model,k] tmp[2] = cos_sin_features(N_model, k, predict_time, omega, bw);

    predict_cosfeatures = tmp[1,:,:];
    predict_sinfeatures = tmp[2,:,:];


  }

  
}
parameters {
  vector[k] beta1;
  vector[k] beta2;

  real<lower=0> bkg1;
  real<lower=0> bkg2;

  real<lower=0> dt;

  real log_amplitude1;
  real log_amplitude2;
  
  
}
transformed parameters {
  vector[N1] fhat1;
  vector[N2] fhat2;

  //  real<lower=0> dt = 10^log_dt;

  fhat1 = filter(time1, tstart, tstop  , 100.) .* exp(cosfeatures1 * beta1 + sinfeatures1*beta2 + log_amplitude1);
  

 

  
  {

    matrix[N2,k] tmp[2] = cos_sin_features(N2, k, time2 - dt, omega, bw);

    fhat2 = filter(time2 - dt, tstart, tstop ,  100.) .* exp(tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2 + log_amplitude2);

  }


}

model {

  beta1 ~ std_normal();
  beta2 ~ std_normal();
  bkg1 ~ normal(50,10);
  bkg2 ~ normal(50,10);

  log_amplitude1 ~ normal(0,1);
  log_amplitude2 ~ normal(0,1);
  
  dt ~ normal(30,10);
  
  counts1 ~ poisson( fhat1 + bkg1);
  counts2 ~ poisson( fhat2 + bkg2);
}



generated quantities {


  vector[N_model] predict = filter(predict_time, tstart, tstop , 100) .* exp(predict_cosfeatures * beta1 + predict_sinfeatures * beta2);

  int ppc1[N1];
  int ppc2[N2];

  for (n in 1:N1) {

    ppc1[n] = poisson_rng( fhat1[n] + bkg1 );

  }
  
  for (n in 1:N2) {

    ppc2[n] = poisson_rng( fhat2[n] + bkg2 );

  }

  
  
  

}
