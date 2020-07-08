functions {
#include functions.stan

}


data {
  int<lower=1> N1;
  int<lower=1> N2;
  vector[N1] time1;
  vector[N2] time2;
  int counts1[N1];
  int counts2[N2];
  vector[N1] exposure1; // exposure of LC 1
  vector[N2] exposure2; // exposure of LC 2
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
  
  //  real dt = 29.64/2.;
  real tstart = 0;
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

  vector[2]  log_bkg;

  real log_scale;
  
  unit_vector[3] grb_xyz;

  vector[2] log_amplitude;
  
}
transformed parameters {
  vector[N1] fhat1;
  vector[N2] fhat2;

  //  real<lower=0> dt = 10^log_dt;

  fhat1 = exp(cosfeatures1 * beta1 + sinfeatures1*beta2);
  

 

  
  {

    matrix[N2,k] tmp[2] = cos_sin_features(N2, k, time2 - dt, omega, bw);

    fhat2 = exp(tmp[1,:,:] * beta1 + tmp[2,:,:] * beta2);

  }


}

model {

  beta1 ~ std_normal();
  beta2 ~ std_normal();
  bkg1 ~ normal(50,10);
  bkg2 ~ normal(50,10);

  
    dt ~ normal(30,10);
  
  counts1 ~ poisson( fhat1 + bkg1);
  counts2 ~ poisson( fhat2 + bkg2);
}



generated quantities {


  vector[N_model] predict =  exp(predict_cosfeatures * beta1 + predict_sinfeatures * beta2);


  
  
  

}
