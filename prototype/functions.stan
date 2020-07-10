real partial_log_like(int[] counts_slice, int start, int end, vector time, vector exposure, row_vector omega1, row_vector omega2, vector beta1, vector beta2, real bw, real dt, real bkg, real scale, real amplitude) {

    int N = size(counts_slice);
    real lp = 0;

    vector[N] time_slice = time[start:end] - dt;
    vector[N] expected_counts_log;
    vector[N] expected_counts;
    
    expected_counts_log = ((cos(bw * time_slice * omega1) + cos(bw * time_slice * omega2)  ) * beta1) + ((sin(bw * time_slice * omega1) + sin(bw * time_slice * omega2)  ) * beta2);

    expected_counts = exposure[start:end] .* (exp(expected_counts_log) + bkg + amplitude);
    
  /*   for (n in 1:N) { */

  /*   lp += poisson_lpmf(counts_slice[n] | expected_counts[n]); */
  /* } */
    return poisson_lpmf(counts_slice | expected_counts);
    
  }



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



  matrix[] cos_sin_features(int N, int k, vector time, vector omega, real bw) {
    /*
      random fourier features based off of
      https://bitbucket.org/flaxter/random-fourier-features-in-stan/src/master/

    */

    // store as an array of matrices so that I can return
    // both at the same time
    matrix[N,k] cos_sin_features_mat[2];
    real scale;
    matrix[N,k] features;

    features = time * omega' * bw;

    scale = sqrt(2.0/k);
    cos_sin_features_mat[1,:,:] = cos(features) * scale;
    cos_sin_features_mat[2,:,:] = sin(features) * scale;
    return cos_sin_features_mat;

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
