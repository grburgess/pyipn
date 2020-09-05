real partial_log_like(int[] counts_slice, int start, int end, vector time, vector exposure, row_vector omega1, row_vector omega2, vector beta1, vector beta2, real bw, real dt, real bkg, real scale, real amplitude) {

  int N = size(counts_slice);

  vector[N] time_slice = time[start:end] - dt;
  vector[N] expected_counts_log;
  vector[N] expected_counts;


  expected_counts_log = ((cos(bw * time_slice * omega1) + cos(bw * time_slice * omega2)  ) * beta1)+ ((sin( bw * time_slice * omega1) + sin(bw * time_slice * omega2)  ) * beta2);

  expected_counts = exposure[start:end] .* (exp(scale * expected_counts_log) * amplitude + bkg);


  return poisson_propto_lpmf(counts_slice | expected_counts);

}





real partial_log_like_bw(int[] counts_slice, int start, int end, vector time, vector exposure, row_vector omega1, row_vector omega2, vector beta1, vector beta2, real dt, real bkg, real scale, real amplitude) {

  int N = size(counts_slice);

  vector[N] time_slice = time[start:end] - dt;
  vector[N] expected_counts_log;
  vector[N] expected_counts;

  expected_counts_log = ((cos( time_slice * omega1) + cos( time_slice * omega2)  ) * beta1) + ((sin( time_slice * omega1) + sin( time_slice * omega2)  ) * beta2);

  expected_counts = exposure[start:end] .* (exp(scale * expected_counts_log) * amplitude + bkg);


  return poisson_propto_lpmf(counts_slice | expected_counts);

}

real partial_log_like_bw_multi_scale(int[] counts_slice, int start, int end, vector time, vector exposure, row_vector omega1, row_vector omega2, vector beta1, vector beta2, real dt, real bkg, real scale1, real scale2, real amplitude, int k) {

  int N = size(counts_slice);

  vector[N] time_slice = time[start:end] - dt;
  vector[N] expected_counts_log;
  vector[N] expected_counts;

  matrix [N,k] tw1 = time_slice * omega1;
  matrix [N,k] tw2 = time_slice * omega2;

  
  expected_counts_log = ((scale1 * cos(tw1) + scale2 * cos(tw2)  ) * beta1) + ((scale1 * sin(tw1) + scale2 * sin(tw2)  ) * beta2);

  expected_counts = exposure[start:end] .* (exp(expected_counts_log) * amplitude + bkg);


  return poisson_propto_lpmf(counts_slice | expected_counts);

}

real partial_log_like_bw_multi_scale_fast(int[] counts_slice, int start, int end, vector time, vector exposure, row_vector omega1, row_vector omega2, vector beta1, vector beta2, real dt, real bkg, real scale1, real scale2, real amplitude, int k) {

  int N = size(counts_slice);

  vector[N] time_slice = time[start:end] - dt;
  
  matrix [N,k] tw1 = time_slice * omega1;
  matrix [N,k] tw2 = time_slice * omega2;


  return poisson_propto_lpmf(counts_slice | exposure[start:end] .* (exp(((scale1 * cos(tw1) + scale2 * cos(tw2)  ) * beta1) + ((scale1 * sin(tw1) + scale2 * sin(tw2)  ) * beta2)) * amplitude + bkg));

}



real partial_log_like_bw_multi_scale_log(int[] counts_slice, int start, int end, vector time, vector log_exposure, row_vector omega1, row_vector omega2, vector beta1, vector beta2, real dt, real log_bkg, real scale1, real scale2, real log_amplitude, int k) {

  int N = size(counts_slice);

  vector[N] time_slice = time[start:end] - dt;
  vector[N] expected_rate_log;
  vector[N] expected_counts_log;
  vector[N] log_exposure_slice = log_exposure[start:end];
  
  
  matrix [N,k] tw1 = time_slice * omega1;
  matrix [N,k] tw2 = time_slice * omega2;

  
  expected_rate_log = ((scale1 * cos(tw1) + scale2 * cos(tw2)  ) * beta1) + ((scale1 * sin(tw1) + scale2 * sin(tw2)  ) * beta2);

  /* expected_rate_log = ((scale1 * cos(time_slice * omega1) + scale2 * cos(time_slice * omega2)  ) * beta1) + ((scale1 * sin(time_slice * omega1) + scale2 * sin(time_slice * omega2)  ) * beta2); */


  for (n in 1:N) {
    
    expected_counts_log[n] = log_sum_exp(expected_rate_log[n] + log_amplitude, log_bkg ) + log_exposure_slice[n];
    
  }
  


  return poisson_log_propto_lpmf(counts_slice | expected_counts_log);

}



real angular_separation(vector grb_xyz, vector sc_pointing_norm) {

  real tmp;

  tmp = dot_product(grb_xyz, sc_pointing_norm) ;

  if (tmp >0) {
  
    return tmp;  
    }

  else {

    return 0;
  }

}



real calculate_horizon_angle(vector sc_position) {

  real earth_radius = 6371.0;
  real altitude = sqrt(sum((sc_position .* sc_position)));
  real horizon_angle = 0.5 * pi() - acos(earth_radius/altitude);

  return horizon_angle;


  

}

real earth_occulation(real horizon_angle, vector sc_position, vector grb_xyz) {

  real tmp;

  real angle = dot_product(grb_xyz, -sc_pointing_norm);

  if (angle < horizon_angle) {

      return 0.
    }

  else {

    return 1.

  }

  

}




real time_delay( vector grb_xyz, vector sc_pos1, vector sc_pos2) {
  // compute the time delay between the signals recieved from the GRB

  real c = 299792.46; // km/s

  real t1 = dot_product(grb_xyz, sc_pos1);
  real t2 = dot_product(grb_xyz, sc_pos2);


  return (t1 - t2)/c;





}
