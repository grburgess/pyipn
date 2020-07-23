

real partial_log_like(int[] counts_slice, int start, int end, vector time, vector exposure, row_vector omega1, row_vector omega2, vector beta1, vector beta2, real bw, real dt, real bkg, real scale, real amplitude) {

  int N = size(counts_slice);

  vector[N] time_slice = time[start:end] - dt;
  vector[N] expected_counts_log;
  vector[N] expected_counts;

  
  expected_counts_log = ((cos(bw * time_slice * omega1) + cos(bw * time_slice * omega2)  ) * beta1)+ ((sin( bw * time_slice * omega1) + sin(bw * time_slice * omega2)  ) * beta2);

  expected_counts = exposure[start:end] .* (exp(scale * expected_counts_log) * amplitude + bkg);


  return poisson_lpmf(counts_slice | expected_counts);

}




real partial_log_like_bw(int[] counts_slice, int start, int end, vector time, vector exposure, row_vector omega1, row_vector omega2, vector beta1, vector beta2, real dt, real bkg, real scale, real amplitude) {

  int N = size(counts_slice);

  vector[N] time_slice = time[start:end] - dt;
  vector[N] expected_counts_log;
  vector[N] expected_counts;

  expected_counts_log = ((cos( time_slice * omega1) + cos( time_slice * omega2)  ) * beta1) + ((sin( time_slice * omega1) + sin( time_slice * omega2)  ) * beta2);

  expected_counts = exposure[start:end] .* (exp(scale * expected_counts_log) * amplitude + bkg);


  return poisson_lpmf(counts_slice | expected_counts);

}






real time_delay( vector grb_xyz, vector sc_pos1, vector sc_pos2) {
  // compute the time delay between the signals recieved from the GRB

  real c = 299792.46; // km/s

  real t1 = dot_product(grb_xyz, sc_pos1);
  real t2 = dot_product(grb_xyz, sc_pos2);


  return (t1 - t2)/c;





}
