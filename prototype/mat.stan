transformed data {

  int N=100;
  int k=10;
  
  matrix [N, k] mat[2];

  matrix[2, k] beta1;

  vector [N] out = mat * beta1;


}
