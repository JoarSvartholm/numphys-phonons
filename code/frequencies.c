#include <math.h>
#include <gsl/gsl_eigen.h>

int map(int i){
  if(i<3){return i;}
  else{return i-3;}
}

void frequencies(double A, double B, double m, double *q, double *omega, double *eps){

  int N = 3;
  double Mii,Mi,omega2;

  gsl_vector *eval = gsl_vector_alloc(N);
  gsl_matrix *evec = gsl_matrix_alloc(N,N);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(N);
  gsl_matrix *M = gsl_matrix_calloc(N,N);

  /*Construct matrix*/
  for(int i = 0; i<N;i++){
    Mii = 1/(2*m)*(A+B)*(8-4*cos(M_PI*q[map(i)])*cos(M_PI*q[map(i+1)])-4*cos(M_PI*q[map(i)])*cos(M_PI*q[map(i+2)]))+B/m*(4-4*cos(M_PI*q[map(i+1)])*cos(M_PI*q[map(i+2)]));
    gsl_matrix_set(M,i,i,Mii);
  }
  for(int i = 0; i<N-1;i++){
    Mi = 1/(2*m)*(A-B)*4*sin(M_PI*q[map(i+1)])*sin(M_PI*q[map(i)]);
    gsl_matrix_set(M,i+1,i,Mi);
  }
  Mi = 1/(2*m)*(A-B)*4*sin(M_PI*q[2])*sin(M_PI*q[0]);
  gsl_matrix_set(M,2,0,Mi);

  /*solve eigenvalue problem*/
  gsl_eigen_symmv(M,eval,evec,w);
  gsl_eigen_symmv_free(w);

  gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_ASC);

  /*save output*/
  for(int i=0;i<N;i++){
    omega2 = gsl_vector_get(eval,i);
    if(omega2<0) omega2=0;
    omega[i] = sqrt(omega2);

    if(eps != NULL){
      eps[i]=gsl_matrix_get(evec,i,0);
      eps[i+3]=gsl_matrix_get(evec,i,1);
      eps[i+6]=gsl_matrix_get(evec,i,2);
    }
  }

}
