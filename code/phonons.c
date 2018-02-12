#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_eigen.h>
#include "frequencies.h"

typedef struct {
  int substance;
  int test;
  int npoints;
}test_case;

typedef struct {
  double sigma;
  double eps;
  double r;
  double m;
}params;

typedef struct {
  double q1[3];
  double q2[3];
}qvecs;

int min(int A, int B){
  if(A<B) {return A;}
  else return B;
}

void check_input(int argc,char const *argv[],test_case *Case){
  //Test substance input
  if(strcmp(argv[1],"Ne")==0){
    Case->substance=1;
  } else if(strcmp(argv[1],"Ar")==0){
    Case->substance=2;
  } else if(strcmp(argv[1],"Kr")==0){
    Case->substance=3;
  } else if(strcmp(argv[1],"Xe")==0){
    Case->substance=4;
  } else {
    fprintf(stderr, "Error in input. Not a valid substance: %s\n",argv[1] );
    exit(1);}

    //Test type input
    if(strcmp(argv[2],"omega")==0){
      Case->test = 1;
    }else if(strcmp(argv[2],"gamma")==0){
      Case->test = 2;
    } else if(strcmp(argv[2],"cv")==0){
      Case->test = 3;
    } else{
      fprintf(stderr, "Error in input. Not a valid test case: %s\n", argv[2]);
      exit(1);
    }
}

void save_qvecs(int argc, char const *argv[],test_case *Case, qvecs *Q){
  if(Case->test == 1 || Case->test == 2){
    if(argc<6 || (argc>6 && argc<9)){
      fprintf(stderr, "Error in input: not enough input parameters\n");
      exit(1);
    }
    if(argc > 10){
      fprintf(stderr, "Error in input: Too many input parameters\n");
      exit(1);
    }
    if(argc == 10 && atoi(argv[9]) < 2){
      fprintf(stderr, "Error in input: Cannot compute at %s points\n",argv[9]);
      exit(1);
    }

    for(int i=3;i<6;i++){
      Q->q1[i-3] = atof(argv[i]);
    }
    if(argc>6){
      for(int i=6;i<9;i++){
        Q->q2[i-6] = atof(argv[i]);
      }
    }
    if(argc ==10){
       Case->npoints = atoi(argv[9]);
     } else if(argc == 6){
       Case->npoints = 1;
     } else Case->npoints = 11;

  }
}

void set_params(int substance, params *par){
  double A = pow(10,-10);
  switch (substance) {
    case 1:
      par->sigma = 3.035*A;
      par->eps = 0.0721*A*A;
      par->r = 3.1562*A;
      par->m = 0.335092 *pow(10,-25);
      break;
    case 2:
      par->sigma = 3.709*A;
      par->eps = 0.236*A*A;
      par->r =  3.7477*A;
      par->m =  0.66335*pow(10,-25);
      break;
    case 3:
      par->sigma = 3.966*A;
      par->eps = 0.325*A*A;
      par->r =   3.9922*A;
      par->m =   1.3915*pow(10,-25);
      break;
    case 4:
      par->sigma = 4.318*A;
      par->eps = 0.458*A*A;
      par->r =    4.3346*A;
      par->m =    2.18017*pow(10,-25);
      break;
  }
}


double compute_A(params *pars){
  double r14 = pow(pars->r,-14);
  double r8 = pow(pars->r,-8);
  return 12*pars->eps*(13*pow(pars->sigma,11)*r14-7*pow(pars->sigma,5)*r8);
}

double compute_B(params *pars){
  double r14 = pow(pars->r,-14);
  double r8 = pow(pars->r,-8);
  return 12*pars->eps*(pow(pars->sigma,5)*r8 - pow(pars->sigma,11)*r14);
}


int main(int argc, char const *argv[]) {

  char fname[64];

  if(argc>9){
    snprintf(fname,64,"../data/phonons_%s_%s_%s%s%s_n=%s.data",argv[1],argv[2],argv[6],argv[7],argv[8],argv[9]);
  }else{
    snprintf(fname,64,"../data/phonons_%s_%s.data",argv[1],argv[2]);
  }
  FILE *f = fopen(fname,"w+");
  test_case Case;
  params pars;
  qvecs Q;
  double omega[3];
  double eps[9];
  double A,B,q[3];

  check_input(argc,argv,&Case);

  set_params(Case.substance,&pars);

  save_qvecs(argc,argv,&Case,&Q);
  A= compute_A(&pars);
  B= compute_B(&pars);

  for(int i=0;i<Case.npoints;i++){
    if(Case.npoints==1){
      frequencies(A,B,pars.m,Q.q1,omega,eps);
      printf("%f %f %f %f %f %f\n", Q.q1[0],Q.q1[1],Q.q1[2],omega[0],omega[1],omega[2]);
      fprintf(f,"%f %f %f %f %f %f\n", Q.q1[0],Q.q1[1],Q.q1[2],omega[0],omega[1],omega[2]);
    }else{
      for(int j=0;j<3;j++){
        q[j] = Q.q1[j]+i*(Q.q2[j]-Q.q1[j])/(Case.npoints-1);
      }
      frequencies(A,B,pars.m,q,omega,eps);
      printf("%f %f %f %f %f %f\n", q[0],q[1],q[2],omega[0],omega[1],omega[2]);
      fprintf(f,"%f %f %f %f %f %f\n", q[0],q[1],q[2],omega[0],omega[1],omega[2]);
    }
  }
  fclose(f);

  return 0;
}
