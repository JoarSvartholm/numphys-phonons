#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_eigen.h>

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

void check_input(int argc,char const *argv[],test_case *Case){
  //Test substance input
  if(strcmp(argv[1],"Ne")==0){
    printf(" Ne funkar\n" );
    Case->substance=1;
  } else if(strcmp(argv[1],"Ar")==0){
    printf("Ar funkar\n" );
    Case->substance=2;
  } else if(strcmp(argv[1],"Kr")==0){
    printf("Kr funkar\n" );
    Case->substance=3;
  } else if(strcmp(argv[1],"Xe")==0){
    printf("Xe funkar\n" );
    Case->substance=4;
  } else {
    fprintf(stderr, "Error in input. Not a valid substance: %s\n",argv[1] );
    exit(1);}

    //Test type input
    if(strcmp(argv[2],"omega")==0){
      printf("omega funkar\n");
      Case->test = 1;
    }else if(strcmp(argv[2],"gamma")==0){
      printf("gamma funkar\n");
      Case->test = 2;
    } else if(strcmp(argv[2],"cv")==0){
      printf("cv funkar\n");
      Case->test = 3;
    } else{
      fprintf(stderr, "Error in input. Not a valid test case: %s\n", argv[2]);
      exit(1);
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

  test_case Case;
  params pars;
  int N = 3;
  double A,B;

  check_input(argc,argv,&Case);

  set_params(Case.substance,&pars);

  A= compute_A(&pars);
  B= compute_B(&pars);
  printf("%e\n",pars.sigma );
  printf("%e\n",pars.eps );
  printf("%e\n",A );
  printf("%e\n",B );

  gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(N);

  return 0;
}
