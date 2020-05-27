//clone and compile the pacakge first
//https://github.com/valandil/wignerSymbols.git
#include "wignerSymbols.h"
#include <iostream>
#include <chrono>
using namespace std;

double *TensorOperator(double J, double j, double m)
{
  int Ndim = (int) 2*J+1;
  double *TensorOperator, m1, m2;
  double prefactor = sqrt((2*j+1)/(2*J+1));
  
  TensorOperator =(double*) calloc(Ndim*Ndim, sizeof(double));
  
  
  for (int k=0; k<Ndim; k++){
    for (int l=0; l<Ndim; l++){
	  m1 = (double)k-J; m2 = (double)l-J;
	  if( ( m1+m==m2 )&&( abs(m2)<=J ) )
	    TensorOperator[k*Ndim + l] = prefactor*WignerSymbols::clebschGordan(J, j, J, m1, m, m2);
      else TensorOperator[k*Ndim + l] = 0;
	}
  }
  
  
  return TensorOperator;

}


int main(){
int Ndim = 50; //(dim = 2J+1)

for (Ndim=2; Ndim<200; Ndim++){

auto start = std::chrono::high_resolution_clock::now();
double J = ( (double) Ndim -1 )/2;
double j,m;
double *tt;
//std::cout.setstate(std::ios_base::failbit); //suppress cout when using CG coefficients
for (int k=0; k<Ndim; k++){
    for (int l=0; l<2*k+1; l++){
        j = (double) k ;
        m = (double) k - (double) l;
        //cout << j << "," << m << "\n";
        tt = (double *) TensorOperator(J, j, m);
        free(tt);
    }
}
auto finish = std::chrono::high_resolution_clock::now();
//std::cout.clear(); //reset cout

std::chrono::duration<double> elapsed = finish - start;
std::cout << "{" << Ndim << "," << elapsed.count() << "},\n";
}


return(0);
}