//clone and compile the pacakge first
//https://github.com/valandil/wignerSymbols.git
#include <wignerSymbols.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <complex>
using namespace std;

double Radius(double J)
{
  return sqrt(J/(2*M_PI));
}

double gammaFactor(double J, double j)
{
  std::cout.setstate(std::ios_base::failbit); //suppress cout when using CG coefficients
  double res = Radius(J)*sqrt((4*M_PI)/(2*J+1))*WignerSymbols::clebschGordan(J, j, J, J, 0., J);
  std::cout.clear(); //reset cout
  return res;
}

double OrderZeroTensorOperatorDiagonalElement(double J, double j, double m)
{
  std::cout.setstate(std::ios_base::failbit); //suppress cout when using CG coefficients
  double res = sqrt((2*j+1)/(2*J+1))*WignerSymbols::clebschGordan(J, j, J, m, 0., m);
  std::cout.clear(); //reset cout
  return res;
}

double SpinParityOperatorElement(double m, double J, double s)
{
  double res=0.;
  for (int j=0; j <= 2*J; j++){
    res += pow( gammaFactor(J, j), -s) * sqrt((2*j+1)/(4*M_PI)) * OrderZeroTensorOperatorDiagonalElement(J, j, m);
  }
  
  return 1/Radius(J)*res;
}


void precalcParity(int Ndim)
{
	std::cout << Ndim << "\n";
	double J = ((double) (Ndim-1))/2;
	const double s = 0.;
	
	std::ostringstream filename;
	filename << "../../Calculated/Parity/ParityD";
	filename << Ndim;
	filename << ".dat";
	
	const std::string tmp = filename.str();
	const char* filetoopen = tmp.c_str();
	
	std::ofstream ofs(filetoopen, ios::out | ios::binary );
	
	//PARITY!!
	//double *paritydouble = parityOperatorFullMatrix(J, s);
	double m;
	std::complex<double> element;
	const int size = sizeof(std::complex<double>);
	for (int i=0;i<Ndim;i++){
	m =  ((double) i) -J;
	element = (std::complex<double>) SpinParityOperatorElement(m, J, s);
	ofs.write(reinterpret_cast<const char*>( &element ), size);
	}

    ofs.close();
	
}




int main(){
	
	for (int l=2; l<=50; l++){
		precalcParity(l);
	}
	return 0;
}








