#include <math.h>
#include <lapacke.h>
#include <complex.h>

//#include "WignerD_fftw.h"

#include <stdio.h>

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
	 
#define min(x,y) (((x) < (y)) ? (x) : (y))

#define MAXFILENAME 100

lapack_complex_double *Iy_Matrix(lapack_int Ndim) 
{
double J = (Ndim-1)/2.;
lapack_complex_double *AB = calloc(Ndim*2,sizeof(lapack_complex_double));

lapack_complex_double zero = lapack_make_complex_double(0.0,0.0);
lapack_complex_double myimag = lapack_make_complex_double(0.0,1.0);
//prepare the upper triangle of Iy
//the array is computed in FORTRAN style (transposed)
   double m;
   for (int i=0; i<Ndim; i++){
     m = i-J;
     AB[2*i] =  sqrt((J+m)*(J-m+1.)) / 2.*(-myimag);
	 AB[2*i+1] =  zero;		  	
	 //printf("Hal %f \n" , cimag(AB[i]));
    }
  
  return AB;
}


lapack_complex_double *Iy_Eigsystem_Lapacke(lapack_int Ndim) 
{

char JOBZ, UPLO;
   lapack_int KD, LDAB, LDZ, INFO;
   
   
   JOBZ='V'; //= 'V':  Compute eigenvalues and eigenvectors.
   UPLO='U'; // = 'U':  Upper triangle of A is stored;
   KD=1;  //number of superdiagonals
   LDAB=KD+1;  //The leading dimension of the array AB = 2
   LDZ=max(1,Ndim); //The leading dimension of the array eigenvectors
   
   //int rworksize=max(1,3*Ndim-2);
   
   //Initialize the array that holds the eigenvalues
    double *Eigenvalue;
	lapack_complex_double *Eigenvec, *WORK; //AB is transposed because of FORTRAN array struct.
	
	Eigenvec = calloc(Ndim*Ndim,sizeof(lapack_complex_double));
	WORK = calloc(Ndim,sizeof(lapack_complex_double));
	Eigenvalue = calloc(Ndim,sizeof(double));
	
   
    lapack_complex_double *AB = Iy_Matrix(Ndim);  //Create the upper triangle matrix
   
    INFO = LAPACKE_zhbev(LAPACK_COL_MAJOR, JOBZ, UPLO, Ndim,
                          KD, AB,
                          LDAB, Eigenvalue, Eigenvec,
                          LDZ );
	
	free(AB);
	free(Eigenvalue);
	free(WORK);
	return Eigenvec;

}








void PrecalculateEigenvec(int Ndim)
{ 
    lapack_complex_double *Eigenvec = Iy_Eigsystem_Lapacke(Ndim);
    
	char filename[MAXFILENAME+1];

	snprintf(filename, MAXFILENAME, "../../Calculated/Eigenvectors/EigenvectorsD%d.dat", Ndim);

	FILE* file = fopen(filename, "wb");
	
	if (file == NULL){
		printf("CANNOT OPEN FILE FOR WRITING: %s\n",filename);
		abort();
	}
	//else printf("WRITING TO FILE: %s\n",filename);
	
    const int size = sizeof(lapack_complex_double);
	const int fftdim = 2*(Ndim-1)+1;
    
    for (int m=0; m<Ndim; m++){
	for (int mu=0; mu<Ndim; mu++){
        if (fwrite(&Eigenvec[m + Ndim*mu], 1, size, file) != size){
			printf("ERROR WHILE WRITING TO FILE: %s\n",filename);
			abort();
		}
	}}
	free(Eigenvec);
	fclose(file);
}

int main(){
    int MAXDIM = 120;
    
	for (int Ndim=2; Ndim<=MAXDIM; Ndim++){
		PrecalculateEigenvec(Ndim);
		printf("Eigenvectors calculated for dimension %d\n", Ndim);
	}	
	return 0;
}