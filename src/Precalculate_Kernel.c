#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define MAXDIM 30
#define TIMING 0
#define MAXFILENAME 100

#if TIMING == 1
	#include <time.h>
#endif


complex *ImportParity(int Ndim){
	char filename[MAXFILENAME+1];

	snprintf(filename, MAXFILENAME, "../Calculated/Parity/ParityD%d.dat", Ndim);

	FILE* file = fopen(filename, "rb");
	
	if (file == NULL){
		printf("CANNOT OPEN FILE: %s\n",filename);
		abort();
	}
	//else printf("FILE OPENED: %s\n",filename);
	
	const int size = sizeof(complex);
	complex* Parity = calloc(Ndim*Ndim,sizeof(complex));
	
	for (int m=0; m<Ndim; m++){
		if (fread(&Parity[m + Ndim*m], 1, size, file) != size){
			printf("ERROR WHILE READING FILE: %s\n",filename);
			abort();
		}
	}
	
	fclose(file);
	return(Parity);
}


complex *ImportWignerDEigenvalues(int Ndim){
	char filename[MAXFILENAME+1];

	snprintf(filename, MAXFILENAME, "../Calculated/Eigenvectors/EigenvectorsD%d.dat", Ndim);

	FILE* file = fopen(filename, "rb");
	
	if (file == NULL){
		printf("CANNOT OPEN FILE: %s\n",filename);
		abort();
	}
	//else printf("FILE OPENED: %s\n",filename);
	
	const int size = sizeof(complex);
	complex* Eigenvec = malloc(Ndim*Ndim*sizeof(complex));
	
	for (int m=0; m<Ndim; m++){
	for (int mu=0; mu<Ndim; mu++){
		if (fread(&Eigenvec[m + Ndim*mu], 1, size, file) != size){
			printf("ERROR WHILE READING FILE: %s\n",filename);
			abort();
		}
	}}
	
	fclose(file);
	return(Eigenvec);
}



complex *WignerD_Fourier_Coefficients(int Ndim) 
{
	//Import eigenvectors
	complex *Eigenvec = ImportWignerDEigenvalues(Ndim); 
	//initialize array for the Fourier coefficients
	complex *Fcoeff = calloc(Ndim*Ndim*Ndim,sizeof(complex));
    
	const int dimsq = Ndim*Ndim;
	
	for (int mu=0; mu<Ndim; mu++){
	  for (int m=0; m<Ndim; m++){
	    for (int n=0; n<Ndim; n++){
                Fcoeff[n + Ndim*m + dimsq*mu]=Eigenvec[m + Ndim*mu]*conj(Eigenvec[n + Ndim*mu]);
	} 
	  } 
	    }
	
	free(Eigenvec);
    return Fcoeff;
}



complex KernelElement(complex* parity, complex *Coeffs, int l, int m, int lambda, int Ndim) 
{
	//This calculates the double sum that corresponds to a single element of the kernel
	//the kernel indexes are given with standard C indexing: 
	// 0 <= l,m < 4J+1 and 0 <= lambda, nu < 2J+1
	int dimsq = Ndim*Ndim;
	//double J = ((double) (Ndim-1))/2;
	int fftdim=2*(Ndim-1)+1;
	complex result = 0.;
	complex b,c;

	for (int nu=0; nu<Ndim; nu++){
		for (int xi=0; xi<Ndim; xi++){
			if ( (Ndim-1 <= nu+l)&&(nu+l < fftdim))
			// the sum only runs for indexes that satisfy (in quantum-mechanical indexing)
			// -J <= nu <= +J and -J <= nu + l <= +J AND
			// -J <= lambda <= +J and -J <= lambda + m < +J
			// now shifting lambda and nu by adding J to their values and l,m by adding 2J to their values
			// We obtain the C style indexing:
			//  0 <= nu < 2J+1 and 2J <= nu + l < 4J+1 AND
			//  0 <= lambda < 2J+1 and 2J <= lambda + m < 4J+1
			// this condition is checked in the if!!!
			// and note that when indexing: the indexes (nu+l) and (lambda+m) are shifted by (Ndim-1)
			// such that they point to array index values between 0 and 2J+1
			{
				b = Coeffs[(lambda+m-(Ndim-1)) + Ndim*(xi) + dimsq*(nu+l-(Ndim-1))];
				c = conj(Coeffs[(lambda) + Ndim*(xi) + dimsq*(nu)]);
				result = result + parity[xi + Ndim*xi] * b * c;
			}
		}
	}

	return result;

}



complex *TransformationKernel(complex* parity, int Ndim) 
{

	const int fftdim=2*(Ndim-1)+1;
	const int fftsqdim = fftdim*fftdim;
	//double J = ((double) (Ndim-1))/2;
	
	//calculate the WignerD Fourier coefficients
	complex *FourierCoeffs = WignerD_Fourier_Coefficients(Ndim);
	
	complex *kernel = calloc(fftdim*fftdim*Ndim,sizeof(complex));
	
	//Calculate the kernel elementwise
	for (int lambda=0; lambda<Ndim; lambda++){
	  for (int m=0; m<fftdim; m++){
	    for (int l=0; l<fftdim; l++){
			if ( (Ndim-1 <= lambda+m)&&(lambda+m < fftdim) )
			//refer to KernelElement for the explanation of this index condition
			{
				kernel[l + fftdim*m + fftsqdim*lambda] = KernelElement(parity, FourierCoeffs, l, m, lambda, Ndim) ;
			}
	} 
	  } 
	    }
	
	free(FourierCoeffs);
    return kernel;
}





void PrecalculateKernel(int Ndim){
	
	complex *parity = ImportParity(Ndim);
	complex *kernel = TransformationKernel(parity, Ndim);
	
	
	#define MAXFILENAME 100
	char filename[MAXFILENAME+1];

	snprintf(filename, MAXFILENAME, "../Calculated/Kernels/KernelD%d.dat", Ndim);

	FILE* file = fopen(filename, "wb");
	
	if (file == NULL){
		printf("CANNOT OPEN FILE FOR WRITING: %s\n",filename);
		abort();
	}
	//else printf("WRITING TO FILE: %s\n",filename);
	
	const int size = sizeof(complex);
	const int fftdim = 2*(Ndim-1)+1;
	
	for (int l=0; l<fftdim; l++){
	for (int m=0; m<fftdim; m++){
	for (int lambda=0; lambda<Ndim; lambda++){
		if (fwrite(&kernel[l + fftdim*m + fftdim*fftdim*lambda], 1, size, file) != size){
			printf("ERROR WHILE WRITING TO FILE: %s\n",filename);
			abort();
		}
	}
	}}
	
	fclose(file);
	free(parity);
	free(kernel);
}



#if TIMING == 1
int main(){
	FILE* file = fopen("../Calculated/KernelTimes.txt", "w");
	fclose(file);
	
	int repetition;

	for (int Ndim=2; Ndim<=MAXDIM; Ndim++){
		FILE* file = fopen("../Calculated/KernelTimes.txt", "a");
		
		if (Ndim < 10) repetition = 500;
		else if (Ndim < 20) repetition = 50;
		else if (Ndim < 30) repetition = 20;
		else if (Ndim < 50) repetition = 5;
		else repetition = 1;
	
		//start measuring time
		clock_t start = clock(), diff;
		
		//Calculate and save kernel
		for (int i=0;i<repetition;i++){
			PrecalculateKernel(Ndim);
		}
		
		//end measuring time
		diff = clock() - start;

		double msec = diff * 1000 / CLOCKS_PER_SEC;
		fprintf(file, "{%d,%f},\n", Ndim, msec/repetition);
		printf("Kernel calculated in %f milliseconds for dimension %d\n", msec/repetition, Ndim);
		printf("(overall %f ms, repetitions: %d)\n", msec, repetition);
		fclose(file);
	}
	
	return 0;
}
#else
int main(){
	for (int Ndim=2; Ndim<=MAXDIM; Ndim++){
		PrecalculateKernel(Ndim);
		printf("Kernel calculated for dimension %d\n", Ndim);
	}	
	return 0;
}
#endif






