#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define TIMING 1
#define MAXFILENAME 120

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#if TIMING == 1
	#include <time.h>
    #define MINDIM 2
    #define MAXDIM 120
    #define STEPDIM 1
#else
    #define MINDIM 4
    #define MAXDIM 4
    #define STEPDIM 10
#endif

complex *ImportKernel(int Ndim){
	char filename[MAXFILENAME+1];

	snprintf(filename, MAXFILENAME, "../Calculated/Kernels/KernelD%d.dat", Ndim);

	FILE* file = fopen(filename, "rb");
	
	if (file == NULL){
		printf("CANNOT OPEN FILE: %s\n",filename);
		abort();
	}
	//else printf("READING FILE: %s\n",filename);
	
	const int size = sizeof(complex);
	const int fftdim = 2*(Ndim-1)+1;
	
	complex *kernel = calloc(fftdim*Ndim*Ndim,size);	
	
	if (fread(kernel, fftdim*Ndim*Ndim, size, file) != size){
		printf("ERROR WHILE READING FROM FILE: %s\n",filename);
		abort();
	}
	
	fclose(file);
	return(kernel);
}


void CalcPSrepresentationL(complex* rho, complex* matrL, complex* PSrepr, complex* prod, int l, int Ndim)
{

	int fftdim=2*(Ndim-1)+1;
    int dimsq = Ndim*Ndim;
    int lind, mshift, ind1, ind2, lowb, upb;
    complex result;
	
    //calculate elementwise prod
    for (int a=0; a<Ndim; a++){
    for (int b=0; b<Ndim; b++){
        prod[a + Ndim*b] = rho[a + Ndim*b] * matrL[a + Ndim*b];
    }
    }
    for (int m=0; m<fftdim; m++){
        result = 0.;
        mshift = m - Ndim +1;
        ind1 = Ndim + 1;
        ind2 = Ndim*mshift;
        lowb = MAX(Ndim - m - 1, 0);
        upb  = MIN(2*Ndim - m -1, Ndim);
        for (int lambda=lowb; lambda<upb; lambda++){
            result += prod[lambda*ind1 + ind2];
        }
        PSrepr[l + fftdim*m] = result;
    }
	
}


int main (){
#if TIMING == 1
	FILE* file = fopen("../Calculated/PhaseSpaceTimes.txt", "w");
	fclose(file);
#endif 
	
	int repetition;

	
for(int Ndim = MINDIM; Ndim <=MAXDIM; Ndim+=STEPDIM){
#if TIMING == 1
	FILE* file = fopen("../Calculated/PhaseSpaceTimes.txt", "a");
#endif    
	
    int fftdim=2*(Ndim-1)+1;
	complex *kernel = ImportKernel(Ndim);
	complex *rho = malloc(Ndim*Ndim*sizeof(complex));
    complex *PSrepr = (complex *) malloc( fftdim*fftdim* sizeof( complex ) );
    //This array will be used to temporarily store the element-wise product
    //of the matrices rho and kernelL
    complex *prod = (complex *) malloc( Ndim*Ndim* sizeof( complex ) );
	
    //Density matrix is a simple full matrix
	for(int i=0;i<Ndim*Ndim;i++){
		rho[i]=(complex) 1.;
	}
    
#if TIMING == 1
	repetition = (int) (100000000. / (Ndim*Ndim*Ndim));
	if (repetition < 1) repetition = 1;
	
	//start measuring time
	clock_t start = clock(), diff;
	//calculate phase-space representation function repeatadely for precise timing
	for (int i=0;i<repetition;i++){
		//CalcPSrepresentation(rho, kernel, PSrepr, Ndim);
        for (int l=0; l<fftdim; l++){
            CalcPSrepresentationL(rho, &kernel[Ndim*Ndim*l], PSrepr, prod, l, Ndim);
        }
	}
	//end measuring time
	diff = clock() - start;

	double msec = diff * 1000 / CLOCKS_PER_SEC;
	fprintf(file,"{%d,%f},\n", Ndim, msec/repetition);
	printf("Phase-space function calculated in %f milliseconds for dimension %d (repetitions: %d)\n", msec/repetition, Ndim,repetition);
#else
    //CalcPSrepresentation(rho, kernel, PSrepr, Ndim);
    for (int l=0; l<fftdim; l++){
        //Put code here if you do not want to import the kernel,
        //but calculate it on-the-fly via:
        //CalcMatrL(int l, complex* ptilde, complex* matr, complex* u, int Ndim)
        //This is prefferred when Ndim is large, e.g, when Ndim > 200
        //and storing the kernel is inefficient
        CalcPSrepresentationL(rho, &kernel[Ndim*Ndim*l], PSrepr, prod, l, Ndim);
    }
    printf("---------Calculated Fourier Coefficients for dimension d=%d----------\n", Ndim);
    for (int l=0; l<fftdim; l++){
        for (int m=0; m<fftdim; m++){
            printf("%f + %f\n",creal(PSrepr[l + fftdim*m]),cimag(PSrepr[l + fftdim*m]));
        }
    }
#endif
	free(kernel);
	free(rho);
    free(PSrepr);
    free(prod);
#if TIMING == 1   
    fclose(file);
#endif
}
	return(0);
}