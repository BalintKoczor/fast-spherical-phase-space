#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define TIMING 0
#define MAXFILENAME 100

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#if TIMING == 1
	#include <time.h>
    #define MINDIM 10
    #define MAXDIM 500
    #define STEPDIM 10
#else
    #define MINDIM 2
    #define MAXDIM 120
    #define STEPDIM 1
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
	complex* Parity = calloc(Ndim,sizeof(complex));
	//Thi is implemented taking into account that parity operators are diagonal
	if (fread(Parity, Ndim, size, file) != size){
			printf("ERROR WHILE READING FILE: %s\n",filename);
			abort();
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
	complex* Eigenvec = malloc(Ndim*Ndim*size);
    if (fread(Eigenvec, Ndim*Ndim, size, file) != size){
		printf("ERROR WHILE READING FILE: %s\n",filename);
		abort();
	}
	
	fclose(file);
	return(Eigenvec);
}



complex *ParityTilde(complex* parity, complex* u, int Ndim) 
{
    complex *P = calloc(Ndim*Ndim,sizeof(complex));
    
    for (int a=0; a<Ndim; a++){
        for (int b=0; b<Ndim; b++){
            for (int xi=0; xi<Ndim; xi++){
                P[a + b*Ndim] += parity[xi] * u[a + xi*Ndim]*conj(u[b + xi*Ndim]);
            }
        }
    }
    return(P);
}

void CalcMatrL(int l, complex* ptilde, complex* matr, complex* u, int Ndim) 
{
    int lshift = l - Ndim +1; 
    int ind1, ind2, ind3, ind4, ind5, lowb, upb;
    complex element;
    
    for (int a=0; a<Ndim; a++){
        for (int b=0; b<Ndim; b++){
            element = 0.;
            // storing the indexes
            ind1 = Ndim+1;
            ind2 = lshift *Ndim;
            ind3 = Ndim*a;
            ind4 = lshift + Ndim*b;
            lowb = MAX(Ndim - l - 1, 0);
            upb  = MIN(2*Ndim - l -1, Ndim);
            for (int nu=lowb; nu<upb; nu++){
                //if ( (Ndim > nu + lshift)&&(nu + lshift >= 0) ){
                element += ptilde[nu*ind1 + ind2]\
                *u[nu + ind3]*conj(u[nu + ind4]);
                //}
            }
            matr[a + b*Ndim] = element;
        }
    }
}


void PrecalculateKernel(complex *parity, complex *u, complex *ptilde, int Ndim){
	 
    const int fftdim = 2*(Ndim-1)+1;
    complex *matr = calloc(Ndim*Ndim,sizeof(complex));
    
     //for testing -- prints out a matrix K_l
    /*CalcMatrL(0, ptilde, matr, u, Ndim);
    for (int a=0; a<Ndim; a++){
        for (int b=0; b<Ndim; b++){
            printf("%f + %f\n",creal(matr[a + Ndim * b]),cimag(matr[a + Ndim * b]));
        }
    }*/
#if TIMING == 0	   
	char filename[MAXFILENAME+1];
	snprintf(filename, MAXFILENAME, "../Calculated/Kernels/KernelD%d.dat", Ndim);
	FILE* file = fopen(filename, "wb");

	if (file == NULL){
		printf("CANNOT OPEN FILE FOR WRITING: %s\n",filename);
		abort();
	}
	//else printf("WRITING TO FILE: %s\n",filename);
	
    const int size = sizeof(complex);

	for (int l=0; l<fftdim; l++){
        CalcMatrL(l, ptilde, matr, u, Ndim);
        if (fwrite(matr, Ndim*Ndim, size, file) != size){
            printf("ERROR WHILE WRITING TO FILE: %s\n",filename);
            abort();
        }
    }
	fclose(file);
#else
	for (int l=0; l<fftdim; l++){
        CalcMatrL(l, ptilde, matr, u, Ndim);
        //Instea of storing the kernel
        //the PS function can be calculated here
    }
#endif
    free(matr);
}



#if TIMING == 1
int main(){
	FILE* file = fopen("../Calculated/KernelTimes.txt", "w");
	fclose(file);
	
	int repetition;

	for (int Ndim=MINDIM; Ndim<=MAXDIM; Ndim+=STEPDIM){
		FILE* file = fopen("../Calculated/KernelTimes.txt", "a");
		
        repetition = (int) (100000000. / (Ndim*Ndim*Ndim*Ndim));
        if (repetition < 1) repetition = 1;
	
        complex *parity = ImportParity(Ndim);
        complex *u      = ImportWignerDEigenvalues(Ndim);
        complex *ptilde = ParityTilde(parity, u, Ndim) ;
        free(parity);
        
		//start measuring time
		clock_t start = clock(), diff;
		
		//Calculate and save kernel
		for (int i=0;i<repetition;i++){
			PrecalculateKernel(parity, u, ptilde, Ndim);
		}
		
		//end measuring time
		diff = clock() - start;
        
        free(u);
        free(ptilde);

		double msec = diff * 1000 / CLOCKS_PER_SEC;
		fprintf(file, "{%d,%f},\n", Ndim, msec/repetition);
		printf("Kernel calculated in %f milliseconds for dimension %d (repetitions: %d)\n", msec/repetition, Ndim,repetition);
		fclose(file);
	}
	
	return 0;
}
#else
int main(){
	for (int Ndim=MINDIM; Ndim<=MAXDIM; Ndim+=STEPDIM){
        complex *parity = ImportParity(Ndim);
        complex *u      = ImportWignerDEigenvalues(Ndim);
        complex *ptilde = ParityTilde(parity, u, Ndim) ;
        free(parity);
		PrecalculateKernel(parity, u, ptilde, Ndim);
        free(u);
        free(ptilde);
		printf("Kernel calculated and stored for dimension %d\n", Ndim);
	}	
	return 0;
}
#endif