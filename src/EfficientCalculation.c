#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>

#define MAXDIM 30
#define MAXFILENAME 100

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
	
	complex *kernel = calloc(fftdim*fftdim*Ndim,sizeof(complex));	
	
	for (int l=0; l<fftdim; l++){
	for (int m=0; m<fftdim; m++){
	for (int lambda=0; lambda<Ndim; lambda++){
		if (fread(&kernel[l + fftdim*m + fftdim*fftdim*lambda], 1, size, file) != size){
			printf("ERROR WHILE READING FROM FILE: %s\n",filename);
			abort();
		}
	}
	}}
	
	fclose(file);
	return(kernel);
}



complex* EfficientPSrepresentation(int Ndim, complex* rho, complex* kernel)
{

	int fftdim=2*(Ndim-1)+1;
	//double J = ((double) (Ndim-1))/2;
	
	
	complex *PSrepr = (complex *) malloc( fftdim*fftdim* sizeof( complex ) );
	
	for (int l=0; l<fftdim; l++){
	for (int m=0; m<fftdim; m++){
	for (int lambda=0; lambda<Ndim; lambda++){
		if ( (Ndim-1 <= lambda+m)&&(lambda+m < fftdim) )
		//refer to KernelElement for the explanation of this index condition
		{
			PSrepr[l + fftdim*m] += rho[lambda + Ndim*(lambda+m-(Ndim-1))] * kernel[l + fftdim*m + fftdim*fftdim*lambda];
		}
	}
	}
	}
	
	return PSrepr;
}



int main (){
	FILE* file = fopen("../Calculated/PhaseSpaceTimes.txt", "w");
	fclose(file);
	
	int repetition;

	
for(int Ndim = 2; Ndim <=MAXDIM; Ndim++){
	FILE* file = fopen("../Calculated/PhaseSpaceTimes.txt", "a");
	
	complex *kernel = ImportKernel(Ndim);
	complex *rho = malloc(Ndim*Ndim*sizeof(complex));
	complex* PSrep;
	//Density matrix is a simple full matrix
	for(int i=0;i<Ndim*Ndim;i++){
		rho[i]=(complex) 1.;
	}
		
	repetition = (int) (100000000. / (Ndim*Ndim*Ndim));
	if (repetition < 1) repetition = 1;
	
	//start measuring time
	clock_t start = clock(), diff;
	//calculate phase-space representation function repeatadely
	for (int i=0;i<repetition;i++){
		PSrep = EfficientPSrepresentation(Ndim, rho, kernel);
		free(PSrep);
	}
	//end measuring time
	diff = clock() - start;

	double msec = diff * 1000 / CLOCKS_PER_SEC;
	fprintf(file,"{%d,%f},\n", Ndim, msec/repetition);
	printf("Phase-space function calculated in %f milliseconds for dimension %d\n", msec/repetition, Ndim);
	printf("(overall %f ms, repetitions: %d)\n", msec, repetition);
	free(kernel);
	free(rho);
	fclose(file);
}
	return(0);
}




