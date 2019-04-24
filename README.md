# fast_spin_phase_space

This repository contains supplementary material
for the manuscript [1]. The scripts [./build.sh](/build.sh) and [./build.bat](/build.bat)
first extract precalculated files necessary for later calculations,
then compile and run the source codes located in [./src](src).
In particular, [./src/Precalculate\_Kernel.c](src/Precalculate_Kernel.c) calculates and stores
the kernels, i.e., the coefficients K\_{\lambda}^{lm}, which will be required by the Mathematica
and Matlab scripts. Alternatively, these kernels are provided as a release for dimensions up to
120 and could be directly used for Mathematica and Matlab scripts.


## Building and testing
- on UNIX or macOS run  ```./build.sh``` -- gcc and tar need to be installed
- on Windows run ```./build.bat``` -- MinGW needs to be installed which provides both gcc
and tar for Windows
- if the compilations and calculations were successful, the precalculated kernels
should be located in [./Calculated/Kernels]().


## Content

The source code in [./src](src) contains the following two programs wirtten in C.
These are compiled and ran by [./build.sh](build.sh) and [./build.bat](build.bat).


### EfficientCalculation.c

The function

- complex\* EfficientPSrepresentation(int Ndim, complex\* rho, complex\* kernel)

in [./src/EfficientCalculation.c](src/EfficientCalculation.c)
calculates the Fourier series decomposition of the phase space function F\_{\rho}
representing the density matrix \rho. It requires the precalculated coefficients
K\_{\lambda}^{lm} as the transformation kernel.

The example provided in [./src/EfficientCalculation.c](src/EfficientCalculation.c)
calculates Fourier decomposition
coefficients of Wigner functions W\_{\rho} up to dimension 30 using the kernels
initially preculculated by [./src/Precalculate_Kernel.c](src/Precalculate_Kernel.c).

### Precalculate\_Kernel.c

The function

- complex\* TransformationKernel(complex\* parity, int Ndim) 

in [./src/Precalculate_Kernel.c](src/Precalculate_Kernel.c) calculates the
coefficients K\_{\lambda}^{lm} required for calculating the Fourier series
representation of phase-space functions. It takes precalulated parity operators as an argument.

The example provided in [./src/Precalculate_Kernel.c](src/Precalculate_Kernel.c) calculates the
kernels for Wigner functions (s=0) up to dimension 30 -- using the precalculated parity operators
provided by ./Parity.tar.gz up to dimension 120. The resulting precalculated kernels are used by
 [./src/EfficientCalculation.c](src/EfficientCalculation.c) and by the Mathematica and Matlab codes.


### Matlab example

A simple example is provided in [./Matlab/example.m](Matlab/example.m) which calculates the Wigner function
of a Schroedinger cat spinstate for a given dimension. This program requires the
precalculated kernels in [./Calculated/Kernels](Calculated/Kernels) -- these are extracted from
./Kernels.tar.gz for dimensions up to 30 when running ```./build.sh``` or ```./build.bat```.

Run the example by typing in the command line ```matlab -r "run ./Matlab/example.m" ``` or in the Matlab
terminal ```run ./Matlab/example.m```.

### Mathematica examples
