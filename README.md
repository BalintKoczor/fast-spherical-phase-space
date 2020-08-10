# fast_spin_phase_space

This repository contains implementations of the efficient computation scheme introduced in
the manuscript [1]. These can be used to compute phase-space representations of spin systems,
such as symmetric qubits states. Using a set of precalculated transformation kernels,
our approach computes Fourier decompositions of phase spaces and applies Fast Fourier Transforms
to them. We provide example code in scripting languages MATLAB, Mathematica and Python which read
in precalculated kernels and compute phase spaces via Fast Fourier Transforms.

Precalculated kernels are provided in
[Calculated\Kernels](Calculated\Kernels) up to dimension 30 (up to 29 quibts in a symmetric state).
Kernels are also provided up to dimension 120 (up to 119 qubits in a symmetric state) in the release.

For larger dimensions the scripts [build.sh](build.sh) and [build.bat](build.bat)
can be used to compile and run the source code located in [src](src), as discussed below.
The C code [src/Precalculate\_Kernel.c](src/Precalculate_Kernel.c) is used to calculate and store
the kernels which will be used by the MATLAB, Mathematica and Python scripts.


## Example programs


### Matlab example

A simple example is provided in [Matlab/example.m](Matlab/example.m) which calculates the Wigner function
of a Schrödinger cat state, the well-knwon GHZ state, for a given dimension. This program uses the
precalculated kernels located in [Calculated/Kernels](Calculated/Kernels). Note that if the desired
dimension is larger than 30, then the corresponding kernel needs to be downloaded from the release and
extracted to [Calculated\Kernels](Calculated\Kernels). Alternatively, run the script [build.sh](build.sh)
or [build.bat](build.bat) to compute the kernels with our C code. 

Run the example by entering the directory [Matlab/](Matlab/) and
type in the command line ```matlab -r "run ./Matlab/example.m" ``` or in the Matlab
terminal ```run ./Matlab/example.m```. This will open a plot of the GHZ state as shown
below.





### Mathematica example

A simple example as a Mathematica Notebook is provided in [Mathematica/example.nb](Mathematica/example.nb)
which calculates the Wigner function of a Schrödinger cat state, the well-knwon GHZ state, for a given dimension.
This program uses the precalculated kernels located in [Calculated/Kernels](Calculated/Kernels). Note that if the desired
dimension is larger than 30, then the corresponding kernel needs to be downloaded from the release and
extracted to [Calculated\Kernels](Calculated\Kernels). Alternatively, run the script [build.sh](build.sh)
or [build.bat](build.bat) to compute the kernels with our C code. 

Run the example by openeing the file [Mathematica/example.nb](Mathematica/example.nb)
in Mathematica and by selecting "Evaluate Notebook" in the "Evaluation" drop-down menu. 



### Python example


## The source code

The source code in [src](src) contains the following two programs wirtten in C.
These are compiled and ran by [build.sh](build.sh) and [build.bat](build.bat).

### Building and testing
- on UNIX or macOS run  ```./build.sh``` -- gcc and tar need to be installed
- on Windows run ```./build.bat``` -- MinGW needs to be installed which provides both gcc
and tar for Windows
- if the compilations and calculations were successful, the precalculated kernels
should be located in [Calculated/Kernels]().


### EfficientCalculation.c

The function

- complex\* EfficientPSrepresentation(int Ndim, complex\* rho, complex\* kernel)

in [src/EfficientCalculation.c](src/EfficientCalculation.c)
calculates the Fourier series decomposition of the phase space function F\_{\rho}
representing the density matrix \rho. It requires the precalculated coefficients
K\_{\lambda}^{lm} as the transformation kernel.

The example provided in [src/EfficientCalculation.c](src/EfficientCalculation.c)
calculates Fourier decomposition
coefficients of Wigner functions W\_{\rho} up to dimension 30 using the kernels
initially preculculated by [src/Precalculate_Kernel.c](src/Precalculate_Kernel.c).

### Precalculate\_Kernel.c

The function

- complex\* TransformationKernel(complex\* parity, int Ndim) 

in [src/Precalculate_Kernel.c](src/Precalculate_Kernel.c) calculates the
coefficients K\_{\lambda}^{lm} required for calculating the Fourier series
representation of phase-space functions. It takes precalulated parity operators as an argument.

The example provided in [src/Precalculate_Kernel.c](src/Precalculate_Kernel.c) calculates the
kernels for Wigner functions (s=0) up to dimension 30 -- using the precalculated parity operators
provided by [Parity.tar.gz](Parity.tar.gz) up to dimension 120. The resulting precalculated kernels are used by
 [src/EfficientCalculation.c](src/EfficientCalculation.c) and by the Mathematica and Matlab codes.



