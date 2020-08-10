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
type in the command line ```matlab -r "run ./example.m" ``` or in the Matlab
terminal ```run ./example.m```. This will open a plot of the GHZ state as shown
below.

![GHZ state computed in Matlab](/Matlab/example_plot_D10.png)



### Mathematica example

A simple example as a Mathematica Notebook is provided in [Mathematica/example.nb](Mathematica/example.nb)
which calculates the Wigner function of a Schrödinger cat state, the well-knwon GHZ state, for a given dimension.
This program uses the precalculated kernels located in [Calculated/Kernels](Calculated/Kernels). Note that if the desired
dimension is larger than 30, then the corresponding kernel needs to be downloaded from the release and
extracted to [Calculated\Kernels](Calculated\Kernels). Alternatively, run the script [build.sh](build.sh)
or [build.bat](build.bat) to compute the kernels with our C code. 

Run the example by openeing the file [Mathematica/example.nb](Mathematica/example.nb)
in Mathematica and by selecting "Evaluate Notebook" in the "Evaluation" drop-down menu.
This will generate a plot of the GHZ state as shown below.

![GHZ state computed in Mathematica](/Mathematica/example_plot_D10.png)



### Python example

A simple example is provided in [Python/example.py](Python/example.py) which calculates the Wigner function
of a Schrödinger cat state, the well-knwon GHZ state, for a given dimension. This program uses the
precalculated kernels located in [Calculated/Kernels](Calculated/Kernels). Note that if the desired
dimension is larger than 30, then the corresponding kernel needs to be downloaded from the release and
extracted to [Calculated\Kernels](Calculated\Kernels). Alternatively, run the script [build.sh](build.sh)
or [build.bat](build.bat) to compute the kernels with our C code. 

Run the example by entering the directory [Python/](Python/) and
type in the command line ```python ./example.py```.
This will create a plot of the GHZ state and save it in the file
"plot_D10.png" as shown below.


![GHZ state computed in Python](/Python/example_plot_D10.png)


## The C source code

The source code in [src](src) contains the following two programs wirtten in C.
These can be compiled and ran by [build.sh](build.sh) and [build.bat](build.bat).
- on UNIX or macOS run  ```./build.sh``` -- gcc needs to be installed
- on Windows run ```./build.bat``` -- MinGW needs to be installed which provides gcc
for Windows
- if the compilations and calculations were successful, the precalculated kernels
should be located in [Calculated/Kernels]().



### Precalculate\_Kernel.c

The function

- void PrecalculateKernel(complex\* parity, complex\* u, complex\* ptilde, int Ndim)

in [src/Precalculate_Kernel.c](src/Precalculate_Kernel.c) calculates and stores the
kernel coefficients K\_{\lambda}^{lm} required for calculating the Fourier series
representation of phase-space functions. It requires both precalulated parity operators 
and eigenvectors of the rotation operator,
which are provided up to dimension 120 in [Calculated/Parity](Calculated/Parity) and in
[Calculated/Eigenvectors](Calculated/Eigenvectors).

The example provided in [src/Precalculate_Kernel.c](src/Precalculate_Kernel.c) calculates the
kernels for Wigner functions (s=0) up to dimension 120 and stores them in [Calculated/Kernels](Calculated/Kernels).
These resulting precalculated kernels are used by
[src/EfficientCalculation.c](src/EfficientCalculation.c) and by the Matlab, Mathematica and Python codes.
 
 
### EfficientCalculation.c

The function

- void CalcPSrepresentationL(complex\* rho, complex\* matrL, complex\* PSrepr, complex\* prod, int l, int Ndim)

in [src/EfficientCalculation.c](src/EfficientCalculation.c)
calculates the Fourier series decomposition of the phase-space function F\_{\rho}
representing the density matrix \rho for a fixed index l. It requires the precalculated transformation kernel
as the coefficients K\_{\lambda}^{lm}.

The example provided in [src/EfficientCalculation.c](src/EfficientCalculation.c)
calculates Fourier decomposition
coefficients of Wigner functions W\_{\rho} up to dimension 120 using the kernels
initially preculculated by [src/Precalculate_Kernel.c](src/Precalculate_Kernel.c).


## Further source code

We provide further source code for precalculating the following
for arbitrary large dimensions. Note that these depend on external
libraries, such as LAPACK, etc.

- computing eigenvectors of the rotation operator via [src/Precalculate_Parity/Precalc_Parity.cpp](src/Precalculate_Parity/Precalc_Parity.cpp)
- computing tensor operators via [src/Precalculate_Tensor_Operators/Tensor_Operator.cpp](src/Precalculate_Tensor_Operators/Tensor_Operator.cpp)
- computing parity operators via [src/Precalculate_Parity/Precalc_Parity.cpp](src/Precalculate_Parity/Precalc_Parity.cpp)



