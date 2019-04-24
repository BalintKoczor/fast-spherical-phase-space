#!/bin/sh

if ! [ -x "$(command -v gcc)" ]; then
  echo 'Error: gcc is not installed.' >&2
  exit 1
fi
if ! [ -x "$(command -v tar)" ]; then
  echo 'Error: tar is not installed.' >&2
  exit 1
fi
if [ ! -d "./Calculated" ]; then
  mkdir Calculated
  cd Calculated
  echo "===================================="
  echo "Exporting precalculated eigenvectors"
  echo "===================================="
  tar -xzvf ../Eigenvectors.tar.gz
  echo "========================================"
  echo "Exporting precalculated parity operators"
  echo "========================================"
  tar -xzvf ../Parity.tar.gz
  echo "======================================="
  echo "Exporting precalculated example Kernels"
  echo "======================================="
  tar -xzvf ../Kernels.tar.gz
  cd ..
fi
if [ ! -d "./Executables" ]; then
	mkdir Executables
fi
cd Executables
echo "==========================================="
echo "Compiling sourcecode: Precalculate_Kernel.c"
echo "==========================================="
gcc -O3 -std=c99 ../src/Precalculate_Kernel.c -o ./Precalculate_Kernel.a
echo "======================"
echo "Precalculating kernels"
echo "======================"
./Precalculate_Kernel.a
echo "===================================================="
echo "Compiling sourcecode example: EfficientCalculation.c"
echo "===================================================="
gcc -O3 -std=c99 ../src/EfficientCalculation.c -o ./EfficientCalculation.a
echo "================"
echo "Runnning example"
echo "================"
./EfficientCalculation.a