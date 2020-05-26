#!/bin/sh
if ! [ -x "$(command -v gcc)" ]; then
  echo 'Error: gcc is not installed.' >&2
  exit 1
fi
if [ ! -d "./Executables" ]; then
	mkdir Executables
fi
cd Executables
echo "==========================================="
echo "Compiling sourcecode: Precalculate_Kernel.c"
echo "==========================================="
gcc -Ofast -march=native -std=c99 ../src/Precalculate_Kernel.c -o ./Precalculate_Kernel.a
echo "======================"
echo "Precalculating kernels"
echo "======================"
./Precalculate_Kernel.a
echo "===================================================="
echo "Compiling sourcecode example: EfficientCalculation.c"
echo "===================================================="
gcc -Ofast -march=native -std=c99 ../src/EfficientCalculation.c -o ./EfficientCalculation.a
echo "================"
echo "Runnning example"
echo "================"
./EfficientCalculation.a