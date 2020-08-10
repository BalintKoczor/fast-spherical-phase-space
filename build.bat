@echo off
cls
set MinGW_DIR=C:\MinGW
echo ---------------------------------------------------------------
echo MinGW needs to be installed for running this program on Windows
echo ---------------------------------------------------------------
echo.
echo ---------------------------------------------------------------
echo Assuming MinGW installation folder: %MinGW_DIR% 
echo ---------------------------------------------------------------
echo.

set PATH=%MinGW_DIR%\bin;%PATH%
cd Executables
echo "==========================================="
echo "Compiling sourcecode: Precalculate_Kernel.c"
echo "==========================================="
gcc -Ofast -march=native -std=c99 ../src/Precalculate_Kernel.c -o ./Precalculate_Kernel.exe
echo "======================"
echo "Precalculating kernels"
echo "======================"
Precalculate_Kernel.exe
echo "===================================================="
echo "Compiling sourcecode example: EfficientCalculation.c"
echo "===================================================="
gcc -Ofast -march=native -std=c99 ../src/EfficientCalculation.c -o ./EfficientCalculation.exe
echo "================"
echo "Runnning example"
echo "================"
EfficientCalculation.exe