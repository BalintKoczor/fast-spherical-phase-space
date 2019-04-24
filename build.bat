@echo off
set MinGW_DIR=C:\MinGW
echo ---------------------------------------------------------------
echo MinGW needs to be installed for running this program on Windows
echo Assuming MinGW installation folder: %MinGW_DIR%
echo ---------------------------------------------------------------
set PATH=%MinGW_DIR%\msys\1.0\bin;%MinGW_DIR%\bin;%PATH%
%MinGW_DIR%\msys\1.0\bin\bash.exe .\build.sh