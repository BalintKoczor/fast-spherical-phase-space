% import and show the precalculated K matrices from
%../Calculated/Kernel -- these correspond to the Wigner function
Ndim = 2; %set the dimension Ndim=2J+1 of the system
Kcoeffs = precalculatedKcoeffs(Ndim);
disp(['Example of precalculatedtransofmation matrices for dimession ', num2str(Ndim)])
for l=1:(2*Ndim-1)
    disp(['l = ', num2str(l-Ndim)])
    disp(squeeze(Kcoeffs(l,:,:)))
end

%set the dimension Ndim=2J+1 of the system
Ndim = 10;
%import the precalculated K coefficients
%from ../Calculated/Kernel which correspond to
%the Wigner function
Kcoeffs = precalculatedKcoeffs(Ndim);

%set the number of points along
%THETA and PHI after the Fourier transformation
%increasing the number of points results in a
%smoother plot
finalpoints = 128;

%prepare the deisred density matrix
%this is a Schroedinger cat state now
vec = zeros(Ndim);
vec(1) = 1/sqrt(2);
vec(Ndim) = 1/sqrt(2);
rho = mtimes(vec,ctranspose(vec));

%calculate the Wigner function of the density matrix
psrep = PSrepresentationFromFourier(rho, Kcoeffs, Ndim, finalpoints);

%plot the Wigner function
surf(real(psrep))