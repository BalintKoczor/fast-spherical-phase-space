function result = PSrepresentationFromFourier(rho, Kcoeffs, Ndim, finalpoints)
psrep = PSrepresentationFourierCoeff(rho, Kcoeffs, Ndim);
psrepREF = PSrepresentationFourierCoeff(eye(Ndim)/sqrt(Ndim), Kcoeffs, Ndim);
result = fft2(psrep,2*finalpoints, finalpoints);
resultREF = fft2(psrepREF,2*finalpoints, finalpoints);
J=(Ndim-1)/2;
identityPS = 1/sqrt(2*J);
result = result ./ resultREF*identityPS;
result = result(1:finalpoints,:);
end
