function result = PSrepresentationFromFourier(rho, Kcoeffs, Ndim, finalpoints)
psrep = PSrepresentationFourierCoeff(rho, Kcoeffs, Ndim);
psrepREF = PSrepresentationFourierCoeff(eye(Ndim), Kcoeffs, Ndim);
result = fft2(psrep,2*finalpoints, finalpoints);
resultREF = fft2(psrepREF,2*finalpoints, finalpoints);
identityPS = sqrt(((Ndim-1)/2)/(2*pi))*sqrt(4*pi)/sqrt(Ndim);
result = result ./ resultREF/identityPS;
result = result(1:finalpoints,:);
end
