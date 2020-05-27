function result = PSrepresentationFourierCoeff(rho, Kcoeffs, Ndim)
fftdim = (Ndim-1)*2+1;

result = zeros(fftdim,fftdim);

for l=1:fftdim
for m=1:fftdim
    temp = reshape( Kcoeffs(l,:,:), Ndim, Ndim);
    result(l,m) = sum( diag(rho .* temp, m-Ndim) );
end
end

end