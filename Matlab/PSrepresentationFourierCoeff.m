function result = PSrepresentationFourierCoeff(rho, Kcoeffs, Ndim)
fftdim = (Ndim-1)*2+1;

result = zeros(fftdim,fftdim);

for l=1:fftdim
for m=1:fftdim
    temp = complex(0,0);
    lowb = max(Ndim-m+1,1);
    upb = min(2*Ndim-m,Ndim);
    for lambda=lowb:upb
        temp = temp + rho(lambda,lambda+m-Ndim)*Kcoeffs(l, lambda, lambda+m-Ndim);
    end
    result(l,m) = temp;
end
end

end