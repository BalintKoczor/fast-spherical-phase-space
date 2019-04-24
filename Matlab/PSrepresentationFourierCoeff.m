function result = PSrepresentationFourierCoeff(rho, Kcoeffs, Ndim)
fftdim = (Ndim-1)*2+1;

result = zeros(fftdim,fftdim);

for l=1:fftdim
for m=1:fftdim
    temp = complex(0,0);
    for lambda=1:Ndim
        if (Ndim+1 <= lambda+m)&&(lambda+m <= fftdim+1)
            temp = temp + rho(lambda,lambda+m-Ndim)*Kcoeffs(l, m, lambda);
        end
    end
    result(l,m) = temp;
end
end

end