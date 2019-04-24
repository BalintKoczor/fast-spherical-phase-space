function A = precalculatedKcoeffs(Ndim)
fftdim = (Ndim-1)*2+1;
fileID = fopen(num2str( Ndim,'../Calculated/Kernels/KernelD%d.dat') );
A = fread(fileID,[2,Ndim*fftdim*fftdim],'double');
A = complex( flip(A(1,:)), flip(A(2,:)) );
A = reshape(A,[Ndim, fftdim, fftdim]);
A = permute(A,[3,2,1]);
end