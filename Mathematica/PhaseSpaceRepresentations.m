(* ::Package:: *)

SphericalTensorOperator[l_,m_,J_]:=N@Sqrt[(2l+1)/(2J+1)]*Table[
If[(m1+m==m2)&&(Abs[m2]<=J),N@ClebschGordan[{J,m1},{l,m},{J,m2}],0]
,{m2,J,-J,-1},{m1,J,-J,-1}];


(* This is a simple code that decomposes a density matrix into tensor operator coeffeficients *)
tensorDecompositionCoefficients[rho_]:=Module[{J,MatrixScalarProduct},
J=N[((Dimensions@rho)[[1]]-1)/2];
MatrixScalarProduct[A_,B_]:=Tr[ConjugateTranspose[A].B];
Return@Flatten@Table[MatrixScalarProduct[SphericalTensorOperator[j,m,J],N@rho],{j,0,2J},{m,-j,j}];
]


SphericalHarmonic[j_,m_,finalpoints_]:=
Table[N@SphericalHarmonicY[j,m,th,ph],{th,0,\[Pi]-\[Pi]/finalpoints,\[Pi]/finalpoints},{ph,0,2\[Pi]-2\[Pi]/finalpoints,2\[Pi]/finalpoints}];


PSrepresentationFromTensors[rho_,finalpoints_]:=Module[{decompCoeffs,ind,dim},
dim=(Dimensions@rho)[[1]];
decompCoeffs=tensorDecompositionCoefficients[rho];
ind=0;
Return[
Sqrt[(2\[Pi])/((dim-1)/2)]*Sum[
ind+=1;
decompCoeffs[[ind]]*SphericalHarmonic[j,m,finalpoints]
,{j,0,dim-1},{m,-j,j}]
];
]


precalculatedKcoeffs[Ndim_,path2kernels_]:=Module[{filename,aa,fftdim},
fftdim =(Ndim-1)*2+1;
filename=ToString@StringForm["KernelD``.dat",Ndim];
filename=FileNameJoin[{path2kernels,filename}];
aa=Import[filename,"Complex128"];
Return@ArrayReshape[Reverse@aa,{fftdim,Ndim,Ndim}];
]


(* Compiled version of the function *)
PSrepresentationFourierCoeff=Compile[{{rho,_Complex,2},{KPRC,_Complex,3},{dim,_Integer,0}},
Table[
Sum[
rho[[\[Lambda],\[Lambda]+m-dim]]*KPRC[[l,\[Lambda],\[Lambda]+m-dim]]
,{\[Lambda],Max[dim-m+1,1],Min[2dim-m,dim]}]
,{l,1,2dim-1},{m,1,2dim-1}]];


(*PSrepresentationFourierCoeff[rho_,KPRC_,dim_]:=Table[
Total[Diagonal[rho*KPRC[[l]],m-dim]]
,{l,1,2dim-1},{m,1,2dim-1}];*)


PSrepresentationFromFourier[rho_,Kcoeffs_,finalpoints_]:=Module[{fourierCoeffs,fourierCoeffsRef,PSrep,identityPS,PSrepRef,ZeroFill,dim},
dim=(Dimensions@rho)[[1]];
fourierCoeffs=PSrepresentationFourierCoeff[rho,Kcoeffs,dim];
fourierCoeffsRef=PSrepresentationFourierCoeff[IdentityMatrix[dim],Kcoeffs,dim];
ZeroFill[in_,finalpts_]:=PadLeft[in,{2finalpts,finalpts},0.];
PSrep=Fourier[ZeroFill[fourierCoeffs,finalpoints]][[1;;finalpoints]];
identityPS=N@Sqrt[((dim-1)/2)/(2\[Pi])]*Sqrt[4\[Pi]]/Sqrt[dim](* *gammafactor *);
PSrepRef=Fourier[ZeroFill[fourierCoeffsRef,finalpoints]][[1;;finalpoints]];
Return@Chop[PSrep/PSrepRef/identityPS];
]


interpolateDiscreteRepresentation[list_,\[Theta]_,\[Phi]_]:=Module[{maxth,maxph,dat},
dat=ConstantArray[0,{Length@list+1,Length@list+1}];
dat[[1;;Length@list,1;;Length@list]]=list;
dat[[Length@list+1,;;]]=list[[Length@list,;;]];
dat[[;;,Length@list+1]]=list[[;;,Length@list]];
Return[ListInterpolation[list,{{0,\[Pi]},{0,2\[Pi]}}][\[Theta],\[Phi]]];
]
