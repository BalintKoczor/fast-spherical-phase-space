import PhaseSpaceRepresentations as ps
import numpy as np

# import and show the precalculated K matrices from
#../Calculated/Kernel -- these correspond to the Wigner function
Ndim = 2 #set the dimension Ndim=2J+1 of the system
path2kernels = '../Calculated/Kernels/'
Kcoeffs = ps.precalculatedKcoeffs(Ndim, path2kernels)
print('Example of precalculated transofmation matrices\nfor dimession', Ndim)
for l in range(0,2*Ndim-1):
    print('l = ', l-Ndim+1,'\n', Kcoeffs[l],'\n')
    
# set the dimension dim=2J+1 of the system *)
Ndim = 10

# import the precalculated K coefficients from
#../Calculated/Kernel -- these correspond to the Wigner function
Kcoeffs = ps.precalculatedKcoeffs(Ndim, path2kernels)

#Set the number of points along THETA and PHI after the Fourier 
#transformation -- increasing the number of points results in a smoother plot 
finalpoints=256

#Prepare the deisred density matrix -- this is a Schroedinger cat state now 
psi = np.zeros(Ndim)
psi[0] = 1/np.sqrt(2)
psi[Ndim-1] = 1/np.sqrt(2)
rho = np.outer(psi,np.conjugate(psi))

#Calculate the Wigner function of the density matrix
result = ps.PSrepresentationFromFourier(rho, Kcoeffs, finalpoints)
#plot the Wigner function and save to file
filename = 'plot_D'+str(Ndim)+'.png'
ps.PSrepPlot(result, filename)
print('phase-space plot saved to file: ' + filename)
