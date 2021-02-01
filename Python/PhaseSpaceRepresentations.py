import numpy as np
import matplotlib.pyplot as plt

def precalculatedKcoeffs(Ndim, path2kernels):
    fftdim = 2*Ndim-1;
    kernel = np.fromfile(path2kernels+'KernelD'+str(Ndim)+'.dat', dtype=np.complex128)
    kernel = (np.flip(kernel)).reshape((fftdim,Ndim,Ndim))
    return kernel
    
def PSrepresentationFourierCoeff(rho, Kcoeffs):
    Ndim = len(rho)
    fftdim = 2*Ndim-1;
    result = np.zeros((fftdim,fftdim),dtype=np.complex128)
    #This is an implementation using the offset trace operation
    #of the element-wise product of matrices
    for l in range(0, fftdim):
        prod = np.multiply(rho, Kcoeffs[l])
        for m in range(0, fftdim):
            result[l,m] = np.trace(prod,m-Ndim+1)
    return result


def PSrepresentationFourierCoeff_LOOP(rho, Kcoeffs):
    Ndim = len(rho)
    fftdim = 2*Ndim-1;
    result = np.zeros((fftdim,fftdim),dtype=np.complex128)
    #This is straightforwad loop-based implementation which might be slow
    for l in range(0, fftdim):
        for m in range(0, fftdim):
            pres = 0
            lowb = max(Ndim - m - 1, 0)
            upb  = min(2*Ndim - m -1, Ndim)
            for lam in range(lowb, upb):
                pres = pres + rho[lam,lam+m-Ndim+1]*Kcoeffs[l,lam,lam+m-Ndim+1]
            result[l,m] = pres
    return result
    
def PSrepresentationFromFourier(rho, Kcoeffs, finalpoints):
    Ndim = len(rho)
    fourierCoeffs = PSrepresentationFourierCoeff(rho, Kcoeffs)
    fourierCoeffsRef = PSrepresentationFourierCoeff(np.identity(Ndim)/np.sqrt(Ndim), Kcoeffs)

    PSrep = np.fft.ifft2(fourierCoeffs,(2*finalpoints,finalpoints))
    PSrepRef = np.fft.ifft2(fourierCoeffsRef,(2*finalpoints,finalpoints))
    Pi=np.pi
    identityPS= 1/np.sqrt(Ndim-1)
    result = np.real((PSrep/PSrepRef*identityPS).round(12))
    result = result[0:finalpoints]
    return result
    
def PSrepPlot(psFunctionArray, filename):
    ny, nx = psFunctionArray.shape

    x = np.linspace(0, 2*np.pi, nx)
    y = np.linspace(0, np.pi, ny)

    xv, yv = np.meshgrid(x, y)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(60, 30)
    ax.plot_surface(xv, yv, psFunctionArray, cmap='inferno', rstride=1, cstride=1, alpha=None, antialiased=True)
    plt.savefig(filename)