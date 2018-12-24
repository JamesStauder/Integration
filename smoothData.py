


import h5py, sys
import numpy as np
import scipy as sp
import scipy.ndimage
import time
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian2DKernel, Tophat2DKernel
def main(argv):
    start = time.time()
    dataFile = 'data/GreenlandInBedCoord_V2.h5'
    f = h5py.File(dataFile, 'r+')
    surface = np.asarray(f['surface'][:])
    thickness = np.asarray(f['thickness'][:])
    numThickness = 6

    sigma = thickness * numThickness
    dx = 900
    
    
    gaussKernel = Gaussian2DKernel(5)
    smoothedData = convolve(surface, gaussKernel)
    
    surfaceGrad = np.gradient(smoothedData, dx, dx)
    '''
    plt.subplot(2,1,1)
    plt.imshow(surfaceGrad[0])
    plt.colorbar()
    plt.subplot(2,1,2)
    plt.imshow(surfaceGrad[1])
    plt.colorbar()
    plt.show() 

    sys.exit()
    '''
    
    f.create_dataset('surfaceGradX', data = -1 *surfaceGrad[1])
    f.create_dataset('surfaceGradY', data = surfaceGrad[0])

    f.close()
           
    print time.time() - start

    

if __name__ == '__main__':
    main(sys.argv)
