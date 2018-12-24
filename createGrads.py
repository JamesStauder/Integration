import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
import time
def main(argv):
    start = time.time()
    f = h5py.File('data/GreenlandInBedCoord_V2.h5', 'r')
    keysToGrad = {'bed': [], 'smb': [], 't2m': [], 'thickness': [], 'surface': [],}
    dx = 900
    for key in keysToGrad:
        keysToGrad[key] = np.gradient(f[key], dx, dx)



    
    plt.subplot(4,1,1)
    plt.imshow(f['bed'][:])
    plt.colorbar()
    plt.subplot(4,1,2)
    plt.imshow((keysToGrad['bed'][1] **2 + keysToGrad['bed'][0]**2)**.5)
    plt.colorbar()
    plt.subplot(4,1,3)
    plt.imshow(f['surface'][:])
    plt.colorbar()
    plt.subplot(4,1,4)
    plt.imshow((keysToGrad['surface'][1] **2 + keysToGrad['surface'][0]**2)**.5)
    plt.colorbar()
    plt.show() 

    sys.exit()
    

        
    print time.time() - start

if __name__ == '__main__':
    main(sys.argv)
