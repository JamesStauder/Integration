import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
import time

def main(argv):
    f = h5py.File('output.h5', 'r')
    totalArray = []

    for key in range(15):
        for array in f[str(key)][:]:
            totalArray.append(array)


    
    plt.imshow(totalArray)
    plt.colorbar()
    plt.show()
    f.close()
    sys.exit()
    

if __name__ == '__main__':
    main(sys.argv)
