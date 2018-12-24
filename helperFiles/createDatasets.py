import time
from classes.Dataset import *
from constants import *

'''
Function: createInitialDataSets
Argument list: 
Purpose: Create dictionary of datasets
Return types, values:
Dependencies:  h5py, Dataset Class
Creator: James Stauder
Date created: 1/31/18
Last edited: 1/31/18
'''


def createInitialDatasets():
    print "Creating data sets"
    t0 = time.time()

    datasetDict = {}

    dataFile = h5py.File(dataFileName, 'r')
    mapCoord['x1'] = len(dataFile['bed'][:][0])
    mapCoord['y1'] = len(dataFile['bed'][:])
    mapCoord['proj_x1'] = dataFile['x'][:][-1]
    mapCoord['proj_y1'] = dataFile['y'][:][-1]



    surfaceX = Dataset('surfaceGradX')
    surfaceY = Dataset('surfaceGradY')
    datasetDict['surfaceGradX'] = surfaceX
    datasetDict['surfaceGradY'] = surfaceY


    velocity = Dataset('velocity')
    datasetDict['velocity'] = velocity


    smb = Dataset('smb')
    datasetDict['smb'] = smb


    bed = Dataset('bed')
    datasetDict['bed'] = bed


    surface = Dataset('surface')
    datasetDict['surface'] = surface

    thickness = Dataset('thickness')
    datasetDict['thickness'] = thickness

    t2m = Dataset('t2m')
    datasetDict['t2m'] = t2m

    datasetDict['x'] = Dataset('x')
    datasetDict['y'] = Dataset('y')

    dataFile.close()

    print "Loaded all data sets in ", time.time() - t0, " seconds"
    return datasetDict
