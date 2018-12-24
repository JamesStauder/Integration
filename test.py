import multiprocessing
import numpy as np
import sys
import time
import h5py

from scipy.interpolate import RectBivariateSpline
from scipy.integrate import ode
from astropy.convolution import convolve, Gaussian2DKernel, Tophat2DKernel

def createInterp(data, xData, yData):
    yData = yData[::-1]
    return RectBivariateSpline(xData, yData, np.flipud(data).transpose())
class FlowIntegrator:
    def __init__(self, xVec, yVec):

        self.xVec = xVec
        self.yVec = yVec
        def rhs(t, u):
            x=u[0]
            y=u[1]
            d=u[2]

            vx = self.xVec(x,y)
            vy = self.yVec(x,y)
            vMag = np.sqrt(vx**2 + vy**2)
            return np.array([-vx / vMag, -vy / vMag, vMag])

        self.integrator = ode(rhs).set_integrator('vode', method = 'adams')

    def integrate(self, flowline):
        x0 = flowline[0][0]
        y0 = flowline[0][1]
        u0 = np.array([x0,y0,0.])
        self.integrator.set_initial_value(u0, 0.0)

        #Skipped some basic code to speed up
        dt = 1000.0

        # vx and vy at current location
        vx = self.xVec(x0, y0)
        vy = self.yVec(x0, y0)
        vMag = np.sqrt(vx**2 + vy**2)

        # x and y positions along flow line
        xs = [x0]
        ys = [y0]
        # Distance traveled
        ds = [0.0]
        # times
        ts = [0.]

        count = 0

        while self.integrator.successful() and count < len(flowline)-1:
            count = count + 1

            # Step forward
            u = self.integrator.integrate(self.integrator.t + dt)

            x = u[0]
            y = u[1]
            d = u[2]

            flowline[count] = [x, y]
        return flowline

def integ(index, xArray, returnDict):

    
    fileName = 'data/GreenlandInBedCoord_V2.h5'
    f = h5py.File(fileName, 'r')

    vxData = f['VX'][:]
    vyData = f['VY'][:]
    xData = f['x'][:]
    yData = f['y'][:]

    vxInt = createInterp(vxData, xData, yData)
    vyInt = createInterp(vyData, xData, yData)

    velIntegrator = FlowIntegrator(vxInt, vyInt)


    dx = 900
    surfData = f['surface'][:]
    gaussKernel = Gaussian2DKernel(5)
    smoothedSurf = convolve(surfData, gaussKernel)
    surfGrad = np.gradient(smoothedSurf, dx, dx)
    f.close()


    sxInt = createInterp(surfGrad[1], xData, yData)
    syInt = createInterp(surfGrad[0], xData, yData)
    surfIntegrator = FlowIntegrator(sxInt, syInt)


    num_cpus = multiprocessing.cpu_count()
    distMatrix = np.zeros((len(xArray), len(yData)))

    for row in range(0, len(xArray)):
        for col in range(0, len(yData)):
            print row
            velLine = [None] * 75
            velLine[0] = [xArray[row], yData[col]]

            surfLine = [None] * 75
            surfLine[0] = [xArray[row], yData[col]]

            try:
                velLine = velIntegrator.integrate(velLine)
                surfLine = surfIntegrator.integrate(surfLine)
                if None in velLine or None in surfLine:
                    distMatrix[row][col] = None
                    print "None"
                else:
                    distMatrix[row][col] = np.sqrt(
                    (velLine[-1][0] - surfLine[-1][0]) **2 +
                    (velLine[-1][1] - surfLine[-1][1]) **2
                        )
                    print "Success"
            except:
                distMatrix[row][col] = None
                print "Error"    
    returnDict[index] = distMatrix

    

if __name__ == '__main__':

    '''
    myArray = np.zeros((36,36))

    print cpu_count()
    for i in range(36):
        for j in range(36):
            myArray[i][j] = i*36 + j

    myArrays = []
    for i in range(3):
        myArrays.append(myArray[i*12:(i+1)*12,])

    p = Pool(5)
    print(p.map(f, 
        myArrays))
    '''

    dataFile = 'data/GreenlandInBedCoord_V2.h5'
    cpus = multiprocessing.cpu_count()
    f = h5py.File(dataFile, 'r')

    xData = f['x'][:]
    f.close()

    splitX = []
    for i in range(0, cpus-1):
        splitX.append(xData[i*(len(xData)/cpus):(i+1)*(len(xData)/cpus),])

    splitX.append(xData[cpus*(len(xData)/cpus):,])



    manager = multiprocessing.Manager()
    returnDict = manager.dict()
    jobs = []
    for i in range(len(splitX)):
        p = multiprocessing.Process(target=integ, args=(i, splitX[i],returnDict))
        jobs.append(p)
        p.start()
    for proc in jobs:
        proc.join()


    g = h5py.File('output.h5', 'w')
    for i in returnDict.keys():
        g.create_dataset(str(i), data = returnDict[i])
    g.close()






