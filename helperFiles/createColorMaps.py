import h5py
from classes.ColorBarAnchorWidget import *
from constants import *

def createColorMap(data):
    colorMapFile = h5py.File(cmFileName, 'r')
    data.colorData = colorMapFile[data.name][:]

    data.imageItem = pg.ImageItem(data.colorData)
    data.imageItem.setOpts(axisOrder='row-major')

    data.plotWidget = pg.PlotWidget()
    data.plotWidget.addItem(data.imageItem)
    data.plotWidget.setAspectLocked(True)
    data.plotWidget.invertY(True)

    '''
    data.colorMap = getCM(data.name)
    data.colorBar = getColorBar(data.name, data.colorMap)

    data.colorBarAnchorWidget = ColorBarAnchorWidget()
    data.colorBarAnchorWidget.hideAxis('left')
    data.colorBarAnchorWidget.hideAxis('bottom')
    data.colorBarAnchorWidget.addItem(data.colorBar)

    data.plotWidget.addItem(data.colorBarAnchorWidget)
    data.colorBarAnchorWidget.setFixedWidth(150)
    data.colorBarAnchorWidget.setFixedHeight(250)
    data.colorBarAnchorWidget.setAspectLocked(True)
    data.colorBarAnchorWidget.getViewBox().setRange(xRange=[-64.0, 114], yRange=[-15, 247], padding=0.0)
    data.colorBarAnchorWidget.invertY(True)
    data.colorBarAnchorWidget.setParentItem(data.plotWidget.getPlotItem())
    data.colorBarAnchorWidget.getViewBox().setMouseEnabled(x=False, y=False)
    data.colorBarAnchorWidget.anchor(itemPos=(1, 0), parentPos=(1, 0), offset=(-10, -10))
    
    colorMapFile.close()


def getCM(dataName):
    if dataName == 'velocity':
        minimum, maximum = 0.0, 12950.7
        c = [[[255, 255, 255], 0.01],  # zero is white
             [[0, 76, 153], 1],  # dark blue 1 - >10
             [[0, 0, 200], 10],
             [[0, 255, 255], 100],
             [[0, 153, 76], 1000],
             [[255, 255, 0], maximum],
             [[200, 100, 0], 0],
             [[255, 0, 0], 0],
             [[153, 0, 0], 0],
             [[120, 0, 0], 0]]

        logMax = np.log10(maximum - 0.01)
        div = logMax / (len(c) - 2)
        pos = [0, 0.01]
        for i in range(1, len(c) - 2):
            pos.append(np.power(10, i * div))
        pos.append(8000)
        return pg.ColorMap(pos, [row[0] for row in c], mode='byte')
    elif dataName == 'thickness':
        colors = [[255, 0, 255],
                  [255, 255, 255],
                  [6, 6, 94],
                  [0, 213, 253],
                  [255, 253, 3],
                  [138, 0, 0],
                  [110, 0, 0]]
        pos = [-0.00001, 0.0, 0.00001, 1000.0, 2000.0, 3000.0, 3500.0]
        return pg.ColorMap(pos, colors, mode='byte')

    # This might need to be redone or explained better
    elif dataName == 'bed':
        bedMin, bedMax = -5052.39510208, 3675.78314839

        a = [[0, 38, 115],
             [0, 38, 115],
             [0, 41, 117],
             [0, 41, 117],
             [0, 43, 122],
             [0, 43, 122],
             [0, 45, 128],
             [0, 45, 128],
             [0, 46, 130],
             [0, 46, 130],
             [0, 47, 135],
             [0, 47, 135],
             [0, 51, 140],
             [0, 51, 140],
             [0, 53, 145],
             [0, 53, 145],
             [0, 54, 148],
             [0, 54, 148],
             [0, 56, 153],
             [0, 56, 153],
             [0, 58, 158],
             [0, 58, 158],
             [0, 63, 163],
             [0, 63, 163],
             [0, 65, 168],
             [0, 65, 168],
             [0, 65, 171],
             [0, 65, 171],
             [0, 67, 176],
             [0, 67, 176],
             [0, 69, 181],
             [0, 69, 181],
             [0, 74, 186],
             [0, 74, 186],
             [0, 77, 191],
             [0, 77, 191],
             [0, 78, 194],
             [0, 78, 194],
             [0, 80, 199],
             [0, 80, 199],
             [0, 82, 204],
             [0, 82, 204],
             [0, 87, 209],
             [0, 87, 209],
             [0, 88, 212],
             [0, 88, 212],
             [0, 90, 217],
             [0, 90, 217],
             [0, 92, 222],
             [0, 92, 222],
             [0, 95, 227],
             [0, 95, 227],
             [0, 101, 232],
             [0, 101, 232],
             [0, 102, 235],
             [0, 102, 235],
             [0, 104, 240],
             [0, 104, 240],
             [0, 106, 245],
             [0, 106, 245],
             [0, 108, 250],
             [0, 108, 250],
             [0, 111, 255],
             [0, 111, 255],
             [5, 118, 255],
             [5, 118, 255],
             [10, 120, 255],
             [10, 120, 255],
             [18, 129, 255],
             [18, 129, 255],
             [25, 133, 255],
             [25, 133, 255],
             [31, 135, 255],
             [31, 135, 255],
             [36, 142, 255],
             [36, 142, 255],
             [41, 144, 255],
             [41, 144, 255],
             [48, 152, 255],
             [48, 152, 255],
             [56, 156, 255],
             [56, 156, 255],
             [61, 161, 255],
             [61, 161, 255],
             [66, 164, 255],
             [66, 164, 255],
             [74, 167, 255],
             [74, 167, 255],
             [82, 174, 255],
             [82, 174, 255],
             [87, 176, 255],
             [87, 176, 255],
             [92, 182, 255],
             [92, 182, 255],
             [97, 184, 255],
             [97, 184, 255],
             [105, 190, 255],
             [105, 190, 255],
             [112, 193, 255],
             [112, 193, 255],
             [117, 195, 255],
             [117, 195, 255],
             [122, 200, 255],
             [122, 200, 255],
             [128, 202, 255],
             [128, 202, 255],
             [138, 208, 255],
             [138, 208, 255],
             [143, 210, 255],
             [143, 210, 255],
             [148, 212, 255],
             [148, 212, 255],
             [153, 216, 255],
             [153, 216, 255],
             [161, 219, 255],
             [161, 219, 255],
             [168, 223, 255],
             [168, 223, 255],
             [173, 225, 255],
             [173, 225, 255],
             [179, 228, 255],
             [179, 228, 255],
             [184, 230, 255],
             [184, 230, 255],
             [191, 233, 255],
             [191, 233, 255]]
        lcArr = a
        div = float(np.abs(bedMin)) / float(len(lcArr) - 1)
        linPos = []
        for i in range(len(lcArr)):
            linPos.append(bedMin + i * div)



        # ELEVATION COLORMAP
        lcArr1 = [
            [0, 153, 0],
            [102, 255, 102],
            [255, 255, 0],
            [255, 0, 0],
            [153, 0, 0],
            [255, 255, 255]
        ]

        linPos1 = [0, 3678 / 5, (2 * 3678) / 5, (3 * 3678) / 5, (4 * 3678) / 5, (5 * 3678) / 5]
        bedColorArr = lcArr1
        bedPosArr = linPos1

        return pg.ColorMap(bedPosArr, bedColorArr, mode='byte')

    elif dataName == 'surface':
        colorArr = [
            [255, 255, 255],
            [0, 132, 53],
            [51, 204, 0],
            [244, 240, 113],
            [244, 189, 69],
            [178, 0, 0],
            [255, 255, 255],
        ]
        mx = 3677

        # meters
        elevationArr = [
            0,
            0.0001,
            .125 * mx,
            .25 * mx,
            .5 * mx,
            .75 * mx,
            mx,
        ]
        return pg.ColorMap(elevationArr, colorArr, mode='byte')

    # TODO redo this
    elif dataName == 'smb':
        colorsSMB = [
            [255, 0, 0],
            [255, 255, 255],
            [0, 0, 255]
        ]
        posSMB = [-11000, 0, 6061]
        return pg.ColorMap(posSMB, colorsSMB)

    elif dataName == 't2m':
        colorArr = [
            [160, 32, 240],  # purple
            [173, 216, 230],  # light blue
            [0, 255, 255],  # teal?
            [0, 255, 0],  # green
            [255, 255, 0],  # yellow
            [255, 165, 0],  # Orange
            [255, 0, 0],  # red
        ]

        # Kelvin
        tempArr = [
            240.0,
            243.0,
            247.0,
            250.0,
            253.0,
            257.0,
            260.0
        ]
        return pg.ColorMap(tempArr, colorArr, mode='byte')

def getColorBar(dataName, cm):
    colorbarHeight = 200
    colorbarWidth = 20
    if dataName == 'velocity':
        cbItem = LogColorBar(cm, colorbarWidth, colorbarHeight,
                             label='Velocity(m/yr)',
                             tick_labels=['0', '10', '100', '1,000', '8,000'],
                             ticks=[0, 10, 100, 1000, 8000])
        return cbItem

    elif dataName == 'bed':
        cbItem = ColorBar(cm, colorbarWidth, colorbarHeight,
                          label='Bed Ele.(m)',
                          tick_labels=['-5,000','-4,000','-3,000','-2,000','-1,000','0', '1,000', '2,000', '3,000', '3,500'],
                          ticks=[-5000,-4000,-3000,-2000,-1000,0, 1000, 2000, 3000, 3500])
        return cbItem
    elif dataName == 'surface':
        cbItem = ColorBar(cm, colorbarWidth, colorbarHeight,
                          label='Surface Ele.(m)',
                          tick_labels=['0', '500', '1,000', '1,500', '2,000', '2,500', '3,000', '3,500'],
                          ticks=[0, 500, 1000, 1500, 2000, 2500, 3000, 3500])
        return cbItem

    elif dataName == 'smb':
        cbItem = ColorBar(cm, colorbarWidth, colorbarHeight,
                          label='SMB(m)',
                          tick_labels=['-11', '0',  '6'],
                          ticks=[-11000, 0, 6000],
                          name='smb')
        return cbItem

    elif dataName == 'thickness':
        cbItem = ColorBar(cm, colorbarWidth, colorbarHeight,
                          label='Thickness(m)',
                          tick_labels=['0', '500', '1000', '1500', '2000', '2500', '3000', '3500'],
                          ticks=[0, 500, 1000, 1500, 2000, 2500, 3000, 3500])
        return cbItem

    elif dataName == 't2m':
        cbItem = ColorBar(cm, colorbarWidth, colorbarHeight,
                          label='T2M(K)',
                          tick_labels=['240', '243', '247', '250', '253', '257', '260'],
                          ticks=[240, 243, 247, 250, 253, 257, 260],
						  name="T2M")
        return cbItem

'''