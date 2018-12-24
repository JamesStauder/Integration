from PyQt4.QtGui import *
from pyqtgraph.Qt import QtGui
import math
from FlowIntegrator import *
from Dataset import *
from Marker import *

from pyproj import Proj
import time

'''
Class: MainWindow
Argument list:
Purpose: Create main window for GUI. Has many functions based on what user does
Return types, values:
Dependencies: pyQT, dolfin, math
Creator: James Stauder
Date created: 1/31/18
Last edited: 5/29/18
'''


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        self.setWindowTitle("Greenland")
        self.setMinimumHeight(1000)
        self.setMinimumWidth(1200)

        self.centralWidget = QtGui.QWidget()
        self.setCentralWidget(self.centralWidget)
        self.mainLayout = QtGui.QHBoxLayout()
        self.centralWidget.setLayout(self.mainLayout)

        # index of current map
        self.currentMap = 0

        # marker selected variables
        self.isMarkerSelected = False
        self.whichMarkerSelected = None
        self.selectedMarkerPosition = None
        self.whichIndexOfFlowlineSelected = None

        # Flowline information
        self.flowlineDistance = 100000
        self.lengthOfFlowline = 1
        self.flowlines = []
        self.flowlineMarkers = []
        self.integratorPerMarker = 10

        '''
        Side widget with button
        '''
        self.maxWidth = 300

        self.buttonBoxWidget = QtGui.QWidget()
        self.buttonBox = QtGui.QVBoxLayout()
        self.buttonBoxWidget.setLayout(self.buttonBox)

        self.mapList = QtGui.QComboBox()
        self.maps = ['Velocity', 'Bed', 'Surface', 'SMB', 'Thickness', 't2m']
        self.mapList.addItems(self.maps)
        self.mapList.setMaximumWidth(self.maxWidth)
        self.buttonBox.addWidget(self.mapList)

        self.spatialResolutionWidget = QtGui.QWidget()
        self.spatialResolutionLayout = QtGui.QHBoxLayout()
        self.spatialResolutionWidget.setLayout(self.spatialResolutionLayout)
        self.spatialResolutionLabel = QtGui.QLabel('Spatial Resolution(m)')
        self.spatialResolutionLineEdit = QtGui.QLineEdit('1000')
        self.spatialResolutionLayout.addWidget(self.spatialResolutionLabel)
        self.spatialResolutionLayout.addWidget(self.spatialResolutionLineEdit)
        self.buttonBox.addWidget(self.spatialResolutionWidget)

        self.distanceWidget = QtGui.QWidget()
        self.distanceLayout = QtGui.QHBoxLayout()
        self.distanceWidget.setLayout(self.distanceLayout)
        self.distanceLabel = QtGui.QLabel('distance(km)')
        self.distanceLineEdit = QtGui.QLineEdit('100')
        self.spatialResolutionLayout.addWidget(self.distanceLabel)
        self.spatialResolutionLayout.addWidget(self.distanceLineEdit)
        self.buttonBox.addWidget(self.distanceWidget)

        self.upButton = QRadioButton('Integrate Up')
        self.downButton = QRadioButton('Integrate Down')
        self.upButton.setChecked(True)
        self.buttonBox.addWidget(self.upButton)
        self.buttonBox.addWidget(self.downButton)

        self.resetButton = QtGui.QPushButton('Reset')
        self.resetButton.setEnabled(True)
        self.resetButton.setMaximumWidth(self.maxWidth)
        self.buttonBox.addWidget(self.resetButton)

        self.distButton = QtGui.QPushButton('calcDist')
        self.distButton.setEnabled(True)
        self.distButton.setMaximumWidth(self.maxWidth)
        self.buttonBox.addWidget(self.distButton)

        self.latLongWidget = QtGui.QWidget()
        self.latLongLayout = QtGui.QHBoxLayout()
        self.latLongWidget.setLayout(self.latLongLayout)
        self.latLabel = QtGui.QLabel('Lat')
        self.latLineEdit = QtGui.QLineEdit('')
        self.longLabel = QtGui.QLabel('Long')
        self.longLineEdit = QtGui.QLineEdit('')
        self.latLongLayout.addWidget(self.latLabel)
        self.latLongLayout.addWidget(self.latLineEdit)
        self.latLongLayout.addWidget(self.longLabel)
        self.latLongLayout.addWidget(self.longLineEdit)
        self.latLongButton = QtGui.QPushButton('Find Lat Long')

        self.buttonBox.addWidget(self.latLongWidget)
        self.buttonBox.addWidget(self.latLongButton)

        self.textOut = QtGui.QTextBrowser()
        self.textOut.setMaximumWidth(self.maxWidth)
        self.buttonBox.addWidget(self.textOut)

        self.leftSideWidget = QtGui.QWidget()
        self.leftSide = QtGui.QVBoxLayout()
        self.leftSideWidget.setLayout(self.leftSide)

        self.imageItemContainer = QtGui.QStackedWidget()

        self.leftSide.addWidget(self.imageItemContainer)

        self.mainLayout.addWidget(self.leftSideWidget)
        self.mainLayout.addWidget(self.buttonBoxWidget)

        self.buttonBoxWidget.setMaximumWidth(self.maxWidth + 12)

        self.connectButtons()

    '''
    Function: addToImageItemContainer
    Argument list: datasetDict
    Purpose: add the different dataset widgets to the imageItemContainer
    Return types, values: 
    Dependencies: 
    Creator: James Stauder
    Date created: 2/25/18
    Last edited: 3/5/18
    '''

    def addToImageItemContainer(self, datasetDict):
        self.imageItemContainer.addWidget(datasetDict['velocity'].plotWidget)
        self.imageItemContainer.setCurrentWidget(datasetDict['velocity'].plotWidget)
        self.imageItemContainer.addWidget(datasetDict['bed'].plotWidget)
        self.imageItemContainer.addWidget(datasetDict['surface'].plotWidget)
        self.imageItemContainer.addWidget(datasetDict['thickness'].plotWidget)
        self.imageItemContainer.addWidget(datasetDict['t2m'].plotWidget)
        self.imageItemContainer.addWidget(datasetDict['smb'].plotWidget)

    '''
    Function: changeMap
    Argument list: 
        index: index of which map to use
    Purpose: Changes the map to a different colormap
    Return types, values: 
    Dependencies: 
    Creator: James Stauder
    Date created: 2/25/18
    Last edited: 3/5/18
    '''

    def changeMap(self, index):
        vr = self.imageItemContainer.currentWidget().getPlotItem().getViewBox().viewRange()
        indexToDatasetDict = {
            0: 'velocity',
            1: 'bed',
            2: 'surface',
            3: 'smb',
            4: 'thickness',
            5: 't2m'}
        if index != self.currentMap:
            oldMap = self.currentMap
            self.currentMap = index

        self.imageItemContainer.setCurrentWidget(self.datasetDict[indexToDatasetDict[self.currentMap]].plotWidget)
        self.datasetDict[indexToDatasetDict[self.currentMap]].imageItem.hoverEvent = self.mouseMove
        self.datasetDict[indexToDatasetDict[self.currentMap]].imageItem.mouseClickEvent = self.mouseClick

        self.datasetDict[indexToDatasetDict[self.currentMap]].plotWidget.getPlotItem().getViewBox().setRange(
            xRange=vr[0],
            yRange=vr[1],
            padding=0.0)
        for line in self.flowlineMarkers:
            for marker in line:
                marker.plotWidget = self.datasetDict[indexToDatasetDict[self.currentMap]]
                self.datasetDict[indexToDatasetDict[oldMap]].plotWidget.removeItem(marker.cross[0])
                self.datasetDict[indexToDatasetDict[oldMap]].plotWidget.removeItem(marker.cross[1])
                self.datasetDict[indexToDatasetDict[self.currentMap]].plotWidget.addItem(marker.cross[0])
                self.datasetDict[indexToDatasetDict[self.currentMap]].plotWidget.addItem(marker.cross[1])

                if marker.lines[0]:
                    self.datasetDict[indexToDatasetDict[self.currentMap]].plotWidget.addItem(marker.lines[0])

    '''
     Function: mouseClick
     Argument list: 
        e: event trigger from mouse being clicked
     Purpose: 
        Create a new flowline or move a previous flowline
     Return types, values: None
     Dependencies: None
     Creator: James Stauder
     Date created: 2/25/18
     Last edited: 5/29/18
     '''

    def mouseClick(self, e):

        # If no marker is selected
        if self.isMarkerSelected is False:

            self.upButton.setEnabled(False)
            self.downButton.setEnabled(False)
            if self.downButton.isChecked():
                self.velocityIntegrator.direction = 1
                self.surfaceIntegrator.direction = 1
            else:
                self.velocityIntegrator.direction = -1
                self.surfaceIntegrator.direction = -1

            # Checks to see only if first marker in each flowline is detected.
            for i in range(len(self.flowlineMarkers)):
                if self.flowlineMarkers[i][0].checkClicked(e.pos()):
                    self.isMarkerSelected = True
                    self.whichMarkerSelected = self.flowlineMarkers[i][0]
                    self.selectedMarkerPosition = [i, 0]
                    self.displayMarkerVariables()
                    self.whichIndexOfFlowlineSelected = [i, 0]
                    break

            # If no marker selected previously or currently create new flowline.
            if self.isMarkerSelected is False:
                self.spatialResolutionLineEdit.setReadOnly(True)
                self.distanceLineEdit.setReadOnly(True)
                self.flowlineDistance = int(self.distanceLineEdit.text()) * 1000
                self.lengthOfFlowline = int(self.flowlineDistance / float(self.spatialResolutionLineEdit.text()))
                self.integratorPerMarker = int(math.ceil(10000 / (float(self.spatialResolutionLineEdit.text()))))
                xClickPosition = e.pos().x()
                yClickPosition = e.pos().y()


                dx, dy = colorToProj(xClickPosition, yClickPosition)

                # Create new flowline
                velFlowline = [None] * self.lengthOfFlowline
                velFlowline[0] = [dx, dy]

                velFlowline = self.velocityIntegrator.integrate(dx, dy, velFlowline, 0,
                                                                float(self.spatialResolutionLineEdit.text()))

                if None in velFlowline:
                    print "Integration Error. Try Again"
                    return

                # Create a flowline of markers spaced out based on the IntegratorPerMarker
                newFlowlineMarkers = velFlowline[::self.integratorPerMarker]

                for i in range(len(newFlowlineMarkers)):
                    dx = newFlowlineMarkers[i][0]
                    dy = newFlowlineMarkers[i][1]
                    cx, cy = colorCoord(dx, dy)
                    newFlowlineMarkers[i] = Marker(cx, cy, dx, dy, self.imageItemContainer.currentWidget())

                self.displayMarkers(newFlowlineMarkers)

                self.flowlines.append(velFlowline)
                self.flowlineMarkers.append(newFlowlineMarkers)

                # Create new flowline for bed
                dx, dy = colorToProj(xClickPosition, yClickPosition)
                surfFlowline = [None] * self.lengthOfFlowline
                surfFlowline[0] = [dx, dy]

                surfFlowline = self.surfaceIntegrator.integrate(dx, dy, surfFlowline, 0,
                                                                float(self.spatialResolutionLineEdit.text()))

                if None in surfFlowline:
                    print "Integration Error. Try Again"
                    return

                # Create a flowline of markers spaced out based on the IntegratorPerMarker
                newFlowlineMarkers = surfFlowline[::self.integratorPerMarker]

                for i in range(len(newFlowlineMarkers)):
                    dx = newFlowlineMarkers[i][0]
                    dy = newFlowlineMarkers[i][1]
                    cx, cy = colorCoord(dx, dy)
                    newFlowlineMarkers[i] = Marker(cx, cy, dx, dy, self.imageItemContainer.currentWidget(), pen=whitePlotPen)

                self.displayMarkers(newFlowlineMarkers, pen=skinnyWhitePlotPen)

                self.flowlines.append(surfFlowline)
                self.flowlineMarkers.append(newFlowlineMarkers)

                lastPoints = velFlowline[-1], surfFlowline[-1]


                velFlowline = np.asarray(velFlowline)
                surfFlowline = np.asarray(surfFlowline)

                distance = sqrt(
                    (velFlowline[:,0] - surfFlowline[:,0]) **2 +
                    (velFlowline[:,1] - surfFlowline[:,1]) **2
                )

                print distance



        # Release the marker that was previously held
        else:
            self.isMarkerSelected = False
            self.whichMarkerSelected = None
            self.textOut.clear()

    '''
    Function: mouseMove
    Argument list: 
    Purpose: This function is used to move the marker that is selected and create a new integration path. 
    Return types, values: 
    Dependencies: 
    Creator: James Stauder
    Date created: 2/25/18
    Last edited: 3/9/18
    TODO: 
        This can be a bit confusing to read. The code is kind of wordy. We could possibly redraw flowline with the 
        display Markers function but that would require some changes to the Markers function to take an index.
    '''

    def mouseMove(self, e):

        if self.isMarkerSelected:

            # change the x , y values of the marker at the selected index
            xPositionOfCursor = e.pos().x()
            yPositionOfCursor = e.pos().y()
            self.whichMarkerSelected.cx = xPositionOfCursor
            self.whichMarkerSelected.cy = yPositionOfCursor
            self.whichMarkerSelected.updateCross()

            # change the x, y values of the flowline at the selected index
            whichFlowlineSelected = self.whichIndexOfFlowlineSelected[0]
            indexSelected = self.whichIndexOfFlowlineSelected[1]
            self.flowlines[whichFlowlineSelected][indexSelected] = [self.whichMarkerSelected.dx,
                                                                    self.whichMarkerSelected.dy]

            self.flowlines[whichFlowlineSelected] = self.velocityIntegrator.integrate(
                self.whichMarkerSelected.dx, self.whichMarkerSelected.dy,
                self.flowlines[whichFlowlineSelected], indexSelected,
                float(self.spatialResolutionLineEdit.text()))

            # Remove every marker past the one we selected
            for i in range(self.selectedMarkerPosition[1] + 1, len(self.flowlineMarkers[0])):
                self.imageItemContainer.currentWidget().removeItem(
                    self.flowlineMarkers[self.selectedMarkerPosition[0]][i])

                # get the flowline position of the new marker
                newPosition = self.flowlines[whichFlowlineSelected][i * self.integratorPerMarker]
                cx, cy = colorCoord(newPosition[0], newPosition[1])

                # Create new marker with new data
                self.flowlineMarkers[self.selectedMarkerPosition[0]][i] = Marker(
                    cx, cy, newPosition[0], newPosition[1],
                    self.imageItemContainer.currentWidget())
            # This section redraws the new markers
            for i in range(self.selectedMarkerPosition[1] + 1, len(self.flowlineMarkers[0])):
                self.imageItemContainer.currentWidget().addItem(
                    self.flowlineMarkers[self.selectedMarkerPosition[0]][i].getCross()[0])
                self.imageItemContainer.currentWidget().addItem(
                    self.flowlineMarkers[self.selectedMarkerPosition[0]][i].getCross()[1])

                xa = [self.flowlineMarkers[self.selectedMarkerPosition[0]][i - 1].cx,
                      self.flowlineMarkers[self.selectedMarkerPosition[0]][i].cx]
                ya = [self.flowlineMarkers[self.selectedMarkerPosition[0]][i - 1].cy,
                      self.flowlineMarkers[self.selectedMarkerPosition[0]][i].cy]
                self.flowlineMarkers[self.selectedMarkerPosition[0]][i].setLine(
                    pg.PlotDataItem(xa, ya, connect='all', pen=skinnyBlackPlotPen), 0)
                self.flowlineMarkers[self.selectedMarkerPosition[0]][i - 1].setLine(
                    self.flowlineMarkers[self.selectedMarkerPosition[0]][i].lines[0], 1)

                self.imageItemContainer.currentWidget().addItem(
                    self.flowlineMarkers[self.selectedMarkerPosition[0]][i].lines[0])

            self.displayMarkerVariables()

            # Connect lines between marker selected and previous marker
            if self.whichMarkerSelected.lines[0] is not None:
                self.whichMarkerSelected.lines[0].setData(
                    [self.whichMarkerSelected.lines[0].getData()[0][0], self.whichMarkerSelected.cx],
                    [self.whichMarkerSelected.lines[0].getData()[1][0], self.whichMarkerSelected.cy])

    '''
    Function: displayMarkers
    Argument list: 
        flowline: flowline in which to display
    Purpose: Takes a flowline of markers and displays them on the gui
    Return types, values: None
    Dependencies: None
    Creator: James Stauder
    Date created: 3/2/18
    Last edited: 3/2/18
    '''
    def displayMarkers(self, flowline, pen=skinnyBlackPlotPen):

        # Add first marker. This needs to be peeled because the for loop
        # connects the markers backwards
        self.imageItemContainer.currentWidget().addItem(flowline[0].getCross()[0])
        self.imageItemContainer.currentWidget().addItem(flowline[0].getCross()[1])

        for i in range(1, len(flowline)):
            self.imageItemContainer.currentWidget().addItem(flowline[i].getCross()[0])
            self.imageItemContainer.currentWidget().addItem(flowline[i].getCross()[1])

            xValuesOfMarkers = [flowline[i - 1].cx, flowline[i].cx]
            yValuesOfMarkers = [flowline[i - 1].cy, flowline[i].cy]

            # Create lines from each marker
            flowline[i].setLine(
                pg.PlotDataItem(xValuesOfMarkers, yValuesOfMarkers, connect='all', pen=pen), 0)
            flowline[i - 1].setLine(flowline[i].lines[0], 1)

            self.imageItemContainer.currentWidget().addItem(flowline[i].lines[0])

    '''
    Function: displayMarkerVariables
    Argument list: None
    Purpose: Displays the marker variables of the marker selected
    Return types, values: None
    Dependencies: Marker to be selected
    Creator: James Stauder
    Date created: 2/25/18
    Last edited: 2/25/18
    '''
    def displayMarkerVariables(self):
        self.textOut.clear()
        selectedMarkerX = self.whichMarkerSelected.dx
        selectedMarkerY = self.whichMarkerSelected.dy

        self.textOut.append(str((self.whichMarkerSelected.dx, self.whichMarkerSelected.dy)))

        for x in self.maps:
            stringOut = str(self.datasetDict[x.lower()].getInterpolatedValue(selectedMarkerX, selectedMarkerY))
            self.textOut.append(x + ": " + stringOut[2:-2])

        self.textOut.append(
            'SurfX: ' + str(self.datasetDict['surfaceGradX'].getInterpolatedValue(selectedMarkerX, selectedMarkerY))[
                        2:-2])
        self.textOut.append(
            'SurfY: ' + str(self.datasetDict['surfaceGradY'].getInterpolatedValue(selectedMarkerX, selectedMarkerY))[
                        2:-2])

    '''
    Function: createIntegrator
    Argument list: Nones
    Purpose: Create integrator class. This will allow us to integrate up the ice flow
    Return types, values: None
    Dependencies: None
    Creator: James Stauder
    Date created: 2/5/18
    Last edited: 2/5/18
    '''
    def createIntegrator(self):
        vx = Dataset('VX')
        vy = Dataset('VY')
        self.velocityIntegrator = FlowIntegrator(vx, vy)
        bx = Dataset('surfaceGradX')
        by = Dataset('surfaceGradY')
        self.surfaceIntegrator = FlowIntegrator(bx,by)

    def reset(self):
        del self.flowlines[:]
        del self.flowlineMarkers[:]
        for x in self.datasetDict:
            self.datasetDict[x].pathData = None
        self.spatialResolutionLineEdit.setReadOnly(False)
        self.distanceLineEdit.setReadOnly(False)
        self.upButton.setEnabled(True)
        self.downButton.setEnabled(True)

    '''
    Function: connectButtons
    Argument list: None
    Purpose: connect the buttons of the gui to various functions
    Return types, values: None
    Dependencies: None
    Creator: James Stauder
    Date created: 2/5/18
    Last edited: 2/5/18
    '''

    def connectButtons(self):
        self.mapList.currentIndexChanged.connect(self.changeMap)
        self.resetButton.clicked.connect(self.reset)
        self.latLongButton.clicked.connect(self.markLatLong)
        self.distButton.clicked.connect(self.calcDistance)


    def calcDistance(self):
        xPoints = self.datasetDict['x'].data
        yPoints = self.datasetDict['y'].data

        intUpDistance = 75 #KM
        distMatrix = np.zeros((len(xPoints), len(yPoints)))


        for x in range(0, len(xPoints)):
            for y in range(0, len(yPoints)):
                try:
                    velIntLine = [None] * intUpDistance
                    velIntLine[0] = [xPoints[x], yPoints[y]]
                    velIntLine = self.velocityIntegrator.integrate(xPoints[x], yPoints[y], velIntLine, 0, intUpDistance)

                    surfIntLine = [None] * intUpDistance
                    surfIntLine[0] = [xPoints[x], yPoints[y]]
                    surfIntLine = self.surfaceIntegrator.integrate(xPoints[x], yPoints[y], surfIntLine, 0, intUpDistance)

                    if None in velIntLine or None in surfIntLine:
                        distMatrix[x][y] = None
                    else:
                        distMatrix[x][y] = sqrt(
                            (velIntLine[-1][0] - surfIntLine[-1][0])**2 +
                            (velIntLine[-1][1] - surfIntLine[-1][1])**2
                        )
                except:
                    distMatrix[x][y] = None
        f = h5py.File('output.h5', 'w')
        f.create_dataset('distance', data=distMatrix)
        f.close()
        print distMatrix
    def markLatLong(self):
        self.rosieButton.setEnabled(True)
        self.downButton.setEnabled(False)
        self.upButton.setEnabled(False)

        if (len(self.flowlines) < 2):
            self.spatialResolutionLineEdit.setReadOnly(True)
            self.distanceLineEdit.setReadOnly(True)
            self.flowlineDistance = int(self.distanceLineEdit.text()) * 1000
            self.lengthOfFlowline = int(self.flowlineDistance / float(self.spatialResolutionLineEdit.text()))
            self.integratorPerMarker = int(math.ceil(10000 / (float(self.spatialResolutionLineEdit.text()))))

            dx, dy = self.translateLatLong(
                float(self.latLineEdit.text()),
                float(self.longLineEdit.text())
            )

            # Create new flowline
            newFlowline = []
            for x in range(0, self.lengthOfFlowline):
                newFlowline.append(None)
            newFlowline[0] = [dx, dy]

            if self.downButton.isChecked():
                self.velocityIntegrator.direction = 1
            newFlowline = self.velocityIntegrator.integrate(dx, dy, newFlowline, 0,
                                                            float(self.spatialResolutionLineEdit.text()))

            if None in newFlowline:
                print "Integration Error. Try Again"
                return

            # Create a flowline of markers spaced out based on the IntegratorPerMarker
            newFlowlineMarkers = newFlowline[::self.integratorPerMarker]

            for i in range(len(newFlowlineMarkers)):
                dx = newFlowlineMarkers[i][0]
                dy = newFlowlineMarkers[i][1]
                cx, cy = colorCoord(dx, dy)
                newFlowlineMarkers[i] = Marker(cx, cy, dx, dy, self.imageItemContainer.currentWidget())

            self.displayMarkers(newFlowlineMarkers)
            self.flowlines.append(newFlowline)
            self.flowlineMarkers.append(newFlowlineMarkers)

    def translateLatLong(self, lat, long):
        proj = Proj(init='epsg:3413')

        return proj(long, lat)
