import director.vtkAll as vtk
import director.visualization as vis
import director.objectmodel as om
from director.debugVis import DebugData
from director.consoleapp import ConsoleApp
from director.timercallback import TimerCallback
from director import applogic
from director import screengrabberpanel
from director import cameracontrolpanel

from director import transformUtils
import numpy as np
import time
import scipy.integrate as integrate
import argparse
import matplotlib.pyplot as plt
import shelve
import math
from PythonQt import QtCore, QtGui

from world import World
from car import CarPlant
from sensor import SensorObj
from sensorApproximator import SensorApproximatorObj
from controller import ControllerObj




class Simulator(object):


    def __init__(self, mode, percentObsDensity=20, endTime=40, nonRandomWorld=False,
                 circleRadius=0.7, worldScale=1.0, autoInitialize=True, verbose=True,
                 numCars=1):
        self.verbose = verbose
        self.startSimTime = time.time()
        self.collisionThreshold = 0.2
        self.randomSeed = 5
        self.Sensor_rayLength = 8

        self.percentObsDensity = percentObsDensity
        self.defaultControllerTime = 1000
        self.nonRandomWorld = nonRandomWorld
        self.circleRadius = circleRadius
        self.worldScale = worldScale
        self.mode = mode
        self.numCars =1 
        # create the visualizer object
        self.app = ConsoleApp()
        self.view = self.app.createView(useGrid=False)
        self.initializeOptions()
        self.initializeColorMap()
        if autoInitialize:
            self.initialize()

    def initializeOptions(self):
        self.options = dict()

        self.options['World'] = dict()
        self.options['World']['obstaclesInnerFraction'] = 0.98
        self.options['World']['randomSeed'] = 40
        self.options['World']['percentObsDensity'] = 30
        self.options['World']['nonRandomWorld'] = True
        self.options['World']['circleRadius'] = 1.0
        self.options['World']['scale'] = 10

        self.options['Sensor'] = dict()
        self.options['Sensor']['rayLength'] = 20
        self.options['Sensor']['numRays'] = 21


        self.options['Car'] = dict()
        self.options['Car']['velocity'] = 20

        self.options['dt'] = 0.05

        self.options['runTime'] = dict()
        self.options['runTime']['defaultControllerTime'] = 100


    def setDefaultOptions(self):

        defaultOptions = dict()

        defaultOptions['World'] = dict()
        defaultOptions['World']['obstaclesInnerFraction'] = 0.98
        defaultOptions['World']['randomSeed'] = 40
        defaultOptions['World']['percentObsDensity'] = 30
        defaultOptions['World']['nonRandomWorld'] = True
        defaultOptions['World']['circleRadius'] = 1.75
        defaultOptions['World']['scale'] = 2.5


        defaultOptions['Sensor'] = dict()
        defaultOptions['Sensor']['rayLength'] = 20
        defaultOptions['Sensor']['numRays'] = 41


        defaultOptions['Car'] = dict()
        defaultOptions['Car']['velocity'] = 20

        defaultOptions['dt'] = 0.05
        defaultOptions['runTime'] = dict()
        defaultOptions['runTime']['defaultControllerTime'] = 100

        for k in defaultOptions:
            self.options.setdefault(k, defaultOptions[k])

        for k in defaultOptions:
            if not isinstance(defaultOptions[k], dict):
                continue
            for j in defaultOptions[k]:
                self.options[k].setdefault(j, defaultOptions[k][j])

    def initializeColorMap(self):
        self.colorMap = dict()
        self.colorMap['default'] = [0,1,0]

    def initialize(self):

        self.dt = self.options['dt']
        self.controllerTypeOrder = ['default']

        self.setDefaultOptions()

        self.Sensor = SensorObj(rayLength=self.options['Sensor']['rayLength'],
                                numRays=self.options['Sensor']['numRays'])

        self.SensorApproximator = SensorApproximatorObj(numRays=self.options['Sensor']['numRays'], circleRadius=self.options['World']['circleRadius'], )

        om.removeFromObjectModel(om.findObjectByName('world'))
        if self.mode=='Obs':
            self.world, self.goalX, self.goalY = World.buildCircleWorld(percentObsDensity=self.options['World']['percentObsDensity'],
                                            circleRadius=self.options['World']['circleRadius'],
                                            nonRandom=self.options['World']['nonRandomWorld'],
                                            scale=self.options['World']['scale'],
                                            randomSeed=self.options['World']['randomSeed'],
                                            obstaclesInnerFraction=self.options['World']['obstaclesInnerFraction'])
            self.EgoCarController = ControllerObj(self.Sensor, self.SensorApproximator, self.mode, self.goalX, self.goalY)
        else:
            print 'simulator mode error'
        
        self.EgoCar = CarPlant(controller=self.EgoCarController,
                            velocity=self.options['Car']['velocity'])
        self.EgoCarController.addingCar(self.EgoCar)
        self.EgoCarController.initializeVelocity(self.EgoCar.v)


        if self.numCars>1:
            # implement agent cars 
            pass

        # create the things needed for simulation
        om.removeFromObjectModel(om.findObjectByName('robot'))

        self.robot, self.frame = World.buildRobot(numCars=self.numCars)

        self.locator = World.buildCellLocator(self.world.visObj.polyData)
        self.Sensor.setLocator(self.locator)
        self.frame = self.robot.getChildFrame()
        self.frame.setProperty('Scale', 3)
        #self.frame.setProperty('Visible', False)
        #self.frame.setProperty('Edit', True)
        self.frame.widget.HandleRotationEnabledOff()
        rep = self.frame.widget.GetRepresentation()
        rep.SetTranslateAxisEnabled(2, False)
        rep.SetRotateAxisEnabled(0, False)
        rep.SetRotateAxisEnabled(1, False)

        self.defaultControllerTime = self.options['runTime']['defaultControllerTime']

        self.EgoCar.setFrame(self.frame)
        print "Finished initialization"


    def runSingleSimulation(self, controllerType='default', simulationCutoff=None):
        self.setRandomCollisionFreeInitialState()

        currentCarState = np.copy(self.EgoCar.state)
        nextCarState = np.copy(self.EgoCar.state)
        currentAngleState = np.copy(self.EgoCar.angles)
        nextAngleState=np.copy(self.EgoCar.angles)
        self.setRobotFrameState(currentCarState[0], currentCarState[1], currentCarState[2])
        currentRaycast = self.Sensor.raycastAll(self.frame)
        nextRaycast = np.zeros(self.Sensor.numRays)

        # record the reward data
        runData = dict()
        startIdx = self.counter
        while (self.counter < self.numTimesteps - 1):
            idx = self.counter
            currentTime = self.t[idx]
            self.stateOverTime[idx,:] = currentCarState
            x = self.stateOverTime[idx,0]
            y = self.stateOverTime[idx,1]
            theta = self.stateOverTime[idx,2]

            self.setRobotFrameState(x,y,theta)
            currentRaycast = self.Sensor.raycastAll(self.frame)
            self.raycastData[idx,:] = currentRaycast
            S_current = (currentCarState, currentRaycast)
          

            if controllerType not in self.colorMap.keys():
                print
                raise ValueError("controller of type " + controllerType + " not supported")


            if controllerType in ["default", "defaultRandom"]:
                controlInput, controlAngle, controlInputIdx = self.EgoCarController.computeControlInput(currentCarState,
                                                                            currentTime, self.frame,
                                                                            raycastDistance=currentRaycast,
                                                                            randomize=False)



            self.angleOverTime[idx,:] = controlAngle

            phi = self.angleOverTime[idx,0]
            alpha = self.angleOverTime[idx,1]
            u = self.angleOverTime[idx,2]


            self.controlInputData[idx] = controlInput

            nextCarState = self.EgoCar.simulateOneStep(controlInput=controlInput, dt=self.dt)        
            x = nextCarState[0]
            y = nextCarState[1]
            theta = nextCarState[2]


            self.setRobotFrameState(x,y,theta)
            nextRaycast = self.Sensor.raycastAll(self.frame)


            # Compute the next control input
            S_next = (nextCarState, nextRaycast)

            if controllerType in ["default", "defaultRandom"]:
                nextControlInput, nextAngle, nextControlInputIdx = self.EgoCarController.computeControlInput(nextCarState,
                                                                            currentTime, self.frame,
                                                                            raycastDistance=nextRaycast,
                                                                            randomize=False)
            

            #bookkeeping

            currentCarState = nextCarState
            nextAngleState=nextAngle
            currentRaycast = nextRaycast
            self.counter+=1

            # break if we are in collision
            if self.checkInCollision(nextRaycast):
                if self.verbose: print "Had a collision, terminating simulation"
                break

            if self.counter >= simulationCutoff:
                break



        # fill in the last state by hand
        self.stateOverTime[self.counter,:] = currentCarState
        self.angleOverTime[self.counter,:] = currentAngleState
        self.raycastData[self.counter,:] = currentRaycast


        # this just makes sure we don't get stuck in an infinite loop.
        if startIdx == self.counter:
            self.counter += 1

        return runData

    def setNumpyRandomSeed(self, seed=1):
        np.random.seed(seed)

    def runBatchSimulation(self, endTime=None, dt=0.05):
        # for use in playback
        self.dt = self.options['dt']
        self.endTime = self.defaultControllerTime # used to be the sum of the other times as well
        self.t = np.arange(0.0, self.endTime, dt)
        maxNumTimesteps = np.size(self.t)
        self.stateOverTime = np.zeros((maxNumTimesteps+1, 3))
        self.angleOverTime = np.zeros((maxNumTimesteps+1, 3))
        
        self.raycastData = np.zeros((maxNumTimesteps+1, self.Sensor.numRays))
        self.controlInputData = np.zeros(maxNumTimesteps+1)
        self.numTimesteps = maxNumTimesteps

        self.controllerTypeOrder = ['default']
        self.counter = 0
        self.simulationData = []
    
        self.initializeStatusBar()

        self.idxDict = dict()
        numRunsCounter = 0


        self.idxDict['default'] = self.counter
        loopStartIdx = self.counter
        simCutoff = min(loopStartIdx + self.defaultControllerTime/dt, self.numTimesteps)
        
        while ((self.counter - loopStartIdx < self.defaultControllerTime/dt) and self.counter < self.numTimesteps-1):
            self.printStatusBar()
            startIdx = self.counter
            runData = self.runSingleSimulation(controllerType='default',
                                               simulationCutoff=simCutoff)
            runData['startIdx'] = startIdx
            runData['controllerType'] = "default"
            runData['duration'] = self.counter - runData['startIdx']
            runData['endIdx'] = self.counter
            runData['runNumber'] = numRunsCounter
            numRunsCounter+=1
            self.simulationData.append(runData)

        # BOOKKEEPING
        # truncate stateOverTime, raycastData, controlInputs to be the correct size
        self.numTimesteps = self.counter + 1
        self.angleOverTime = self.angleOverTime[0:self.counter+1, :]
        self.stateOverTime = self.stateOverTime[0:self.counter+1, :]
        self.raycastData = self.raycastData[0:self.counter+1, :]
        self.controlInputData = self.controlInputData[0:self.counter+1]
        self.endTime = 1.0*self.counter/self.numTimesteps*self.endTime



    def initializeStatusBar(self):
        self.numTicks = 10
        self.nextTickComplete = 1.0 / float(self.numTicks)
        self.nextTickIdx = 1
        print "Simulation percentage complete: (", self.numTicks, " # is complete)"
    
    def printStatusBar(self):
        fractionDone = float(self.counter) / float(self.numTimesteps)
        if fractionDone > self.nextTickComplete:

            self.nextTickIdx += 1
            self.nextTickComplete += 1.0 / self.numTicks

            timeSoFar = time.time() - self.startSimTime 
            estimatedTimeLeft_sec = (1 - fractionDone) * timeSoFar / fractionDone
            estimatedTimeLeft_minutes = estimatedTimeLeft_sec / 60.0

            print "#" * self.nextTickIdx, "-" * (self.numTicks - self.nextTickIdx), "estimated", estimatedTimeLeft_minutes, "minutes left"


    def setCollisionFreeInitialState(self):
        tol = 5

        while True:
            
            x = 0.0
            y =   -5.0
            theta = 0 #+ np.random.uniform(0,2*np.pi,1)[0] * 0.01
            
            self.EgoCar.setCarState(x,y,theta)
            self.setRobotFrameState(x,y,theta)

            print "In loop"

            if not self.checkInCollision():
                break
                
        return x,y,theta


    def setRandomCollisionFreeInitialState(self):
        tol = 5

        while True:
            
            x = np.random.uniform(self.world.Xmin+tol, self.world.Xmax-tol, 1)[0]
            y = np.random.uniform(self.world.Ymin+tol, self.world.Ymax-tol, 1)[0]
            theta = np.random.uniform(0,2*np.pi,1)[0]
            
            self.EgoCar.setCarState(x,y,theta)
            self.setRobotFrameState(x,y,theta)

            if not self.checkInCollision():
                break

        return x,y,theta

    def setupPlayback(self):

        self.timer = TimerCallback(targetFps=30)
        self.timer.callback = self.tick

        playButtonFps = 1.0/self.dt
        print "playButtonFPS", playButtonFps
        self.playTimer = TimerCallback(targetFps=playButtonFps)
        self.playTimer.callback = self.playTimerCallback
        self.sliderMovedByPlayTimer = False

        panel = QtGui.QWidget()
        l = QtGui.QHBoxLayout(panel)

        playButton = QtGui.QPushButton('Play/Pause')
        playButton.connect('clicked()', self.onPlayButton)

        slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        slider.connect('valueChanged(int)', self.onSliderChanged)
        self.sliderMax = self.numTimesteps
        slider.setMaximum(self.sliderMax)
        self.slider = slider

        l.addWidget(playButton)
        l.addWidget(slider)

        w = QtGui.QWidget()
        l = QtGui.QVBoxLayout(w)
        l.addWidget(self.view)
        self.view.orientationMarkerWidget().Off()
        l.addWidget(panel)
        w.showMaximized()

        self.frame.connectFrameModified(self.updateDrawIntersection)
        # self.frame.connectFrameModified(self.updateDrawPolyApprox)
        self.updateDrawIntersection(self.frame)
        # self.updateDrawPolyApprox(self.frame)

        camera = self.view.camera()
        camera_control_panel = cameracontrolpanel.CameraControlPanel(self.view)


        if self.mode=='Obs' or self.mode=='2in1':
            robot_center = self.robot.getChildFrame().transform.TransformPoint([0,0,0])
            robot_camera = self.robot.getChildFrame().transform.TransformPoint([-20,0,10])
            camera.SetPosition(robot_camera)
            camera.SetFocalPoint(robot_center)
            panel = screengrabberpanel.ScreenGrabberPanel(self.view)
            panel.widget.show()
            robot = om.findObjectByName('AgentCar1') # or whatever you need to do to get the object
            # hacky way to be compatible with director update
            # TODO: clean up once director provides class with legit getTargetFrame() method
            # camera_control_panel.trackerManager.setTarget(TargetFrameConverter(robot))

            camera_control_panel.trackerManager.target = robot;
            camera_control_panel.trackerManager.targetFrame = robot.getChildFrame();
            camera_control_panel.trackerManager.callbackId = camera_control_panel.trackerManager.targetFrame.connectFrameModified(camera_control_panel.trackerManager.onTargetFrameModified)
            camera_control_panel.trackerManager.initTracker()
            camera_control_panel.trackerManager.setTrackerMode('Smooth Follow')
        elif self.mode=='Goal':
            applogic.resetCamera(viewDirection=[0.2,0,-1])
            # camera_control_panel.trackerManager.setTarget(self.world)
            # camera_control_panel.trackerManager.setTrackerMode('Position')


            # applogic.resetCamera(viewDirection=[0.2,0,-1])
            # self.view.showMaximized()
            self.view.raise_()
        else:
            print 'view camera mode error'
        

        # cameracontrolpanel.CameraControlPanel(self.view).widget.show()

        # robot = om.findObjectByName('robot') # or whatever you need to do to get the object
        # cameracontrolpanel.CameraControlPanel.trackerManager.setTarget(robot)
        # cameracontrolpanel.CameraControlPanel.trackerManager.setTrackerMode('Position & Orientation')
      
        camera_control_panel.widget.show()

        elapsed = time.time() - self.startSimTime
        simRate = self.counter/elapsed
        print "Total run time", elapsed
        print "Ticks (Hz)", simRate
        print "Number of steps taken", self.counter
        self.app.start()

    def run(self, launchApp=True):
        self.counter = 1
        self.runBatchSimulation()
        if launchApp:
            self.setupPlayback()


    def updateDrawPolyApprox(self, frame):
        distances = self.Sensor.raycastAll(frame)
        polyCoefficients = self.SensorApproximator.polyFitConstrainedLP(distances)
    
        d = DebugData()
        
        x = self.SensorApproximator.approxThetaVector
        y = x * 0.0
        for index,val in enumerate(y):
            y[index] = self.horner(x[index],polyCoefficients)
        
        origin = np.array(frame.transform.GetPosition())
        origin[2] = -0.001

        for i in xrange(self.SensorApproximator.numApproxPoints):
            if y[i] > 0:
                ray = self.SensorApproximator.approxRays[:,i]
                rayTransformed = np.array(frame.transform.TransformNormal(ray))
                intersection = origin + rayTransformed * y[i]
                intersection[2] = -0.001
                d.addLine(origin, intersection, color=[0,0.1,1])

        vis.updatePolyData(d.getPolyData(), 'polyApprox', colorByName='RGB255')

    def horner(self, x, weights):
        coefficients = weights[::-1]
        result = 0
        for i in coefficients:
            result = result * x + i
        return result
        

    def updateDrawIntersection(self, frame):

        origin = np.array(frame.transform.GetPosition())
        #print "origin is now at", origin
        d = DebugData()

        sliderIdx = self.slider.value

        controllerType = self.getControllerTypeFromCounter(sliderIdx)
        colorMaxRange = self.colorMap[controllerType]

        for i in xrange(self.Sensor.numRays):
            ray = self.Sensor.rays[:,i]
            rayTransformed = np.array(frame.transform.TransformNormal(ray))
            #print "rayTransformed is", rayTransformed
            intersection = self.Sensor.raycast(self.locator, origin, origin + rayTransformed*self.Sensor.rayLength)

            if intersection is not None:
                d.addLine(origin, intersection, color=[1,0,0])
            else:
                d.addLine(origin, origin+rayTransformed*self.Sensor.rayLength, color=colorMaxRange)
        vis.updatePolyData(d.getPolyData(), 'rays', colorByName='RGB255')

        #camera = self.view.camera()
        #camera.SetFocalPoint(frame.transform.GetPosition())
        #camera.SetPosition(frame.transform.TransformPoint((-30,0,10)))


    def getControllerTypeFromCounter(self, counter):
        name = self.controllerTypeOrder[0]
        for controllerType in self.controllerTypeOrder[1:]:
            if counter >= self.idxDict[controllerType]:
                name = controllerType
        return name


    def setRobotFrameState(self, x, y, theta):
        t = vtk.vtkTransform()
        t.Translate(x,y,0.0)
        t.RotateZ(np.degrees(theta))
        self.robot.getChildFrame().copyFrame(t)

    # returns true if we are in collision
    def checkInCollision(self, raycastDistance=None):
        if raycastDistance is None:
            self.setRobotFrameState(self.EgoCar.state[0],self.EgoCar.state[1],self.EgoCar.state[2])
            raycastDistance = self.Sensor.raycastAll(self.frame)
        if np.min(raycastDistance) < self.collisionThreshold:
            return True
        else:
            return False

    def tick(self):
        #print timer.elapsed
        #simulate(t.elapsed)
        x = np.sin(time.time())
        y = np.cos(time.time())
        self.setRobotFrameState(x,y,0.0)
        if (time.time() - self.playTime) > self.endTime:
            self.playTimer.stop()

    # just increment the slider, stop the timer if we get to the end
    def playTimerCallback(self):
        self.sliderMovedByPlayTimer = True
        currentIdx = self.slider.value
        nextIdx = currentIdx + 1
        self.slider.setSliderPosition(nextIdx)
        if currentIdx >= self.sliderMax:
            print "reached end of tape, stopping playTime"
            self.playTimer.stop()

    def onSliderChanged(self, value):
        if not self.sliderMovedByPlayTimer:
            self.playTimer.stop()
        numSteps = len(self.stateOverTime)
        idx = int(np.floor(numSteps*(1.0*value/self.sliderMax)))
        idx = min(idx, numSteps-1)
        x,y,theta = self.stateOverTime[idx]
        phi,alpha,u = self.angleOverTime[idx]
        ray=self.raycastData[idx]
        if not self.sliderMovedByPlayTimer:
            # print 'ray cast'
            print ray
        self.setRobotFrameState(x,y,theta)
        self.sliderMovedByPlayTimer = False

    def onPlayButton(self):
        if self.playTimer.isActive():
            self.onPauseButton()
            return
        print 'play'
        self.playTimer.start()
        self.playTime = time.time()

    def onPauseButton(self):
        print 'pause'
        self.playTimer.stop()

    def saveToFile(self, filename):
        # should also save the run data if it is available, i.e. stateOverTime, rewardOverTime
        filename = 'data/' + filename + ".txt"
        my_shelf = shelve.open(filename,'n')
        my_shelf['options'] = self.options
        my_shelf['simulationData'] = self.simulationData
        my_shelf['stateOverTime'] = self.stateOverTime
        my_shelf['raycastData'] = self.raycastData
        my_shelf['controlInputData'] = self.controlInputData
        my_shelf['numTimesteps'] = self.numTimesteps
        my_shelf['idxDict'] = self.idxDict
        my_shelf['counter'] = self.counter
        my_shelf.close()
    @staticmethod
    def loadFromFile(filename):
        filename = 'data/' + filename + ".out"
        sim = Simulator(autoInitialize=False, verbose=False)

        my_shelf = shelve.open(filename)
        sim.options = my_shelf['options']
        sim.initialize()

        sim.simulationData = my_shelf['simulationData']
        sim.stateOverTime = np.array(my_shelf['stateOverTime'])
        sim.raycastData = np.array( my_shelf['raycastData'])
        sim.controlInputData = np.array(my_shelf['controlInputData'])
        sim.numTimesteps = my_shelf['numTimesteps']
        sim.idxDict = my_shelf['idxDict']
        sim.counter = my_shelf['counter']
        my_shelf.close()
        return sim

if __name__ == "__main__":
    # main(sys.argv[1:])
    parser = argparse.ArgumentParser(description='interpret simulation parameters')
    parser.add_argument('--percentObsDensity', type=float, nargs=1, default=[10])
    parser.add_argument('--endTime', type=int, nargs=1, default=[40])
    parser.add_argument('--nonRandomWorld', action='store_true', default=False)
    parser.add_argument('--circleRadius', type=float, nargs=1, default=0.7)
    parser.add_argument('--worldScale', type=float, nargs=1, default=1.0)
    parser.add_argument('--numCars', type=int, nargs=1, default=1)

    argNamespace = parser.parse_args()
    percentObsDensity = argNamespace.percentObsDensity[0]
    endTime = argNamespace.endTime[0]
    nonRandomWorld = argNamespace.nonRandomWorld
    circleRadius = argNamespace.circleRadius[0]
    worldScale = argNamespace.worldScale[0]
    numCars = argNamespace.numCars[0]
    
    sim = Simulator(percentObsDensity=percentObsDensity, endTime=endTime,
                    nonRandomWorld=nonRandomWorld, circleRadius=circleRadius, worldScale=worldScale,
                    numCars=numCars)
    sim.run()


