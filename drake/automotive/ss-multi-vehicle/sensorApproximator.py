import numpy as np
import scipy.optimize as opt
# from linear_regression import LinearRegression
# import cvxopt
import math
from director.debugVis import DebugData
import director.visualization as vis

class SensorApproximatorObj(object):

    def __init__(self, numRays, circleRadius):
        self.N = 1
        self.numRays = numRays
        self.circleRadius = circleRadius

    def initializeThetaVector(self,thetaVector):
        self.thetaVector = thetaVector

    def initializeApproxThetaVector(self, angleMin, angleMax):
        self.numApproxPoints = 200
        self.approxThetaVector = np.linspace(angleMin, angleMax, self.numApproxPoints)
        self.approxRays = np.zeros((3,self.numApproxPoints))
        self.approxRays[0,:] = np.cos(self.approxThetaVector)
        self.approxRays[1,:] = -np.sin(self.approxThetaVector)
        

    def polyFitConstrainedLP(self, distances):
        self.laserDepths = np.array(distances) - np.ones((np.shape(distances)))*self.circleRadius # decrease each sensor by the circle radius (i.e., inflate all obstacles)
        #self.laserDepths[0] = 0.1
        #self.laserDepths[-1] = 0.1
        self.setUpOptimization()
        self.constrainedLP()
        return self.polyCoefficientsLP


    
    