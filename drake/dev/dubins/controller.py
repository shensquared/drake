import numpy as np
import scipy.integrate as integrate
import director.objectmodel as om
import math
import shelve
import time
import os
import matplotlib.pyplot as plt

class ControllerObj(object):

    def __init__(self, sensor, sensor_approximator, mode, goalX, goalY, u_max=4, epsilonRand=0.4):
        self.Sensor = sensor
        self.SensorApproximator = sensor_approximator
        self.SensorApproximator.initializeThetaVector(self.Sensor.angleGrid)
        self.SensorApproximator.initializeApproxThetaVector(self.Sensor.angleMin, self.Sensor.angleMax)
        self.numRays = self.Sensor.numRays
        self.actionSet = np.array([u_max,0,-u_max])
        self.epsilonRand = epsilonRand
        self.actionSetIdx = np.arange(0,np.size(self.actionSet))
        self.u_max = u_max
        self.slackParam = 0.1
        self.k = 5
        self.kTurn = 50000000
        self.goalX=goalX
        self.goalY=goalY
        self.goalTheta=0
        self.mode=mode
        self.history=[]
        self.angleState = np.array([0,0,0])
        self.time=time.time()

    def addingCar (self,Car):
        self.Car=Car

    def initializeVelocity(self,velocity):
        self.velocity = velocity
        
    def computeControlInput(self, state, t, frame, raycastDistance=None, randomize=False):
        self.Car.state = state
        if raycastDistance is None:
            self.distances = self.Sensor.raycastAll(frame)
        else:
            self.distances = raycastDistance
        if self.mode=='Goal':
            u = self.mothController()
        elif self.mode=='Obs':
            u = self.naiveController()
        elif self.mode=='2in1':
            u = self.twoin1Controller() 
        elif self.mode=='Manual':
            u = self.manualController() 
        else:
            print 'controller-mode setup error'

        return u

    def manualController(self, steering):
        return steering

    def mothController(self, arrived=False):
        distances=self.distances
        state=self.Car.state
        v=self.velocity
        x=state[0]
        y=state[1]
        theta=state[2]

        # if theta>2*math.pi:
        #     print 'theta is in degree'

        u_max=self.u_max
        goalY=self.goalY
        goalX=self.goalX

        r=np.sqrt((goalY-y)**2+(goalX-x)**2)
        if not arrived:
            # degree from rad
            phi=(math.atan2((goalY-y),(goalX-x)))
            # still the range would be -360 to 360
            alpha=(phi-theta)
            sine=np.sin(alpha)
            cosine=np.cos(alpha)
            if cosine<0:
                u=u_max
            else:
                u=sine
            if u>u_max:
                print 'u_max reached'
                u=u_max
            elif u<-u_max:
                print '-u_max reached'
                u=-u_max
        angleState = np.array([phi, alpha, u])

            # self.phi.append(phi)
            # self.history.append([phi,theta,alpha,sine,cosine,r])
            # path = '/Users/shenshen/ofcar/ofcar_debug'
            # filename='debug'+str(np.int(self.time))+'.txt'
            # filename = os.path.join(path, filename)
            # f = open(filename, 'w+')
            # f.write(str(self.history))
            # f.close
        return u


    def naiveController(self):
        distances=self.distances
        state=self.Car.state
        x=state[0]
        y=state[1]
        theta=state[2]
        degree=180/math.pi

        goalY=self.goalY
        goalX=self.goalX
  
        fod=np.diff(distances)
        fod=np.absolute(fod)

        total_index=len(self.distances)
        mid_index = (total_index-1)/2 
        maxindex=np.argmax(fod)
            
        if np.max(distances)<=15:
            if maxindex<mid_index and distances[maxindex]<distances[maxindex+1]: 
                u=-(28*(mid_index-1-maxindex))/total_index
            elif maxindex<mid_index and distances[maxindex]>=distances[maxindex+1]: 
                u=-(28*(mid_index-maxindex))/total_index
            elif maxindex>mid_index and distances[maxindex]<distances[maxindex+1]: 
                u=(28*(maxindex-mid_index+1))/total_index
            elif maxindex>mid_index and distances[maxindex]>=distances[maxindex+1]:
                u=(28*(maxindex-mid_index))/total_index
            else:
                u=0
               
        else:
            for i in range(0,mid_index+1):
                if distances[mid_index-i]>15:
                    u=-(28*i)/total_index
                    break
                elif distances[mid_index+i]>15:
                    u=(28*i)/total_index
                    break

        if u > self.u_max:
            u = self.u_max
        if u < -self.u_max:
            u = -self.u_max
        angleState=self.angleState
        return -u

    def twoin1Controller(self, arrived=False):
        distances=self.distances
        state=self.Car.state
        x=state[0]
        y=state[1]
        theta=state[2]
        degree=180/math.pi


        goalY=self.goalY
        goalX=self.goalX
    
        fod=np.diff(distances)
        fod=np.absolute(fod)

        total_index=len(self.distances)
        mid_index = (total_index-1)/2 
        maxindex=np.argmax(fod)
            

        if np.min(distances)>=15:
            u, angels, actionIdx= self.mothController()
        else:
            u, angels, actionIdx= self.naiveController()
        

        if u > self.u_max:
            u = self.u_max
        if u < -self.u_max:
            u = -self.u_max
        return u 

