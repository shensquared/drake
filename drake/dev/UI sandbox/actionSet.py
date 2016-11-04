import director.vtkAll as vtk
import director.visualization as vis

import numpy as np
import director.objectmodel as om

from director.debugVis import DebugData


class ActionSetObj(object):

    def __init__(self):
        self.a_max = 25 # this is approximately for a 70 degree pitch angle

        self.num_x_bins = 5
        self.a_x = np.linspace(-self.a_max, self.a_max, self.num_x_bins)

        self.num_y_bins = 5
        self.a_y = np.linspace(-self.a_max, self.a_max, self.num_y_bins)

        self.t_f = 0.500 # 500 ms simulate forward time

        self.numPointsToDraw = 10
        self.t_vector = np.linspace(0,self.t_f,self.numPointsToDraw)
        self.t_vector_squared = 1.0*self.t_vector
        for index, value in enumerate(self.t_vector_squared):
            self.t_vector_squared[index] = value**2


    def computeFinalPositions(self, v_x_initial, v_y_initial):
        self.p_x_final = 1.0/2.0 * self.a_x * self.t_f**2 + np.ones(self.num_x_bins) * v_x_initial *self.t_f
        self.p_y_final = 1.0/2.0 * self.a_y * self.t_f**2 + np.ones(self.num_y_bins) * v_y_initial *self.t_f


    def computeAllPositions(self, x_initial, y_initial, v_x_initial, v_y_initial):
        self.p_x_trajectories = 1.0/2.0 * np.outer(self.a_x, self.t_vector_squared) + np.outer(np.ones(self.num_x_bins) * v_x_initial, self.t_vector) + np.ones((self.numPointsToDraw,self.num_x_bins)).T*x_initial
        self.p_y_trajectories = 1.0/2.0 * np.outer(self.a_y, self.t_vector_squared) + np.outer(np.ones(self.num_y_bins) * v_y_initial, self.t_vector) + np.ones((self.numPointsToDraw,self.num_y_bins)).T*y_initial

        # print "SHAPES"
        # print np.shape(self.p_x_trajectories)
        # print np.shape(np.ones((self.numPointsToDraw,self.num_x_bins)).T)

    def drawActionSetFinal(self):
        #print "I am drawing the action set"

        d = DebugData()

        for x_index, x_value in enumerate(self.p_x_final):
            for y_index, y_value in enumerate(self.p_y_final):
        
                firstEndpt = (0.0,0.0,0.0)
                secondEndpt = (x_value,y_value,0.0)

                d.addLine(firstEndpt, secondEndpt, radius=0.02, color=[0.8,0,0.8])


        obj = vis.updatePolyData(d.getPolyData(), 'action_set', colorByName='RGB255')

    def drawActionSetFull(self):
        #print "I am drawing the action set"

        d = DebugData()

        for x_index, x_value in enumerate(self.a_x):
            for y_index, y_value in enumerate(self.a_y):

                for time_step_index in xrange(self.numPointsToDraw-1):
        
                    firstX = self.p_x_trajectories[x_index, time_step_index]
                    firstY = self.p_y_trajectories[y_index, time_step_index]

                    secondX = self.p_x_trajectories[x_index, time_step_index+1]
                    secondY = self.p_y_trajectories[y_index, time_step_index+1]

                    firstEndpt = (firstX,firstY,0.0)
                    secondEndpt = (secondX,secondY,0.0)

                    d.addLine(firstEndpt, secondEndpt, radius=0.02, color=[0.8,0,0.8])


        obj = vis.updatePolyData(d.getPolyData(), 'action_set', colorByName='RGB255')
        
    