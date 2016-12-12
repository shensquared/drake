from CarSimulator import Simulator
import argparse

parser = argparse.ArgumentParser(description='interpret simulation parameters')
parser.add_argument('--mode', type=str, nargs= 1, default='Obs')
parser.add_argument('--numCars', type=int, nargs=1, default=2)


argNamespace = parser.parse_args()
mode = argNamespace.mode[0]
numCars = argNamespace.numCars[0]

sim = Simulator(mode=mode ,autoInitialize=False, verbose=False, numCars = numCars)

sim.Sarsa_numInnerBins = 4
sim.Sarsa_numOuterBins = 4
sim.Sensor_rayLength = 10


sim.randomSeed = 7
sim.randomizeControl       = True
sim.percentObsDensity      = 4
sim.nonRandomWorld         = True
sim.circleRadius           = 2.5
sim.worldScale             = 1


sim.initialize()
sim.defaultControllerTime  = 100
sim.run()
