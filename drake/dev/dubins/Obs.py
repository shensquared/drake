from CarSimulator import Simulator
import argparse

parser = argparse.ArgumentParser(description='interpret simulation parameters')
parser.add_argument('--mode', type=str, nargs= 1, default='Obs')

argNamespace = parser.parse_args()
mode = argNamespace.mode[0]

sim = Simulator(mode=mode ,autoInitialize=False, verbose=False)

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
