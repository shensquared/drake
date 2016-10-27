Multiple Vehicle Simulation
===========================

This README file provides instructions on how to run Drake's car simulations.

 
Start the Drake Visualizer
--------------------------

The Drake Visualizer displays the current state of the simulation. It is a
separate process that communicates with the Drake simulation process via the
[Lightweight Communications and Marshalling (LCM)](https://lcm-proj.github.io/)
middleware.

To run the Drake Visualizer, open a terminal and execute the following commands:

```
$ cd [drake-distro]/drake/automotive/ss-multi-vehicle
$ ../../../build/install/bin/directorPython Obs.py
``` 
There are several python package dependencies; most of them are actually non-essential for the demo to run and I'll need to clean things up. For the moment, it might be the easist to follow the error and pip install the dependencies. 

Testing Different Controllers
--------------------------
Add the controller as a method in the controller.py file, and update the controller mode in the simulation launch python script

Use Ctrl-C in your terminal to stop and close the demo.
