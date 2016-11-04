Multiple Vehicle Simulation
===========================

This README file provides instructions on how to run Drake's car simulations.

 
Start the Drake Visualizer
--------------------------

The Drake Visualizer displays the current state of the simulation. It is a
separate process that communicates with the Drake simulation process via the
[Lightweight Communications and Marshalling (LCM)](https://lcm-proj.github.io/)
middleware.

To run a demo in which the dubins car uses a naive controller to avoid obstacle, open a terminal and execute the following commands:

```
$ cd [drake-distro]/drake/automotive/ss-multi-vehicle
$ ../../../build/install/bin/directorPython Obs.py
``` 
There may be unnecessary python package dependencies; I will clean those up later. For the moment, it might be the easiest to follow the error msg, if any, and pip install the dependencies. 

Testing Different Controllers
--------------------------
Add your desired controller as a method in the controller.py class, and update the controller mode in the simulation launch python script.

Use Ctrl-C in your terminal to stop and close the demo.
