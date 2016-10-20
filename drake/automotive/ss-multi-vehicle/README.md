Multiple Vehicle Simulation
===========================

This README file provides instructions on how to run Drake's car simulations.

The instructions are written for Ubuntu Linux and OS X users. Windows users will
need to adjust the instructions slightly. See the notes at the end of this
section.
 
Start the Drake Visualizer
--------------------------

The Drake Visualizer displays the current state of the simulation. It is a
separate process that communicates with the Drake simulation process via the
[Lightweight Communications and Marshalling (LCM)](https://lcm-proj.github.io/)
middleware.

To run the Drake Visualizer, open a terminal and execute the following commands:

```
$ cd [drake-distro]/drake/automotive/ss-multi-vehicle
$ ../../../build/install/bin/directorPython runMultiDemo.py
``` 


Testing Different Controller
--------------------------
Change the controller.py file, and update the controller mode in the simulation
start-up python script


Use Ctrl-C in your terminal to stop and close the demo.
