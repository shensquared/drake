Drake / ROS Car Simulation README
=================================

This README provides details about and instructions on how to run Drake's car
demonstrations. It assumes you have [installed Drake within a ROS Catkin
workspace](http://drake.mit.edu/from_source_ros.html).

Demo 1: Single Car in MIT Stata Garage
======================================

To run this demo, open two termainals. In the first terminal, execute the
following command to launch Drake and RViz. The simulation runs in
Drake while RViz serves as the visualizer.

```
$ roslaunch drake_cars_examples single_car_in_stata_garage.launch
```

In the second terminal, execute the following command to launch an application
that allows you to issue driving commands to the simulated vehicle:

```
$ rosrun ackermann_drive_teleop ackermann_drive_keyop.py 1.0 0.7 /drake/ackermann_cmd
```

You can now type arrow keys in your second terminal to issue driving commands to
the simulated vehicle.

Demo 2: Five Cars on a Plane
============================

To run this demo, open six terminals. In the first terminal, execute the
following command to launch Drake and RViz. The simulation runs in
Drake while RViz serves as the visualizer.

```
$ roslaunch drake_cars_examples multi_cars_on_plane.launch
```

Next, execute the following commands, one in each of the remaining five
terminals.

```
$ rosrun ackermann_drive_teleop ackermann_drive_keyop.py 1.0 0.7 /drake/Prius_1/ackermann_cmd
$ rosrun ackermann_drive_teleop ackermann_drive_keyop.py 1.0 0.7 /drake/Prius_2/ackermann_cmd
$ rosrun ackermann_drive_teleop ackermann_drive_keyop.py 1.0 0.7 /drake/Prius_3/ackermann_cmd
$ rosrun ackermann_drive_teleop ackermann_drive_keyop.py 1.0 0.7 /drake/Prius_4/ackermann_cmd
$ rosrun ackermann_drive_teleop ackermann_drive_keyop.py 1.0 0.7 /drake/Prius_5/ackermann_cmd
```

You can now issue driving commands to the vehicles using the terminals running
`ackermann_drive_keyop.py`.

Demo 3: Maximum Acceleration of Prius
=====================================

A ROS node that makes the vehicle in the single car demo execute a maximum
acceleration trajectory is available in
`nodes/max_accel_2016_toyota_prius.py`.

To run the demo, first start Demo 1 as described above. You can omit the second
command that runs `ackermann_drive_keyop.py`. Instead, execute the following
command, which starts a node that publishes a maximum acceleration trajectory:

```
$ roslaunch drake_cars_examples max_accel_2016_toyota_prius.launch
```

You can plot the reference and actual speed of the vehicle by executing:

```
$ rqt_plot /drake/ackermann_cmd/drive/speed
```

Tuning Tips
-----------

There are several parameters that impact the stability of the vehicle and
simulation. These parameters are loaded onto the ROS parameter server and can
be changed in
`drake_cars_examples/launch/single_car_in_stata_garage.launch`. Below are
descriptions of these parameters.

Contacts are modeled using virtual springs. The stiffness and damping gains of
these springs can be set by modifying parameters "penetration_stiffness" and
"penetration_damping".

The steering and throttle of the vehicle are PD controlled. The gains of these
controllers can be tweaked using parameters "steering_kp", "steering_kd", and
"throttle_k".

The initial time step using by the simulation's "solver" is set by parameter
"initial_step_size".