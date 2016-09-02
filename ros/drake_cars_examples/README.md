Drake / ROS Car Simulation README
=================================

This README provides details about and instructions on how to run Drake's car
demonstrations. It assumes you have [installed Drake within a ROS Catkin
workspace](http://drake.mit.edu/from_source_ros.html).

Demo 1: Single Car in State Garage
==================================

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
================================

To run this demo, open six termainals. In the first terminal, execute the
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
