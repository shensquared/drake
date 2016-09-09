#!/usr/bin/env python

'''
A ROS node that publishes the maximum acceleration curve for a 2016 Toyota
Prius. The acceleration curve was obtained from the following graph:

http://hiwaay.net/~bzwilson/prius/2016_metrics_010.jpg

The polynomial curve fitting the graph data is not correct. See comments
embedded in the code for more details on how the maximum acceleration trajectory
was determined.
'''

__author__ = 'Chien-Liang Fok'
__license__ = 'MIT'
__maintainer__ = 'Chien-Liang Fok'
__email__ = 'drake-users@mit.edu'

import math
import rospy
from ackermann_msgs.msg import AckermannDriveStamped

def ComputeReferenceSpeed(elapsed_time):
  # The following equation was taken from the following graph but does not match
  # the data in the graph probably because of insufficient numbers of decimal
  # places in the coefficient.
  #
  # http://hiwaay.net/~bzwilson/prius/2016_metrics_010.jpg
  #
  # speed_mph = 3e-6 * math.pow(elapsed_time, 6) - \
  #             0.0002 * math.pow(elapsed_time, 5) + \
  #             0.0044 * math.pow(elapsed_time, 4) - \
  #             0.0338 * math.pow(elapsed_time, 3) - \
  #             0.2874 * math.pow(elapsed_time, 2) + \
  #             9.4805 * elapsed_time
  #
  # Since the equation above did not work, I manually created a new graph in
  # LibreOffice Calc based on visual inspection of the points in the original
  # graph linked to above and used curve fitting to derive the equation below.
  # See drake-distro/ros/drake_cars_examples/docs/max_accel.ods.
  #
  speed_mph = -5.64212882544153e-5 * math.pow(elapsed_time, 5) \
              + 0.0037713708 * math.pow(elapsed_time, 4) \
              - 0.0870244015 * math.pow(elapsed_time, 3) \
              + 0.6660764631 * math.pow(elapsed_time, 2) \
              + 4.5202731988 * elapsed_time

  # Converts the speed from MPH to m/s.
  speed = speed_mph * 0.44704
  return speed

if __name__ == '__main__':
  rospy.init_node('max_accel_curve_2016_toyota_prius', anonymous=True)

  # Defines the total duration in seconds of the speed trajectory.
  kTrajectoryDuration = 23.5

  # Instantiates a ROS topic publisher that publishes the AckermannDriveStamped
  # command message.
  command_pub = rospy.Publisher("/drake/ackermann_cmd", AckermannDriveStamped,
      queue_size=1)

  # Loop waiting for the first time to be received
  start_time = rospy.get_time()
  while start_time == 0:
    rospy.sleep(0.01)
    start_time = rospy.get_time()

  elapsed_time = rospy.get_time() - start_time;
  print "Initial elapsed Time = {0}".format(elapsed_time)

  while elapsed_time < kTrajectoryDuration:
    speed = ComputeReferenceSpeed(elapsed_time)
    print "Elapsed Time = {0}, speed = {1}".format(elapsed_time, speed)

    # Saves the reference speed in an AckermannDriveStamped message.
    ackermann_cmd_msg = AckermannDriveStamped()
    ackermann_cmd_msg.header.stamp = rospy.Time.now()
    ackermann_cmd_msg.drive.speed = speed
    ackermann_cmd_msg.drive.steering_angle = 0

    # Publishes the message.
    command_pub.publish(ackermann_cmd_msg)
    rospy.sleep(0.01)

    elapsed_time = rospy.get_time() - start_time;
    print "Elapsed Time = {0}, speed = {1}".format(elapsed_time, speed)