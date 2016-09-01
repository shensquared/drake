#pragma once

#include <cmath>

#include "ros/ros.h"
#include "rosgraph_msgs/Clock.h"

#include "drake/systems/simulation_options.h"

namespace drake {
namespace ros {

/**
 * Adds a custom stop function that (1) checks whether the simulation should
 * abort based on a call to ros::ok(), and (2) publishes a clock message so
 * other systems within ROS can be time synchronized with simulation time.
 */
void AddAbortFunction(drake::SimulationOptions* options,
    ::ros::Publisher* clock_publisher) {
  options->should_stop = [clock_publisher](double sim_time) {
    // Computes the whole-second portion and the fractional-second portion of
    // sim_time.
    double whole_part, fractional_part;
    fractional_part = modf(sim_time, &whole_part);

    // Saves the time in the clock message.
    rosgraph_msgs::Clock clock_msg;
    clock_msg.clock.sec = static_cast<int>(whole_part);
    clock_msg.clock.nsec = static_cast<int>(fractional_part * 1e9);
    clock_publisher->publish(clock_msg);

    // Publishes the clock message.
    return !::ros::ok();
  };
}

}  // namespace ros
}  // namespace drake
