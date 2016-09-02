#include <cmath>

#include "ros/ros.h"
#include "ros/console.h"

#include "rosgraph_msgs/Clock.h"

#include "drake/examples/Cars/car_simulation.h"
#include "drake/examples/Cars/gen/driving_command.h"
#include "drake/ros/simulation_abort_function.h"
#include "drake/ros/systems/ros_clock_publisher.h"
#include "drake/ros/systems/ros_sensor_publisher_joint_state.h"
#include "drake/ros/systems/ros_sensor_publisher_lidar.h"
#include "drake/ros/systems/ros_sensor_publisher_odometry.h"
#include "drake/ros/systems/ros_tf_publisher.h"
#include "drake/ros/systems/ros_vehicle_system.h"
#include "drake/systems/LCMSystem.h"
#include "drake/systems/LinearSystem.h"
#include "drake/systems/pd_control_system.h"
#include "drake/systems/plants/BotVisualizer.h"
#include "drake/systems/plants/RigidBodySystem.h"
#include "drake/util/drakeAppUtil.h"

using drake::BotVisualizer;
using drake::SimulationOptions;
using drake::cascade;

using Eigen::VectorXd;

namespace drake {
namespace ros {
namespace cars {
namespace {

using drake::examples::cars::CreateRigidBodySystem;
using drake::examples::cars::CreateVehicleSystem;
using drake::examples::cars::GetCarSimulationDefaultOptions;
using drake::examples::cars::ParseDuration;

using drake::ros::systems::RosTfPublisher;
using drake::ros::systems::RosClockPublisher;
using drake::ros::systems::RosSensorPublisherJointState;
using drake::ros::systems::RosSensorPublisherLidar;
using drake::ros::systems::RosSensorPublisherOdometry;

using drake::ros::systems::run_ros_vehicle_sim;

/**
 * This implements the main method of the single car simulation. The vehicle
 * resides within the Stata garage.
 */
int DoMain(int argc, const char* argv[]) {
  ::ros::init(argc, const_cast<char**>(argv), "single_car_in_stata_garage");

  // Instantiates a ROS node handle. For more information, see:
  // http://wiki.ros.org/roscpp/Overview/NodeHandles.
  ::ros::NodeHandle node_handle;

  // Sets the log level to be INFO.
  if(::ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME,
      ::ros::console::levels::Info) ) {
    ::ros::console::notifyLoggerLevelsChanged();
  }

  // Initializes the communication layer.
  std::shared_ptr<lcm::LCM> lcm = std::make_shared<lcm::LCM>();

  // Instantiates a duration variable that will be set by the call to
  // CreateRigidBodySystem() below.
  double duration = std::numeric_limits<double>::infinity();

  // Instantiates a data structure that maps model instance names to their model
  // instance IDs.
  drake::parsers::ModelInstanceIdTable model_instances;

  // Initializes the rigid body system.
  auto rigid_body_sys = CreateRigidBodySystem(argc, argv, &duration,
      &model_instances);

  // Sets the desired contact penetration stiffness and damping in the
  // RigidBodySystem.
  rigid_body_sys->penetration_stiffness =
      GetROSParameter<double>(node_handle, "penetration_stiffness");

  rigid_body_sys->penetration_damping =
      GetROSParameter<double>(node_handle, "penetration_damping");

  auto const& tree = rigid_body_sys->getRigidBodyTree();

  // Instantiates a map that converts model instance IDs to model instance
  // names.
  std::map<int, std::string> model_instance_name_table;
  // TODO(liang.fok): Once #3088 is resolved, include the model instance ID and
  // name of the world in model_instance_name_table.
  model_instance_name_table[model_instances["prius_1"]] = "prius";
  model_instance_name_table[model_instances["P1"]] = "stata_garage";

  std::map<int, std::string> model_instance_name_table_odometry;
  model_instance_name_table_odometry[model_instances["prius_1"]] = "prius";

  // Obtains the gains to be used by the steering and throttle controllers.
  double steering_kp = GetROSParameter<double>(node_handle, "steering_kp");
  double steering_kd = GetROSParameter<double>(node_handle, "steering_kd");
  double throttle_k = GetROSParameter<double>(node_handle, "throttle_k");

  // Initializes and cascades all of the other systems.

  // Wraps the RigidBodySystem within a PD control system that adds PD
  // controllers for each actuator within the RigidBodySystem. It then cascades
  // the PD control system block behind a gain block and returns the resulting
  // cascade.
  auto vehicle_sys = CreateVehicleSystem(rigid_body_sys, steering_kp,
      steering_kd, throttle_k);

  auto visualizer =
      std::make_shared<BotVisualizer<RigidBodySystem::StateVector>>(lcm, tree);

  auto lidar_publisher = std::make_shared<
      RosSensorPublisherLidar<RigidBodySystem::StateVector>>(
      rigid_body_sys, model_instance_name_table);

  auto odometry_publisher = std::make_shared<
      RosSensorPublisherOdometry<RigidBodySystem::StateVector>>(
      rigid_body_sys, model_instance_name_table_odometry);

  auto tf_publisher = std::make_shared<
      RosTfPublisher<RigidBodySystem::StateVector>>(tree,
          model_instance_name_table);

  auto joint_state_publisher = std::make_shared<
      RosSensorPublisherJointState<RigidBodySystem::StateVector>>(
      rigid_body_sys, model_instance_name_table);

  auto clock_publisher =
      std::make_shared<RosClockPublisher<RigidBodySystem::StateVector>>();

  auto sys =
      cascade(
        cascade(
          cascade(
            cascade(
              cascade(
                cascade(
                  vehicle_sys, visualizer),
                lidar_publisher),
              odometry_publisher),
            tf_publisher),
          joint_state_publisher),
        clock_publisher);

  // Instantiates a ROS topic publisher for publishing clock information. For
  // more information, see: http://wiki.ros.org/Clock.
  // ::ros::Publisher clock_publisher =
  //     node_handle.advertise<rosgraph_msgs::Clock>("/clock", 1);

  // Initializes the simulation options.
  SimulationOptions options = GetCarSimulationDefaultOptions();
  AddAbortFunction(&options);

  options.initial_step_size =
      GetROSParameter<double>(node_handle, "initial_step_size");

  ROS_INFO_STREAM("Using:" << std::endl
      << " - penetration_stiffness = " << rigid_body_sys->penetration_stiffness
      << std::endl
      << " - penetration_damping = " << rigid_body_sys->penetration_damping
      << std::endl
      << "  - steering_kp = " << steering_kp
      << std::endl
      << "  - steering_kd = " << steering_kd
      << std::endl
      << "  - throttle_k = " << throttle_k
      << std::endl
      << "  - initial_step_size = " << options.initial_step_size);

  // Obtains a valid zero configuration for the vehicle.
  VectorXd x0 = VectorXd::Zero(rigid_body_sys->getNumStates());
  x0.head(tree->number_of_positions()) = tree->getZeroConfiguration();

  // Defines the start time of the simulation.
  const double kStartTime = 0;

  // Starts the simulation.
  run_ros_vehicle_sim(sys, kStartTime, duration, x0, options);

  return 0;
}

}  // namespace
}  // namespace cars
}  // namespace ros
}  // namespace drake

int main(int argc, const char* argv[]) {
  return drake::ros::cars::DoMain(argc, argv);
}
