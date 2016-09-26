#include "drake/automotive/car_simulation.h"
#include "drake/systems/LCMSystem.h"
#include "drake/systems/LinearSystem.h"
#include "drake/systems/pd_control_system.h"
#include "drake/systems/plants/BotVisualizer.h"
#include "drake/systems/plants/parser_model_instance_id_table.h"
#include "drake/systems/plants/RigidBodySystem.h"
#include "drake/util/drakeAppUtil.h"
#include "lcmtypes/drake/lcmt_driving_command_t.hpp"
#include <lcm/lcm.h>
#include <bot_lcmgl_client/lcmgl.h>

using Eigen::VectorXd;
// bot_lcmgl_t *lcmgl;


namespace drake {
namespace automotive {
namespace {

int do_main(int argc, const char* argv[]) {
  // Initializes the communication layer.
  std::shared_ptr<lcm::LCM> lcm = std::make_shared<lcm::LCM>();

  // Instantiates a duration variable that will be set by the call to
  // CreateRigidBodySystem() below.
  double duration = std::numeric_limits<double>::infinity();

  drake::parsers::ModelInstanceIdTable model_instances;

  double penetration_stiffness = atof(argv[argc-6]);
  double penetration_damping = atof(argv[argc-5]);
  double friction_coefficient = atof(argv[argc-4]); 

  // Initializes the rigid body system.
  auto rigid_body_sys = CreateRigidBodySystem(argc, argv, &duration,
      &model_instances, penetration_stiffness, penetration_damping, friction_coefficient);

  // Initializes and cascades all of the other systems.

  double steering_kp = atof(argv[argc-3]);
  double steering_kd = atof(argv[argc-2]);
  double throttle_k = atof(argv[argc-1]); 


  auto vehicle_sys = CreateVehicleSystem(rigid_body_sys, steering_kp, steering_kd, throttle_k);

  auto const& tree = rigid_body_sys->getRigidBodyTree();
  auto visualizer =
      std::make_shared<BotVisualizer<RigidBodySystem::StateVector>>(lcm, tree);

  auto sys = cascade(vehicle_sys, visualizer);

  // Initializes the simulation options.
  SimulationOptions options =
      GetCarSimulationDefaultOptions();

  // Defines the start time of the simulation.
  const double kStartTime = 0;

  // Starts the simulation.
  drake::runLCM(sys, lcm, kStartTime, duration,
                GetInitialState(*(rigid_body_sys.get())),
                options);

  // visulize contact forces and friction forces

  // lcm_t *lcm;
  // lcm = lcm_create(nullptr);
  // if (!lcm) return 1;

  // lcmgl = bot_lcmgl_init(lcm, "zmp-based CoM estimate");
  // drake_lcmt_zmp_com_observer_state_subscribe(lcm, "ZMP_COM_OBSERVER_STATE",
  //                                             &handleMessage, nullptr);

  // while (true) lcm_handle(lcm);

  // bot_lcmgl_destroy(lcmgl);
  // lcm_destroy(lcm);

  return 0;


}

}  // namespace
}  // namespace automotive
}  // namespace drake

int main(int argc, const char* argv[]) {
  return drake::automotive::do_main(argc, argv);
}

// static void handleMessage(const lcm_recv_buf_t *rbuf, const char *channel,
//                           const drake_lcmt_zmp_com_observer_state *msg,
//                           void *user) {


//   bot_lcmgl_color3f(lcmgl, 1.0, 0.0, 0.0);          // red
//   double v0 = 
//   double v1 = 
//   bot_lcmgl_vertex2d(lcmgl, v0, v1)
//   bot_lcmgl_switch_buffer(lcmgl);
// }


// drake_lcmt_zmp_com_observer_state_subscription_t* drake_lcmt_zmp_com_observer_state_subscribe (lcm_t *lcm,
//                     const char *channel,
//                     drake_lcmt_zmp_com_observer_state_handler_t f, void *userdata)
// {
//     drake_lcmt_zmp_com_observer_state_subscription_t *n = (drake_lcmt_zmp_com_observer_state_subscription_t*)
//                        malloc(sizeof(drake_lcmt_zmp_com_observer_state_subscription_t));
//     n->user_handler = f;
//     n->userdata = userdata;
//     n->lc_h = lcm_subscribe (lcm, channel,
//                                  drake_lcmt_zmp_com_observer_state_handler_stub, n);
//     if (n->lc_h == NULL) {
//         fprintf (stderr,"couldn't reg drake_lcmt_zmp_com_observer_state LCM handler!\n");
//         free (n);
//         return NULL;
//     }
//     return n;
// }
