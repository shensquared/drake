#include "drake/automotive/gen/driving_command_translator.h"

// This file is generated by a script.  Do not edit!
// See drake/automotive/lcm_vector_gen.py.

#include <stdexcept>

#include "drake/common/drake_assert.h"

namespace drake {
namespace automotive {

std::unique_ptr<systems::BasicVector<double>>
DrivingCommandTranslator::AllocateOutputVector() const {
  return std::make_unique<DrivingCommand<double>>();
}

void DrivingCommandTranslator::TranslateVectorBaseToLcm(
    const systems::VectorBase<double>& vector_base,
    std::vector<uint8_t>* lcm_message_bytes) const {
  const auto* const vector =
      dynamic_cast<const DrivingCommand<double>*>(&vector_base);
  DRAKE_DEMAND(vector != nullptr);
  drake::lcmt_driving_command_t message;
  message.steering_angle = vector->steering_angle();
  message.throttle = vector->throttle();
  message.brake = vector->brake();
  const int lcm_message_length = message.getEncodedSize();
  lcm_message_bytes->resize(lcm_message_length);
  message.encode(lcm_message_bytes->data(), 0, lcm_message_length);
}

void DrivingCommandTranslator::TranslateLcmToVectorBase(
    const void* lcm_message_bytes, int lcm_message_length,
    systems::VectorBase<double>* vector_base) const {
  DRAKE_DEMAND(vector_base != nullptr);
  auto* const my_vector = dynamic_cast<DrivingCommand<double>*>(vector_base);
  DRAKE_DEMAND(my_vector != nullptr);

  drake::lcmt_driving_command_t message;
  int status = message.decode(lcm_message_bytes, 0, lcm_message_length);
  if (status < 0) {
    throw std::runtime_error("Failed to decode LCM message driving_command.");
  }
  my_vector->set_steering_angle(message.steering_angle);
  my_vector->set_throttle(message.throttle);
  my_vector->set_brake(message.brake);
}

}  // namespace automotive
}  // namespace drake
