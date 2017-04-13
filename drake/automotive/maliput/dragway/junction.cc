#include "drake/automotive/maliput/dragway/junction.h"

#include "drake/automotive/maliput/dragway/road_geometry.h"
#include "drake/automotive/maliput/dragway/segment.h"

namespace drake {
namespace maliput {
namespace dragway {

Junction::Junction(RoadGeometry* road_geometry,
	int num_segments,
    int num_lanes,
    double length,
    double lane_width,
    double shoulder_width)
  : id_({"Dragway Junction"}),
    road_geometry_(road_geometry),
    num_segments_(num_segments){

  DRAKE_DEMAND(road_geometry != nullptr);

    for (int i = 0; i < num_segments; ++i) {
    auto segment = std::make_unique<Segment>(this, 
    	api::SegmentId({"Dragway_Segment_" + std::to_string(i)}),
      i,
    	num_lanes, length, lane_width, shoulder_width);
    segments_.push_back(move(segment));
  };
}

const api::Segment* Junction::do_segment(int index) const {
    DRAKE_DEMAND(index < num_segments());
    return segments_.at(index).get();
}


const api::RoadGeometry* Junction::do_road_geometry() const {
  return road_geometry_;
}

}  // namespace dragway
}  // namespace maliput
}  // namespace drake
