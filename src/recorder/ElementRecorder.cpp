#include "fem/recorder/ElementRecorder.hpp"

#include <stdexcept>

namespace fem {

ElementResponse ElementRecorder::record(
    const Model& model,
    ElementId element_id,
    const StaticResult& result) const {
  return model.elementResponse(element_id, result.displacement, 0.0);
}

ElementResponseHistory ElementRecorder::record(
    const Model& model,
    ElementId element_id,
    const NewmarkResult& result) const {
  if (result.displacement.cols() != static_cast<Eigen::Index>(result.time.size())) {
    throw std::invalid_argument(
        "Newmark result time and displacement history sizes do not match");
  }

  const Eigen::Index step_count = result.displacement.cols();
  if (step_count == 0) {
    throw std::invalid_argument("Cannot record an empty displacement history");
  }

  const ElementResponse first =
      model.elementResponse(element_id, result.displacement.col(0), result.time[0]);
  const Eigen::Index ndof = first.local_end_force.size();

  ElementResponseHistory history;
  history.element_id = element_id;
  history.time = result.time;
  history.local_end_force = Eigen::MatrixXd::Zero(ndof, step_count);
  history.global_end_force = Eigen::MatrixXd::Zero(ndof, step_count);
  history.local_end_force.col(0) = first.local_end_force;
  history.global_end_force.col(0) = first.global_end_force;

  for (Eigen::Index step = 1; step < step_count; ++step) {
    const ElementResponse response =
        model.elementResponse(
            element_id,
            result.displacement.col(step),
            result.time[static_cast<std::size_t>(step)]);
    history.local_end_force.col(step) = response.local_end_force;
    history.global_end_force.col(step) = response.global_end_force;
  }

  return history;
}

}  // namespace fem
