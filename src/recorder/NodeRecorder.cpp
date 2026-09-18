#include "fem/recorder/NodeRecorder.hpp"

#include <stdexcept>

namespace fem {
namespace {

const Eigen::MatrixXd& quantityMatrix(
    const NodeResponseHistory& history,
    NodeResponseQuantity quantity) {
  switch (quantity) {
    case NodeResponseQuantity::Displacement:
      return history.displacement;
    case NodeResponseQuantity::Velocity:
      return history.velocity;
    case NodeResponseQuantity::Acceleration:
      return history.acceleration;
  }
  throw std::logic_error("Unknown NodeResponseQuantity");
}

void validateHistoryMatrix(
    const Eigen::MatrixXd& matrix,
    std::size_t time_size,
    const char* name) {
  if (matrix.rows() != static_cast<Eigen::Index>(kDofsPerFrameNode) ||
      matrix.cols() != static_cast<Eigen::Index>(time_size)) {
    throw std::invalid_argument(
        std::string("Node response ") + name +
        " history must be 6 x time.size()");
  }
}

}  // namespace

ScalarHistory NodeResponseHistory::series(
    NodeResponseQuantity quantity,
    Dof dof) const {
  const Eigen::MatrixXd& matrix = quantityMatrix(*this, quantity);
  validateHistoryMatrix(matrix, time.size(), "component");

  ScalarHistory result;
  result.time = time;
  result.value.resize(time.size());

  const Eigen::Index row =
      static_cast<Eigen::Index>(dofOffset(dof));
  for (std::size_t step = 0; step < time.size(); ++step) {
    result.value[step] =
        matrix(row, static_cast<Eigen::Index>(step));
  }
  return result;
}

NodeResponse NodeRecorder::record(
    NodeId node_id,
    const StaticResult& result) const {
  NodeResponse response;
  response.node_id = node_id;

  for (std::size_t offset = 0; offset < kDofsPerFrameNode; ++offset) {
    const Dof dof = dofFromOffset(offset);
    const Eigen::Index row = static_cast<Eigen::Index>(offset);
    response.displacement[row] =
        result.displacementAt(node_id, dof);
    response.reaction[row] =
        result.reactionAt(node_id, dof);
  }
  return response;
}

NodeResponseHistory NodeRecorder::record(
    NodeId node_id,
    const NewmarkResult& result) const {
  const std::size_t step_count = result.time.size();
  const Eigen::Index global_dof_count =
      static_cast<Eigen::Index>(result.dofs.size());

  for (const auto* matrix :
       {&result.displacement, &result.velocity, &result.acceleration}) {
    if (matrix->rows() != global_dof_count ||
        matrix->cols() != static_cast<Eigen::Index>(step_count)) {
      throw std::invalid_argument(
          "Newmark result matrices must match DOF count and time size");
    }
  }

  NodeResponseHistory history;
  history.node_id = node_id;
  history.time = result.time;
  history.displacement =
      Eigen::MatrixXd::Zero(
          static_cast<Eigen::Index>(kDofsPerFrameNode),
          static_cast<Eigen::Index>(step_count));
  history.velocity = history.displacement;
  history.acceleration = history.displacement;

  for (std::size_t offset = 0; offset < kDofsPerFrameNode; ++offset) {
    const Dof dof = dofFromOffset(offset);
    const Eigen::Index local_row =
        static_cast<Eigen::Index>(offset);
    const Eigen::Index global_row =
        static_cast<Eigen::Index>(
            result.dofs.equation(node_id, dof));

    history.displacement.row(local_row) =
        result.displacement.row(global_row);
    history.velocity.row(local_row) =
        result.velocity.row(global_row);
    history.acceleration.row(local_row) =
        result.acceleration.row(global_row);
  }

  return history;
}

}  // namespace fem
