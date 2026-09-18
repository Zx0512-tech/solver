#pragma once

#include "fem/core/Types.hpp"
#include "fem/recorder/ResultHistory.hpp"
#include "fem/solver/LinearStaticSolver.hpp"
#include "fem/solver/NewmarkBetaSolver.hpp"

#include <Eigen/Core>

#include <vector>

namespace fem {

enum class NodeResponseQuantity {
  Displacement,
  Velocity,
  Acceleration,
};

struct NodeResponse {
  NodeId node_id{};
  Eigen::Matrix<double, 6, 1> displacement{
      Eigen::Matrix<double, 6, 1>::Zero()};
  Eigen::Matrix<double, 6, 1> reaction{
      Eigen::Matrix<double, 6, 1>::Zero()};
};

struct NodeResponseHistory {
  NodeId node_id{};
  std::vector<double> time;
  Eigen::MatrixXd displacement;
  Eigen::MatrixXd velocity;
  Eigen::MatrixXd acceleration;

  ScalarHistory series(
      NodeResponseQuantity quantity,
      Dof dof) const;
};

class NodeRecorder {
 public:
  NodeResponse record(
      NodeId node_id,
      const StaticResult& result) const;

  NodeResponseHistory record(
      NodeId node_id,
      const NewmarkResult& result) const;
};

}  // namespace fem
