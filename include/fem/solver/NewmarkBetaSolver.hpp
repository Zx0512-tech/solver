#pragma once

#include "fem/core/Model.hpp"
#include "fem/dynamics/RayleighDamping.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <functional>
#include <vector>

namespace fem {

struct NodalTimeLoad {
  NodeId node_id;
  Dof dof;
  std::function<double(double)> value;
};

struct NewmarkSettings {
  double time_step{0.01};
  std::size_t step_count{1};
  double beta{0.25};
  double gamma{0.5};
};

struct DynamicInitialState {
  Eigen::VectorXd displacement;
  Eigen::VectorXd velocity;
};

struct NewmarkResult {
  std::vector<double> time;
  Eigen::MatrixXd displacement;
  Eigen::MatrixXd velocity;
  Eigen::MatrixXd acceleration;
  DofManager dofs;

  double displacementAt(std::size_t step, NodeId node_id, Dof dof) const {
    return displacement(
        static_cast<Eigen::Index>(dofs.equation(node_id, dof)),
        static_cast<Eigen::Index>(step));
  }

  double velocityAt(std::size_t step, NodeId node_id, Dof dof) const {
    return velocity(
        static_cast<Eigen::Index>(dofs.equation(node_id, dof)),
        static_cast<Eigen::Index>(step));
  }

  double accelerationAt(std::size_t step, NodeId node_id, Dof dof) const {
    return acceleration(
        static_cast<Eigen::Index>(dofs.equation(node_id, dof)),
        static_cast<Eigen::Index>(step));
  }
};

class NewmarkBetaSolver {
 public:
  NewmarkResult solve(const Model& model,
                      const NewmarkSettings& settings,
                      const RayleighDamping& damping = RayleighDamping{},
                      const std::vector<NodalTimeLoad>& time_loads = {},
                      const DynamicInitialState& initial = {}) const;
};

}  // namespace fem
