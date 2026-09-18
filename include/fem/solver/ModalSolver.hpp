#pragma once

#include "fem/core/Model.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace fem {

struct ModalResult {
  std::vector<double> angular_frequencies;
  std::vector<double> frequencies_hz;
  Eigen::MatrixXd mode_shapes;
  DofManager dofs;

  double modeValue(std::size_t mode_index, NodeId node_id, Dof dof) const {
    return mode_shapes(
        static_cast<Eigen::Index>(dofs.equation(node_id, dof)),
        static_cast<Eigen::Index>(mode_index));
  }
};

class ModalSolver {
 public:
  ModalResult solve(const Model& model, std::size_t requested_modes = 0) const;
};

}  // namespace fem
