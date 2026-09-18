#pragma once

#include "fem/core/Model.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace fem {

enum class ModalSolverBackend {
  SparseLanczosShiftInvert,
  DirectSingleDof,
};

struct ModalResult {
  std::vector<double> angular_frequencies;
  std::vector<double> frequencies_hz;
  Eigen::MatrixXd mode_shapes;
  DofManager dofs;
  ModalSolverBackend backend{ModalSolverBackend::SparseLanczosShiftInvert};
  std::size_t iterations{0};
  std::size_t operations{0};

  double modeValue(std::size_t mode_index, NodeId node_id, Dof dof) const {
    return mode_shapes(
        static_cast<Eigen::Index>(dofs.equation(node_id, dof)),
        static_cast<Eigen::Index>(mode_index));
  }
};

class ModalSolver {
 public:
  // requested_modes == 0 selects an automatic low-mode count (up to 6).
  // For multi-DOF systems, the sparse Lanczos backend computes a partial
  // spectrum, so an explicit requested mode count must be smaller than the
  // number of free DOFs.
  ModalResult solve(const Model& model, std::size_t requested_modes = 0) const;
};

}  // namespace fem
