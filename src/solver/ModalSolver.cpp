#include "fem/solver/ModalSolver.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace fem {
namespace {

std::vector<Eigen::Index> freeDofs(const Model& model, const DofManager& dofs) {
  std::vector<Eigen::Index> result;
  result.reserve(dofs.size());

  for (const NodeId node_id : dofs.nodeOrder()) {
    const Node& node = model.node(node_id);
    for (std::size_t offset = 0; offset < kDofsPerFrameNode; ++offset) {
      const Dof dof = dofFromOffset(offset);
      if (!node.isFixed(dof)) {
        result.push_back(static_cast<Eigen::Index>(dofs.equation(node_id, dof)));
      }
    }
  }
  return result;
}

Eigen::MatrixXd selectMatrix(const Eigen::MatrixXd& matrix,
                             const std::vector<Eigen::Index>& indices) {
  const Eigen::Index n = static_cast<Eigen::Index>(indices.size());
  Eigen::MatrixXd reduced(n, n);
  for (Eigen::Index i = 0; i < n; ++i) {
    for (Eigen::Index j = 0; j < n; ++j) {
      reduced(i, j) = matrix(indices[static_cast<std::size_t>(i)],
                             indices[static_cast<std::size_t>(j)]);
    }
  }
  return reduced;
}

}  // namespace

ModalResult ModalSolver::solve(const Model& model, std::size_t requested_modes) const {
  AssembledSystem system = model.assemble();
  const auto free = freeDofs(model, system.dofs);
  if (free.empty()) {
    throw std::runtime_error("Model has no free degrees of freedom");
  }

  const Eigen::MatrixXd kff = selectMatrix(system.stiffness, free);
  const Eigen::MatrixXd mff = selectMatrix(system.mass, free);

  Eigen::LLT<Eigen::MatrixXd> mass_factor(mff);
  if (mass_factor.info() != Eigen::Success) {
    throw std::runtime_error(
        "Reduced mass matrix is not positive definite; check material density and constraints");
  }

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigen(kff, mff);
  if (eigen.info() != Eigen::Success) {
    throw std::runtime_error("Generalized eigenvalue solution failed");
  }

  std::vector<Eigen::Index> positive_modes;
  const double scale = std::max(1.0, eigen.eigenvalues().cwiseAbs().maxCoeff());
  const double tolerance = 1.0e-10 * scale;
  for (Eigen::Index i = 0; i < eigen.eigenvalues().size(); ++i) {
    if (eigen.eigenvalues()[i] > tolerance) {
      positive_modes.push_back(i);
    }
  }

  if (positive_modes.empty()) {
    throw std::runtime_error(
        "No positive structural modes found; check supports and stiffness");
  }

  std::size_t count = positive_modes.size();
  if (requested_modes != 0U) {
    count = std::min(count, requested_modes);
  }

  constexpr double kTwoPi = 6.28318530717958647692;
  std::vector<double> omega;
  std::vector<double> hz;
  omega.reserve(count);
  hz.reserve(count);

  Eigen::MatrixXd full_modes =
      Eigen::MatrixXd::Zero(system.stiffness.rows(), static_cast<Eigen::Index>(count));

  for (std::size_t mode = 0; mode < count; ++mode) {
    const Eigen::Index source = positive_modes[mode];
    const double w = std::sqrt(eigen.eigenvalues()[source]);
    omega.push_back(w);
    hz.push_back(w / kTwoPi);

    for (std::size_t i = 0; i < free.size(); ++i) {
      full_modes(free[i], static_cast<Eigen::Index>(mode)) =
          eigen.eigenvectors()(static_cast<Eigen::Index>(i), source);
    }
  }

  return {std::move(omega), std::move(hz), std::move(full_modes), std::move(system.dofs)};
}

}  // namespace fem
