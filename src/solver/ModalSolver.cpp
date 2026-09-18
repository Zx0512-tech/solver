#include "fem/solver/ModalSolver.hpp"
#include "fem/solver/SparseDofReducer.hpp"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SymShiftInvert.h>
#include <Spectra/SymGEigsShiftSolver.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace fem {
namespace {

std::vector<Eigen::Index> freeDofs(
    const Model& model,
    const DofManager& dofs) {
  std::vector<Eigen::Index> result;
  result.reserve(dofs.size());

  for (const NodeId node_id : dofs.nodeOrder()) {
    const Node& node = model.node(node_id);
    for (std::size_t offset = 0;
         offset < kDofsPerFrameNode;
         ++offset) {
      const Dof dof = dofFromOffset(offset);
      if (!node.isFixed(dof)) {
        result.push_back(
            static_cast<Eigen::Index>(
                dofs.equation(node_id, dof)));
      }
    }
  }
  return result;
}

double generalizedScale(
    const Eigen::SparseMatrix<double>& stiffness,
    const Eigen::SparseMatrix<double>& mass) {
  double scale = 1.0;
  for (Eigen::Index i = 0; i < stiffness.rows(); ++i) {
    const double mii = std::abs(mass.coeff(i, i));
    if (mii > 0.0) {
      scale = std::max(
          scale,
          std::abs(stiffness.coeff(i, i)) / mii);
    }
  }
  return scale;
}

std::size_t automaticModeCount(std::size_t free_dofs) {
  if (free_dofs <= 1U) {
    return free_dofs;
  }
  return std::min<std::size_t>(6U, free_dofs - 1U);
}

ModalResult solveSingleDof(
    Eigen::SparseMatrix<double> stiffness,
    Eigen::SparseMatrix<double> mass,
    AssembledSystem system,
    const std::vector<Eigen::Index>& free,
    std::size_t requested_modes) {
  if (requested_modes > 1U) {
    throw std::invalid_argument(
        "A single-free-DOF model has only one structural mode");
  }

  const double m = mass.coeff(0, 0);
  if (!(m > 0.0)) {
    throw std::runtime_error(
        "Reduced mass matrix is not positive definite; "
        "check material density and constraints");
  }

  const double k = stiffness.coeff(0, 0);
  const double lambda = k / m;
  if (!(lambda > 0.0) || !std::isfinite(lambda)) {
    throw std::runtime_error(
        "No positive structural modes found; "
        "check supports and stiffness");
  }

  constexpr double kTwoPi = 6.28318530717958647692;
  const double omega = std::sqrt(lambda);

  Eigen::MatrixXd full_modes =
      Eigen::MatrixXd::Zero(system.stiffness.rows(), 1);
  // Mass-normalize the one-dimensional eigenvector.
  full_modes(free.front(), 0) = 1.0 / std::sqrt(m);

  return {
      {omega},
      {omega / kTwoPi},
      std::move(full_modes),
      std::move(system.dofs),
      ModalSolverBackend::DirectSingleDof,
      0U,
      0U};
}

}  // namespace

ModalResult ModalSolver::solve(
    const Model& model,
    std::size_t requested_modes) const {
  AssembledSystem system = model.assemble();
  const auto free = freeDofs(model, system.dofs);
  if (free.empty()) {
    throw std::runtime_error(
        "Model has no free degrees of freedom");
  }

  Eigen::SparseMatrix<double> kff =
      SparseDofReducer::matrix(system.stiffness, free);
  Eigen::SparseMatrix<double> mff =
      SparseDofReducer::matrix(system.mass, free);
  kff.makeCompressed();
  mff.makeCompressed();

  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>
      mass_factor;
  mass_factor.compute(mff);
  if (mass_factor.info() != Eigen::Success) {
    throw std::runtime_error(
        "Reduced mass matrix is not positive definite; "
        "check material density and constraints");
  }

  const std::size_t n =
      static_cast<std::size_t>(kff.rows());
  std::size_t target_modes = requested_modes;
  if (target_modes == 0U) {
    target_modes = automaticModeCount(n);
  }

  if (n == 1U) {
    return solveSingleDof(
        std::move(kff),
        std::move(mff),
        std::move(system),
        free,
        target_modes);
  }

  if (target_modes == 0U || target_modes >= n) {
    throw std::invalid_argument(
        "Sparse Lanczos modal solve requires requested_modes "
        "to satisfy 1 <= requested_modes < free DOF count");
  }

  // A small negative shift avoids singular factorization for free-free
  // structures while preserving the ordering of nonnegative eigenvalues:
  // the eigenvalues closest to the negative shift are the rigid-body modes
  // followed by the lowest positive structural modes.
  const double scale = generalizedScale(kff, mff);
  const double shift =
      -std::max(1.0, scale) * 1.0e-8;
  const double positive_tolerance =
      std::max(1.0e-12, scale * 1.0e-10);

  constexpr std::size_t kRigidBodyAllowance = 6U;
  const std::size_t candidate_modes =
      std::min(
          n - 1U,
          target_modes + kRigidBodyAllowance);

  const std::size_t minimum_subspace =
      candidate_modes + 1U;
  const std::size_t preferred_subspace =
      std::max(
          minimum_subspace,
          2U * candidate_modes + 1U);
  const std::size_t ncv =
      std::min(n, preferred_subspace);

  using ShiftInvertOp =
      Spectra::SymShiftInvert<
          double,
          Eigen::Sparse,
          Eigen::Sparse>;
  using MassOp = Spectra::SparseSymMatProd<double>;
  using Solver =
      Spectra::SymGEigsShiftSolver<
          ShiftInvertOp,
          MassOp,
          Spectra::GEigsMode::ShiftInvert>;

  ShiftInvertOp shift_invert(kff, mff);
  MassOp mass_op(mff);

  Solver eigen(
      shift_invert,
      mass_op,
      static_cast<Eigen::Index>(candidate_modes),
      static_cast<Eigen::Index>(ncv),
      shift);

  eigen.init();
  const Eigen::Index converged =
      eigen.compute(
          Spectra::SortRule::LargestMagn,
          1000,
          1.0e-10,
          Spectra::SortRule::SmallestAlge);

  if (eigen.info() != Spectra::CompInfo::Successful ||
      converged <= 0) {
    throw std::runtime_error(
        "Sparse Lanczos generalized eigenvalue solution failed");
  }

  const Eigen::VectorXd eigenvalues =
      eigen.eigenvalues();
  const Eigen::MatrixXd eigenvectors =
      eigen.eigenvectors();

  std::vector<Eigen::Index> positive_modes;
  positive_modes.reserve(target_modes);
  for (Eigen::Index i = 0;
       i < eigenvalues.size();
       ++i) {
    const double lambda = eigenvalues[i];
    if (std::isfinite(lambda) &&
        lambda > positive_tolerance) {
      positive_modes.push_back(i);
      if (positive_modes.size() == target_modes) {
        break;
      }
    }
  }

  if (positive_modes.size() < target_modes) {
    throw std::runtime_error(
        "Sparse Lanczos did not recover the requested number "
        "of positive structural modes; check supports, "
        "mechanisms, or request fewer modes");
  }

  constexpr double kTwoPi = 6.28318530717958647692;
  std::vector<double> omega;
  std::vector<double> hz;
  omega.reserve(target_modes);
  hz.reserve(target_modes);

  Eigen::MatrixXd full_modes =
      Eigen::MatrixXd::Zero(
          system.stiffness.rows(),
          static_cast<Eigen::Index>(target_modes));

  for (std::size_t mode = 0;
       mode < target_modes;
       ++mode) {
    const Eigen::Index source = positive_modes[mode];
    const double lambda = eigenvalues[source];
    const double w = std::sqrt(lambda);

    Eigen::VectorXd reduced_mode =
        eigenvectors.col(source);
    const double mass_norm_squared =
        reduced_mode.dot(mff * reduced_mode);
    if (!(mass_norm_squared > 0.0) ||
        !std::isfinite(mass_norm_squared)) {
      throw std::runtime_error(
          "Lanczos mode has invalid mass norm");
    }
    reduced_mode /= std::sqrt(mass_norm_squared);

    omega.push_back(w);
    hz.push_back(w / kTwoPi);

    for (std::size_t i = 0;
         i < free.size();
         ++i) {
      full_modes(
          free[i],
          static_cast<Eigen::Index>(mode)) =
          reduced_mode[
              static_cast<Eigen::Index>(i)];
    }
  }

  return {
      std::move(omega),
      std::move(hz),
      std::move(full_modes),
      std::move(system.dofs),
      ModalSolverBackend::SparseLanczosShiftInvert,
      static_cast<std::size_t>(eigen.num_iterations()),
      static_cast<std::size_t>(eigen.num_operations())};
}

}  // namespace fem
