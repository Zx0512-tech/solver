#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/solver/ModalSolver.hpp"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void fail(const std::string& message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(EXIT_FAILURE);
}

void expectNear(double actual,
                double expected,
                double relative_tolerance,
                const std::string& message) {
  const double scale = std::max(1.0, std::abs(expected));
  if (std::abs(actual - expected) > relative_tolerance * scale) {
    std::cerr << "FAIL: " << message
              << " actual=" << actual
              << " expected=" << expected << '\n';
    std::exit(EXIT_FAILURE);
  }
}

void fixAll(fem::Node& node) {
  for (std::size_t i = 0; i < fem::kDofsPerFrameNode; ++i) {
    node.fix(fem::dofFromOffset(i));
  }
}

void constrainPlanarYBending(fem::Node& node) {
  node.fix(fem::Dof::UX);
  node.fix(fem::Dof::UZ);
  node.fix(fem::Dof::RX);
  node.fix(fem::Dof::RY);
}

std::vector<Eigen::Index> freeDofs(
    const fem::Model& model,
    const fem::DofManager& dofs) {
  std::vector<Eigen::Index> free;
  free.reserve(dofs.size());
  for (const fem::NodeId node_id : dofs.nodeOrder()) {
    const fem::Node& node = model.node(node_id);
    for (std::size_t offset = 0;
         offset < fem::kDofsPerFrameNode;
         ++offset) {
      const fem::Dof dof = fem::dofFromOffset(offset);
      if (!node.isFixed(dof)) {
        free.push_back(static_cast<Eigen::Index>(
            dofs.equation(node_id, dof)));
      }
    }
  }
  return free;
}

Eigen::MatrixXd denseReduced(
    const Eigen::SparseMatrix<double>& matrix,
    const std::vector<Eigen::Index>& free) {
  const Eigen::Index n =
      static_cast<Eigen::Index>(free.size());
  Eigen::MatrixXd result(n, n);
  for (Eigen::Index i = 0; i < n; ++i) {
    for (Eigen::Index j = 0; j < n; ++j) {
      result(i, j) = matrix.coeff(
          free[static_cast<std::size_t>(i)],
          free[static_cast<std::size_t>(j)]);
    }
  }
  return result;
}

std::vector<double> densePositiveEigenvalues(
    const fem::Model& model) {
  const auto system = model.assemble();
  const auto free = freeDofs(model, system.dofs);
  const Eigen::MatrixXd kff =
      denseReduced(system.stiffness, free);
  const Eigen::MatrixXd mff =
      denseReduced(system.mass, free);

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>
      eigen(kff, mff);
  if (eigen.info() != Eigen::Success) {
    fail("dense reference generalized eigensolve failed");
  }

  std::vector<double> result;
  const double scale =
      std::max(1.0,
               eigen.eigenvalues().cwiseAbs().maxCoeff());
  const double tolerance = 1.0e-9 * scale;
  for (Eigen::Index i = 0;
       i < eigen.eigenvalues().size();
       ++i) {
    if (eigen.eigenvalues()[i] > tolerance) {
      result.push_back(eigen.eigenvalues()[i]);
    }
  }
  return result;
}

fem::Model makePlanarCantilever(int element_count) {
  constexpr double length = 6.0;
  constexpr double e = 210.0e9;
  constexpr double rho = 7850.0;
  constexpr double area = 0.02;
  constexpr double iy = 5.0e-5;
  constexpr double iz = 8.0e-5;

  fem::Model model;
  for (int i = 0; i <= element_count; ++i) {
    auto& node = model.addNode(
        static_cast<fem::NodeId>(i + 1),
        {length * static_cast<double>(i) /
             static_cast<double>(element_count),
         0.0,
         0.0});
    if (i == 0) {
      fixAll(node);
    } else {
      constrainPlanarYBending(node);
    }
  }

  const fem::LinearElasticMaterial material(e, 0.3, rho);
  const fem::BeamSection section(area, iy, iz, 7.0e-5);
  for (int i = 0; i < element_count; ++i) {
    model.addElement<fem::Beam3D>(
        static_cast<fem::ElementId>(i + 1),
        static_cast<fem::NodeId>(i + 1),
        static_cast<fem::NodeId>(i + 2),
        material,
        section);
  }
  return model;
}

void testLanczosMatchesDenseReference() {
  fem::Model model = makePlanarCantilever(20);
  constexpr std::size_t requested = 4;

  const auto result =
      fem::ModalSolver{}.solve(model, requested);
  const auto reference =
      densePositiveEigenvalues(model);

  if (result.backend !=
      fem::ModalSolverBackend::SparseLanczosShiftInvert) {
    fail("modal solver must report sparse Lanczos backend");
  }
  if (result.angular_frequencies.size() != requested) {
    fail("sparse Lanczos returned wrong mode count");
  }
  if (result.iterations == 0U ||
      result.operations == 0U) {
    fail("Lanczos diagnostics must report iterations/operations");
  }

  for (std::size_t i = 0; i < requested; ++i) {
    const double expected_omega =
        std::sqrt(reference.at(i));
    expectNear(
        result.angular_frequencies[i],
        expected_omega,
        1.0e-8,
        "Lanczos frequency vs dense reference");
  }

  const auto system = model.assemble();
  for (std::size_t mode = 0;
       mode < requested;
       ++mode) {
    const Eigen::VectorXd phi =
        result.mode_shapes.col(
            static_cast<Eigen::Index>(mode));
    const double lambda =
        result.angular_frequencies[mode] *
        result.angular_frequencies[mode];
    const Eigen::VectorXd residual =
        system.stiffness * phi -
        lambda * (system.mass * phi);
    const double scale =
        std::max(
            1.0,
            (system.stiffness * phi).norm());
    if (residual.norm() / scale > 1.0e-8) {
      fail("Lanczos mode residual is too large");
    }
  }
}

void testFreeFreeRigidModesAreSkipped() {
  constexpr double length = 3.0;
  fem::Model model;
  model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {length, 0.0, 0.0});
  model.addElement<fem::Beam3D>(
      1,
      1,
      2,
      fem::LinearElasticMaterial(210.0e9, 0.3, 7850.0),
      fem::BeamSection(0.015, 4.0e-5, 6.0e-5, 8.0e-5));

  const auto reference =
      densePositiveEigenvalues(model);
  const auto result =
      fem::ModalSolver{}.solve(model, 2);

  if (result.angular_frequencies.size() != 2U) {
    fail("free-free sparse modal solve must return two positive modes");
  }

  for (std::size_t i = 0; i < 2U; ++i) {
    expectNear(
        result.angular_frequencies[i],
        std::sqrt(reference.at(i)),
        1.0e-7,
        "free-free positive Lanczos mode");
  }
}

}  // namespace

int main() {
  testLanczosMatchesDenseReference();
  testFreeFreeRigidModesAreSkipped();
  std::cout << "All sparse-modal tests passed.\n";
  return EXIT_SUCCESS;
}
