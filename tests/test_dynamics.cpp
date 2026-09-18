#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/dynamics/RayleighDamping.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/solver/ModalSolver.hpp"
#include "fem/solver/NewmarkBetaSolver.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

void expectNear(double actual,
                double expected,
                double relative_tolerance,
                const std::string& message) {
  const double scale = std::max(1.0, std::abs(expected));
  if (std::abs(actual - expected) > relative_tolerance * scale) {
    std::cerr << "FAIL: " << message << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(EXIT_FAILURE);
  }
}

void fixAll(fem::Node& node) {
  for (std::size_t i = 0; i < fem::kDofsPerFrameNode; ++i) {
    node.fix(fem::dofFromOffset(i));
  }
}

void testConsistentMass() {
  const double rho = 7850.0;
  const double area = 0.02;
  const double iy = 3.0e-5;
  const double iz = 5.0e-5;
  const double length = 2.4;

  fem::Beam3D beam(
      1,
      1,
      2,
      fem::LinearElasticMaterial(210.0e9, 0.3, rho),
      fem::BeamSection(area, iy, iz, 7.0e-5));

  const auto m = beam.localMass(length);
  expectNear((m - m.transpose()).cwiseAbs().maxCoeff(), 0.0, 1.0e-12,
             "consistent mass matrix symmetry");

  Eigen::Matrix<double, 12, 1> rigid_x = Eigen::Matrix<double, 12, 1>::Zero();
  rigid_x[0] = 1.0;
  rigid_x[6] = 1.0;
  expectNear((rigid_x.transpose() * m * rigid_x)(0, 0), rho * area * length, 1.0e-12,
             "axial rigid-translation mass");

  Eigen::Matrix<double, 12, 1> rigid_y = Eigen::Matrix<double, 12, 1>::Zero();
  rigid_y[1] = 1.0;
  rigid_y[7] = 1.0;
  expectNear((rigid_y.transpose() * m * rigid_y)(0, 0), rho * area * length, 1.0e-12,
             "transverse-y rigid-translation mass");

  Eigen::Matrix<double, 12, 1> rigid_z = Eigen::Matrix<double, 12, 1>::Zero();
  rigid_z[2] = 1.0;
  rigid_z[8] = 1.0;
  expectNear((rigid_z.transpose() * m * rigid_z)(0, 0), rho * area * length, 1.0e-12,
             "transverse-z rigid-translation mass");

  Eigen::Matrix<double, 12, 1> rigid_rx = Eigen::Matrix<double, 12, 1>::Zero();
  rigid_rx[3] = 1.0;
  rigid_rx[9] = 1.0;
  expectNear((rigid_rx.transpose() * m * rigid_rx)(0, 0),
             rho * (iy + iz) * length,
             1.0e-12,
             "torsional rotary inertia");
}

void testRayleighTargets() {
  const double omega_1 = 10.0;
  const double omega_2 = 30.0;
  const double zeta_1 = 0.02;
  const double zeta_2 = 0.05;

  const auto damping =
      fem::RayleighDamping::fromModalTargets(omega_1, zeta_1, omega_2, zeta_2);

  expectNear(damping.dampingRatio(omega_1), zeta_1, 1.0e-12,
             "Rayleigh damping at first target");
  expectNear(damping.dampingRatio(omega_2), zeta_2, 1.0e-12,
             "Rayleigh damping at second target");
}

struct AxialModel {
  fem::Model model;
  double omega;
};

AxialModel makeAxialModel() {
  const double e = 210.0e9;
  const double rho = 7850.0;
  const double area = 0.01;
  const double length = 2.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);

  tip.fix(fem::Dof::UY);
  tip.fix(fem::Dof::UZ);
  tip.fix(fem::Dof::RX);
  tip.fix(fem::Dof::RY);
  tip.fix(fem::Dof::RZ);

  model.addElement<fem::Beam3D>(
      1,
      1,
      2,
      fem::LinearElasticMaterial(e, 0.3, rho),
      fem::BeamSection(area, 2.0e-5, 3.0e-5, 4.0e-5));

  const double omega = std::sqrt(3.0 * e / (rho * length * length));
  return {std::move(model), omega};
}

void testModalAxialMode() {
  auto fixture = makeAxialModel();
  const auto modal = fem::ModalSolver{}.solve(fixture.model, 1);

  if (modal.angular_frequencies.size() != 1U) {
    std::cerr << "FAIL: expected exactly one requested mode\n";
    std::exit(EXIT_FAILURE);
  }

  expectNear(modal.angular_frequencies[0], fixture.omega, 1.0e-10,
             "single-DOF axial modal frequency");
}

void testNewmarkUndampedFreeVibration() {
  auto fixture = makeAxialModel();
  const auto assembled = fixture.model.assemble();

  fem::DynamicInitialState initial;
  initial.displacement = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(assembled.dofs.size()));
  initial.velocity = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(assembled.dofs.size()));

  const double initial_u = 1.0e-3;
  initial.displacement[
      static_cast<Eigen::Index>(assembled.dofs.equation(2, fem::Dof::UX))] = initial_u;

  constexpr double two_pi = 6.28318530717958647692;
  const double period = two_pi / fixture.omega;

  fem::NewmarkSettings settings;
  settings.time_step = period / 200.0;
  settings.step_count = 200;
  settings.beta = 0.25;
  settings.gamma = 0.5;

  const auto result =
      fem::NewmarkBetaSolver{}.solve(fixture.model, settings, fem::RayleighDamping{}, {}, initial);

  expectNear(result.displacementAt(100, 2, fem::Dof::UX), -initial_u, 5.0e-6,
             "Newmark half-period displacement");
  expectNear(result.displacementAt(200, 2, fem::Dof::UX), initial_u, 5.0e-6,
             "Newmark full-period displacement");
  expectNear(result.velocityAt(200, 2, fem::Dof::UX), 0.0, 5.0e-5,
             "Newmark full-period velocity");
}

}  // namespace

int main() {
  testConsistentMass();
  testRayleighTargets();
  testModalAxialMode();
  testNewmarkUndampedFreeVibration();
  std::cout << "All dynamics tests passed.\n";
  return EXIT_SUCCESS;
}
