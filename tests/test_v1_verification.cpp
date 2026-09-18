#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/dynamics/RayleighDamping.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/loads/BeamUniformLoad3D.hpp"
#include "fem/recorder/EnvelopeRecorder.hpp"
#include "fem/recorder/NodeRecorder.hpp"
#include "fem/recorder/SectionRecorder.hpp"
#include "fem/solver/LinearStaticSolver.hpp"
#include "fem/solver/ModalSolver.hpp"
#include "fem/solver/NewmarkBetaSolver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct VerificationState {
  std::size_t checks{0};
  double worst_relative_error{0.0};

  void checkNear(const std::string& name,
                 double actual,
                 double expected,
                 double relative_tolerance,
                 double absolute_tolerance = 1.0e-12) {
    const double scale = std::max(std::abs(expected), absolute_tolerance);
    const double absolute_error = std::abs(actual - expected);
    const double relative_error = absolute_error / scale;
    worst_relative_error = std::max(worst_relative_error, relative_error);
    ++checks;

    std::cout << std::left << std::setw(44) << name
              << " actual=" << std::setw(16) << actual
              << " expected=" << std::setw(16) << expected
              << " rel_err=" << relative_error
              << " tol=" << relative_tolerance << '\n';

    if (absolute_error > absolute_tolerance &&
        relative_error > relative_tolerance) {
      std::cerr << "FAIL: " << name << '\n';
      std::exit(EXIT_FAILURE);
    }
  }

  void checkTrue(const std::string& name, bool condition) {
    ++checks;
    std::cout << std::left << std::setw(44) << name
              << (condition ? " PASS" : " FAIL") << '\n';
    if (!condition) {
      std::cerr << "FAIL: " << name << '\n';
      std::exit(EXIT_FAILURE);
    }
  }
};

void fixAll(fem::Node& node) {
  for (std::size_t i = 0; i < fem::kDofsPerFrameNode; ++i) {
    node.fix(fem::dofFromOffset(i));
  }
}

void constrainAxialOnly(fem::Node& node) {
  node.fix(fem::Dof::UY);
  node.fix(fem::Dof::UZ);
  node.fix(fem::Dof::RX);
  node.fix(fem::Dof::RY);
  node.fix(fem::Dof::RZ);
}

void constrainPlanarYBending(fem::Node& node) {
  node.fix(fem::Dof::UX);
  node.fix(fem::Dof::UZ);
  node.fix(fem::Dof::RX);
  node.fix(fem::Dof::RY);
}

void verifyCombinedStaticCantilever(VerificationState& state) {
  constexpr double e = 210.0e9;
  constexpr double nu = 0.30;
  constexpr double area = 0.015;
  constexpr double iy = 7.0e-6;
  constexpr double iz = 1.1e-5;
  constexpr double j = 1.4e-5;
  constexpr double length = 3.2;
  constexpr double fx = 45000.0;
  constexpr double fy = -12000.0;
  constexpr double fz = 8000.0;
  constexpr double mx = 3500.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);

  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(e, nu, 7850.0),
      fem::BeamSection(area, iy, iz, j));

  tip.addLoad(fem::Dof::UX, fx);
  tip.addLoad(fem::Dof::UY, fy);
  tip.addLoad(fem::Dof::UZ, fz);
  tip.addLoad(fem::Dof::RX, mx);

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const double g = e / (2.0 * (1.0 + nu));

  state.checkNear(
      "static axial tip displacement",
      result.displacementAt(2, fem::Dof::UX),
      fx * length / (e * area),
      1.0e-10);
  state.checkNear(
      "static local-y tip displacement",
      result.displacementAt(2, fem::Dof::UY),
      fy * std::pow(length, 3) / (3.0 * e * iz),
      1.0e-10);
  state.checkNear(
      "static local-z tip displacement",
      result.displacementAt(2, fem::Dof::UZ),
      fz * std::pow(length, 3) / (3.0 * e * iy),
      1.0e-10);
  state.checkNear(
      "static torsional tip rotation",
      result.displacementAt(2, fem::Dof::RX),
      mx * length / (g * j),
      1.0e-10);

  state.checkNear(
      "static root axial reaction",
      result.reactionAt(1, fem::Dof::UX),
      -fx,
      1.0e-10);
  state.checkNear(
      "static root y reaction",
      result.reactionAt(1, fem::Dof::UY),
      -fy,
      1.0e-10);
  state.checkNear(
      "static root z reaction",
      result.reactionAt(1, fem::Dof::UZ),
      -fz,
      1.0e-10);
  state.checkNear(
      "static root torsional reaction",
      result.reactionAt(1, fem::Dof::RX),
      -mx,
      1.0e-10);
  state.checkNear(
      "static root bending reaction RY",
      result.reactionAt(1, fem::Dof::RY),
      fz * length,
      1.0e-10);
  state.checkNear(
      "static root bending reaction RZ",
      result.reactionAt(1, fem::Dof::RZ),
      -fy * length,
      1.0e-10);
}

void verifyDistributedLoadSectionEquilibrium(VerificationState& state) {
  constexpr double e = 200.0e9;
  constexpr double iz = 9.0e-6;
  constexpr double length = 4.0;
  constexpr double q = -2200.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);

  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(e, 0.3, 7850.0),
      fem::BeamSection(0.012, 6.0e-6, iz, 1.0e-5));
  model.addElementLoad<fem::BeamUniformLoad3D>(
      1, Eigen::Vector3d(0.0, q, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);

  state.checkNear(
      "UDL cantilever tip deflection",
      result.displacementAt(2, fem::Dof::UY),
      q * std::pow(length, 4) / (8.0 * e * iz),
      1.0e-10);
  state.checkNear(
      "UDL cantilever tip rotation",
      result.displacementAt(2, fem::Dof::RZ),
      q * std::pow(length, 3) / (6.0 * e * iz),
      1.0e-10);

  constexpr double x = 1.25;
  const auto forces =
      model.beamSectionForces(1, x, result.displacement);
  state.checkNear(
      "UDL section shear",
      forces.Vy,
      q * (length - x),
      1.0e-10);
  state.checkNear(
      "UDL section bending moment",
      forces.Mz,
      q * std::pow(length - x, 2) / 2.0,
      1.0e-10);
}

void verifyModalFrequency(VerificationState& state) {
  constexpr double e = 210.0e9;
  constexpr double rho = 7850.0;
  constexpr double area = 0.02;
  constexpr double iz = 8.0e-5;
  constexpr double length = 4.0;
  constexpr int elements = 8;

  fem::Model model;
  for (int i = 0; i <= elements; ++i) {
    auto& node = model.addNode(
        static_cast<fem::NodeId>(i + 1),
        {length * static_cast<double>(i) / elements, 0.0, 0.0});
    if (i == 0) {
      fixAll(node);
    } else {
      constrainPlanarYBending(node);
    }
  }

  const fem::LinearElasticMaterial material(e, 0.3, rho);
  const fem::BeamSection section(area, 5.0e-5, iz, 7.0e-5);
  for (int i = 0; i < elements; ++i) {
    model.addElement<fem::Beam3D>(
        static_cast<fem::ElementId>(i + 1),
        static_cast<fem::NodeId>(i + 1),
        static_cast<fem::NodeId>(i + 2),
        material,
        section);
  }

  const auto modal = fem::ModalSolver{}.solve(model, 1);
  constexpr double beta1 = 1.875104068711961;
  const double expected =
      beta1 * beta1 *
      std::sqrt(e * iz / (rho * area * std::pow(length, 4)));

  state.checkNear(
      "cantilever first bending omega",
      modal.angular_frequencies.at(0),
      expected,
      3.0e-5);
}

struct AxialDynamicFixture {
  fem::Model model;
  double omega{};
  double initial_displacement{};
};

AxialDynamicFixture makeAxialDynamicFixture() {
  constexpr double e = 210.0e9;
  constexpr double rho = 7850.0;
  constexpr double area = 0.01;
  constexpr double length = 2.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);
  constrainAxialOnly(tip);

  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(e, 0.3, rho),
      fem::BeamSection(area, 2.0e-5, 3.0e-5, 4.0e-5));

  const double omega =
      std::sqrt(3.0 * e / (rho * length * length));
  return {std::move(model), omega, 1.0e-3};
}

void verifyNewmarkAndRecorderEnvelope(VerificationState& state) {
  auto fixture = makeAxialDynamicFixture();
  const auto assembled = fixture.model.assemble();

  fem::DynamicInitialState initial;
  initial.displacement =
      Eigen::VectorXd::Zero(
          static_cast<Eigen::Index>(assembled.dofs.size()));
  initial.velocity =
      Eigen::VectorXd::Zero(
          static_cast<Eigen::Index>(assembled.dofs.size()));
  initial.displacement[
      static_cast<Eigen::Index>(
          assembled.dofs.equation(2, fem::Dof::UX))] =
      fixture.initial_displacement;

  const double two_pi = 2.0 * std::acos(-1.0);
  const double period = two_pi / fixture.omega;

  fem::NewmarkSettings settings;
  settings.time_step = period / 1200.0;
  settings.step_count = 1200;
  settings.beta = 0.25;
  settings.gamma = 0.5;

  const auto result =
      fem::NewmarkBetaSolver{}.solve(
          fixture.model, settings, {}, {}, initial);

  state.checkNear(
      "Newmark full-period displacement",
      result.displacementAt(1200, 2, fem::Dof::UX),
      fixture.initial_displacement,
      2.0e-5,
      1.0e-12);
  state.checkNear(
      "Newmark half-period displacement",
      result.displacementAt(600, 2, fem::Dof::UX),
      -fixture.initial_displacement,
      2.0e-5,
      1.0e-12);

  const auto node_history =
      fem::NodeRecorder{}.record(2, result);
  const auto envelope =
      fem::EnvelopeRecorder{}.record(
          node_history,
          fem::NodeResponseQuantity::Displacement,
          fem::Dof::UX);

  state.checkNear(
      "NodeRecorder/Envelope max abs displacement",
      envelope.maximum_absolute.magnitude,
      fixture.initial_displacement,
      2.0e-5,
      1.0e-12);

  const auto section_history =
      fem::SectionRecorder{}.record(
          fixture.model, 1, 0.0, result);
  const auto axial_envelope =
      fem::EnvelopeRecorder{}.record(
          section_history,
          fem::BeamSectionForceComponent::N);

  const double expected_force =
      210.0e9 * 0.01 / 2.0 * fixture.initial_displacement;
  state.checkNear(
      "SectionRecorder axial-force envelope",
      axial_envelope.maximum_absolute.magnitude,
      expected_force,
      2.0e-5);
}

void verifyRayleighTargets(VerificationState& state) {
  constexpr double omega1 = 12.0;
  constexpr double omega2 = 40.0;
  constexpr double zeta1 = 0.02;
  constexpr double zeta2 = 0.05;

  const auto damping =
      fem::RayleighDamping::fromModalTargets(
          omega1, zeta1, omega2, zeta2);

  state.checkNear(
      "Rayleigh damping target 1",
      damping.dampingRatio(omega1),
      zeta1,
      1.0e-12);
  state.checkNear(
      "Rayleigh damping target 2",
      damping.dampingRatio(omega2),
      zeta2,
      1.0e-12);
}

void verifySectionStress(VerificationState& state) {
  constexpr double length = 3.0;
  constexpr double axial = 60000.0;
  constexpr double fy = -10000.0;

  const fem::RectangleSection section(0.20, 0.40);

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);

  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(210.0e9, 0.3, 7850.0),
      section);

  tip.addLoad(fem::Dof::UX, axial);
  tip.addLoad(fem::Dof::UY, fy);

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const double y = section.size_y / 2.0;
  const double root_mz = fy * length;
  const double expected =
      axial / section.area -
      root_mz * y / section.iz;

  state.checkNear(
      "rectangle root edge normal stress",
      model.beamNormalStressAt(
          1, 0.0, y, 0.0, result.displacement),
      expected,
      1.0e-10);
}

}  // namespace

int main() {
  std::cout << std::setprecision(12);
  std::cout << "solver v1.0 verification suite\n";
  std::cout << "--------------------------------------------\n";

  VerificationState state;
  verifyCombinedStaticCantilever(state);
  verifyDistributedLoadSectionEquilibrium(state);
  verifyModalFrequency(state);
  verifyNewmarkAndRecorderEnvelope(state);
  verifyRayleighTargets(state);
  verifySectionStress(state);

  std::cout << "--------------------------------------------\n";
  std::cout << "checks=" << state.checks
            << " worst_relative_error="
            << state.worst_relative_error << '\n';
  std::cout << "V1 verification suite passed.\n";
  return EXIT_SUCCESS;
}
