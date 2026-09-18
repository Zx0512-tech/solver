#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/recorder/EnvelopeRecorder.hpp"
#include "fem/recorder/NodeRecorder.hpp"
#include "fem/recorder/SectionRecorder.hpp"
#include "fem/solver/LinearStaticSolver.hpp"
#include "fem/solver/NewmarkBetaSolver.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void near(double actual, double expected, double tol, const std::string& message) {
  const double scale = std::max(1.0, std::abs(expected));
  if (std::abs(actual - expected) > tol * scale) {
    std::cerr << "FAIL: " << message << " actual=" << actual
              << " expected=" << expected << '\n';
    std::exit(EXIT_FAILURE);
  }
}

void require(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

void fixAll(fem::Node& node) {
  for (std::size_t i = 0; i < fem::kDofsPerFrameNode; ++i) {
    node.fix(fem::dofFromOffset(i));
  }
}

fem::LinearElasticMaterial material() {
  return fem::LinearElasticMaterial(200.0e9, 0.3, 7850.0);
}

void testNodeRecorderStatic() {
  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {2.0, 0.0, 0.0});
  fixAll(root);
  model.addElement<fem::Beam3D>(
      1, 1, 2, material(), fem::BeamSection(0.01, 1e-5, 1e-5, 1e-5));

  tip.addLoad(fem::Dof::UX, 10000.0);
  tip.addLoad(fem::Dof::UY, -5000.0);

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto node = fem::NodeRecorder{}.record(2, result);

  require(node.node_id == 2, "static node recorder id");
  for (std::size_t i = 0; i < fem::kDofsPerFrameNode; ++i) {
    const fem::Dof dof = fem::dofFromOffset(i);
    near(node.displacement[static_cast<Eigen::Index>(i)],
         result.displacementAt(2, dof), 1e-12,
         "static node displacement");
    near(node.reaction[static_cast<Eigen::Index>(i)],
         result.reactionAt(2, dof), 1e-12,
         "static node reaction");
  }
}

void testNodeRecorderDynamicAndScalarSeries() {
  fem::Model model;
  model.addNode(1, {0.0, 0.0, 0.0});
  fem::DofManager dofs(model.nodes());

  fem::NewmarkResult result;
  result.time = {0.0, 0.1, 0.2, 0.3};
  result.dofs = dofs;
  result.displacement = Eigen::MatrixXd::Zero(6, 4);
  result.velocity = Eigen::MatrixXd::Zero(6, 4);
  result.acceleration = Eigen::MatrixXd::Zero(6, 4);

  result.displacement.row(0) << 0.0, 1.0, -2.0, 0.5;
  result.velocity.row(1) << 0.0, 3.0, 4.0, -1.0;
  result.acceleration.row(2) << 2.0, -5.0, 1.0, 0.0;

  const auto history = fem::NodeRecorder{}.record(1, result);

  require(history.node_id == 1, "dynamic node recorder id");
  require(history.time == result.time, "dynamic node recorder time");
  require(history.displacement.rows() == 6 && history.displacement.cols() == 4,
          "dynamic node recorder shape");
  near(history.displacement(0, 2), -2.0, 1e-12,
       "dynamic node displacement sample");
  near(history.velocity(1, 2), 4.0, 1e-12,
       "dynamic node velocity sample");
  near(history.acceleration(2, 1), -5.0, 1e-12,
       "dynamic node acceleration sample");

  const auto series = history.series(
      fem::NodeResponseQuantity::Displacement, fem::Dof::UX);
  require(series.time == result.time, "node scalar series time");
  near(series.value[2], -2.0, 1e-12, "node scalar series value");
}

fem::Model makeAxialBeam() {
  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {2.0, 0.0, 0.0});
  fixAll(root);
  tip.fix(fem::Dof::UY);
  tip.fix(fem::Dof::UZ);
  tip.fix(fem::Dof::RX);
  tip.fix(fem::Dof::RY);
  tip.fix(fem::Dof::RZ);
  model.addElement<fem::Beam3D>(
      1, 1, 2, material(),
      fem::GeneralSection(0.01, 1e-5, 1e-5, 1e-5));
  return model;
}

void testSectionRecorderStatic() {
  fem::Model model = makeAxialBeam();
  model.node(2).addLoad(fem::Dof::UX, 12000.0);
  const auto result = fem::LinearStaticSolver{}.solve(model);

  const auto section =
      fem::SectionRecorder{}.record(model, 1, 0.7, result);

  require(section.element_id == 1, "static section element id");
  near(section.x, 0.7, 1e-12, "static section x");
  near(section.forces.N, 12000.0, 1e-10, "static section axial force");
}

void testSectionRecorderDynamicAndNormalStressHistory() {
  fem::Model model = makeAxialBeam();
  const fem::DofManager dofs(model.nodes());
  const Eigen::Index tip_ux =
      static_cast<Eigen::Index>(dofs.equation(2, fem::Dof::UX));

  fem::NewmarkResult result;
  result.time = {0.0, 0.1, 0.2, 0.3};
  result.dofs = dofs;
  result.displacement =
      Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(dofs.size()), 4);
  result.velocity = result.displacement;
  result.acceleration = result.displacement;

  result.displacement(tip_ux, 0) = 0.0;
  result.displacement(tip_ux, 1) = 1.0e-6;
  result.displacement(tip_ux, 2) = -2.0e-6;
  result.displacement(tip_ux, 3) = 0.5e-6;

  const auto history =
      fem::SectionRecorder{}.record(model, 1, 1.0, result);

  require(history.time == result.time, "section history time");
  require(history.forces.rows() == 6 && history.forces.cols() == 4,
          "section force history shape");

  const double ea_over_l = 200.0e9 * 0.01 / 2.0;
  near(history.forceAt(1, fem::BeamSectionForceComponent::N),
       ea_over_l * 1.0e-6, 1e-10, "section history N positive");
  near(history.forceAt(2, fem::BeamSectionForceComponent::N),
       ea_over_l * -2.0e-6, 1e-10, "section history N negative");

  const auto n_series =
      history.series(fem::BeamSectionForceComponent::N);
  near(n_series.value[2], ea_over_l * -2.0e-6, 1e-10,
       "section force scalar series");

  const auto stress_history =
      fem::SectionRecorder{}.recordNormalStress(
          model, 1, 1.0, 0.0, 0.0, result);
  near(stress_history.value[1],
       ea_over_l * 1.0e-6 / 0.01, 1e-10,
       "section normal stress history positive");
  near(stress_history.value[2],
       ea_over_l * -2.0e-6 / 0.01, 1e-10,
       "section normal stress history negative");
}

void testEnvelopeRecorderReportsValuesTimesAndSteps() {
  fem::ScalarHistory history;
  history.time = {0.0, 0.1, 0.2, 0.3, 0.4};
  history.value = {2.0, -5.0, 4.0, -3.0, 5.0};

  const auto envelope = fem::EnvelopeRecorder{}.record(history);

  near(envelope.minimum.value, -5.0, 1e-12, "envelope minimum value");
  near(envelope.minimum.time, 0.1, 1e-12, "envelope minimum time");
  require(envelope.minimum.step == 1, "envelope minimum step");

  near(envelope.maximum.value, 5.0, 1e-12, "envelope maximum value");
  near(envelope.maximum.time, 0.4, 1e-12, "envelope maximum time");
  require(envelope.maximum.step == 4, "envelope maximum step");

  // Tie on |value| is deterministic: first occurrence wins.
  near(envelope.maximum_absolute.magnitude, 5.0, 1e-12,
       "envelope absolute maximum magnitude");
  near(envelope.maximum_absolute.value, -5.0, 1e-12,
       "envelope absolute maximum signed value");
  near(envelope.maximum_absolute.time, 0.1, 1e-12,
       "envelope absolute maximum time");
  require(envelope.maximum_absolute.step == 1,
          "envelope absolute maximum first tie step");
}

void testEnvelopeConvenienceForNodeAndSectionHistories() {
  fem::NodeResponseHistory node;
  node.node_id = 1;
  node.time = {0.0, 0.1, 0.2};
  node.displacement = Eigen::MatrixXd::Zero(6, 3);
  node.velocity = Eigen::MatrixXd::Zero(6, 3);
  node.acceleration = Eigen::MatrixXd::Zero(6, 3);
  node.displacement.row(1) << 1.0, -4.0, 3.0;

  const auto node_env =
      fem::EnvelopeRecorder{}.record(
          node, fem::NodeResponseQuantity::Displacement, fem::Dof::UY);
  near(node_env.minimum.value, -4.0, 1e-12,
       "node envelope minimum");
  near(node_env.maximum_absolute.magnitude, 4.0, 1e-12,
       "node envelope absolute maximum");

  fem::SectionResponseHistory section;
  section.element_id = 3;
  section.x = 1.0;
  section.time = {0.0, 0.1, 0.2};
  section.forces = Eigen::MatrixXd::Zero(6, 3);
  section.forces.row(5) << -2.0, 8.0, -10.0;

  const auto section_env =
      fem::EnvelopeRecorder{}.record(
          section, fem::BeamSectionForceComponent::Mz);
  near(section_env.maximum.value, 8.0, 1e-12,
       "section envelope maximum");
  near(section_env.minimum.value, -10.0, 1e-12,
       "section envelope minimum");
  near(section_env.minimum.time, 0.2, 1e-12,
       "section envelope minimum time");
}

void testRecorderValidation() {
  {
    fem::ScalarHistory bad;
    bad.time = {0.0, 0.1};
    bad.value = {1.0};

    bool threw = false;
    try {
      (void)fem::EnvelopeRecorder{}.record(bad);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    require(threw, "envelope should reject mismatched history sizes");
  }

  {
    fem::ScalarHistory empty;
    bool threw = false;
    try {
      (void)fem::EnvelopeRecorder{}.record(empty);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    require(threw, "envelope should reject empty history");
  }
}

}  // namespace

int main() {
  testNodeRecorderStatic();
  testNodeRecorderDynamicAndScalarSeries();
  testSectionRecorderStatic();
  testSectionRecorderDynamicAndNormalStressHistory();
  testEnvelopeRecorderReportsValuesTimesAndSteps();
  testEnvelopeConvenienceForNodeAndSectionHistories();
  testRecorderValidation();
  std::cout << "All result-recorder tests passed.\n";
  return EXIT_SUCCESS;
}
