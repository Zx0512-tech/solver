#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/dynamics/RayleighDamping.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/loads/BeamPartialLinearLoad3D.hpp"
#include "fem/loads/BeamPointLoad3D.hpp"
#include "fem/loads/TimeDependentElementLoad.hpp"
#include "fem/recorder/ElementRecorder.hpp"
#include "fem/solver/NewmarkBetaSolver.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

void near(double actual, double expected, double tol, const std::string& message) {
  const double scale = std::max(1.0, std::abs(expected));
  if (std::abs(actual - expected) > tol * scale) {
    std::cerr << "FAIL: " << message << " actual=" << actual
              << " expected=" << expected << '\n';
    std::exit(EXIT_FAILURE);
  }
}

void fixAll(fem::Node& node) {
  for (std::size_t i = 0; i < fem::kDofsPerFrameNode; ++i) {
    node.fix(fem::dofFromOffset(i));
  }
}

fem::Model makeAxialModel(double length) {
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
      1, 1, 2,
      fem::LinearElasticMaterial(210.0e9, 0.3, 7850.0),
      fem::BeamSection(0.01, 2.0e-5, 3.0e-5, 4.0e-5));
  return model;
}

void testTimeScaledPointLoadVector() {
  constexpr double length = 2.0;
  constexpr double p0 = 10000.0;
  fem::Model model = makeAxialModel(length);

  model.addTimeDependentElementLoad<fem::BeamPointLoad3D>(
      [](double t) { return 2.0 * t; },
      1, length, Eigen::Vector3d(p0, 0.0, 0.0));

  const auto f0 = model.loadVector(0.0);
  const auto f025 = model.loadVector(0.25);
  const auto f050 = model.loadVector(0.50);
  const fem::DofManager dofs(model.nodes());
  const Eigen::Index tip_ux =
      static_cast<Eigen::Index>(dofs.equation(2, fem::Dof::UX));

  near(f0[tip_ux], 0.0, 1.0e-12, "time point load at t=0");
  near(f025[tip_ux], 0.5 * p0, 1.0e-12, "time point load at t=0.25");
  near(f050[tip_ux], p0, 1.0e-12, "time point load at t=0.5");
}

void testNewmarkTimeElementLoadMatchesNodalTimeLoad() {
  constexpr double length = 2.0;
  constexpr double p0 = 8000.0;
  constexpr double omega = 30.0;

  fem::Model element_model = makeAxialModel(length);
  element_model.addTimeDependentElementLoad<fem::BeamPointLoad3D>(
      [omega](double t) { return std::sin(omega * t); },
      1, length, Eigen::Vector3d(p0, 0.0, 0.0));

  fem::Model nodal_model = makeAxialModel(length);

  fem::NewmarkSettings settings;
  settings.time_step = 0.0005;
  settings.step_count = 200;

  const auto element_result =
      fem::NewmarkBetaSolver{}.solve(element_model, settings);

  const std::vector<fem::NodalTimeLoad> nodal_loads = {
      {2, fem::Dof::UX,
       [p0, omega](double t) { return p0 * std::sin(omega * t); }}};
  const auto nodal_result =
      fem::NewmarkBetaSolver{}.solve(
          nodal_model, settings, fem::RayleighDamping{}, nodal_loads);

  for (std::size_t step : {0U, 40U, 100U, 200U}) {
    near(element_result.displacementAt(step, 2, fem::Dof::UX),
         nodal_result.displacementAt(step, 2, fem::Dof::UX),
         1.0e-10,
         "time-dependent element load matches nodal load");
  }
}

void testRecorderUsesCurrentElementLoadTime() {
  constexpr double length = 2.0;
  constexpr double p0 = 5000.0;
  fem::Model model = makeAxialModel(length);

  model.addTimeDependentElementLoad<fem::BeamPointLoad3D>(
      [](double t) { return t; },
      1, length, Eigen::Vector3d(p0, 0.0, 0.0));

  fem::NewmarkSettings settings;
  settings.time_step = 0.01;
  settings.step_count = 4;
  const auto result = fem::NewmarkBetaSolver{}.solve(model, settings);
  const auto history = fem::ElementRecorder{}.record(model, 1, result);

  const std::size_t step = 4;
  const auto expected =
      model.elementResponse(
          1,
          result.displacement.col(static_cast<Eigen::Index>(step)),
          result.time[step]);

  near(history.local_end_force(6, static_cast<Eigen::Index>(step)),
       expected.local_end_force[6],
       1.0e-12,
       "recorder evaluates element load at current time");
}

void testTimeScaledPartialDistributedLoad() {
  constexpr double length = 5.0;
  fem::Model model;
  model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {length, 0.0, 0.0});
  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(200.0e9, 0.3, 7850.0),
      fem::BeamSection(0.01, 6.0e-6, 8.0e-6, 1.0e-5));

  constexpr double a = 1.0;
  constexpr double b = 4.0;
  constexpr double q = -1000.0;
  model.addTimeDependentElementLoad<fem::BeamPartialLinearLoad3D>(
      [](double t) { return 1.0 + t; },
      1, a, b,
      Eigen::Vector3d(0.0, q, 0.0),
      Eigen::Vector3d(0.0, q, 0.0));

  const auto f0 = model.loadVector(0.0);
  const auto f2 = model.loadVector(2.0);
  near(f2.sum(), 3.0 * f0.sum(), 1.0e-12,
       "time scaling applies to partial distributed load");
}

}  // namespace

int main() {
  testTimeScaledPointLoadVector();
  testNewmarkTimeElementLoadMatchesNodalTimeLoad();
  testRecorderUsesCurrentElementLoadTime();
  testTimeScaledPartialDistributedLoad();
  std::cout << "All time-dependent element-load tests passed.\n";
  return EXIT_SUCCESS;
}
