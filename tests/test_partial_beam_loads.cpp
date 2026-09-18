#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/loads/BeamPartialLinearLoad3D.hpp"
#include "fem/recorder/ElementRecorder.hpp"
#include "fem/solver/LinearStaticSolver.hpp"

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

void testPartialUniformLoad() {
  constexpr double length = 6.0;
  constexpr double a = 2.0;
  constexpr double b = 5.0;
  constexpr double q = -1200.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);

  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(200.0e9, 0.3, 7850.0),
      fem::BeamSection(0.01, 6.0e-6, 8.0e-6, 1.0e-5));
  model.addElementLoad<fem::BeamPartialLinearLoad3D>(
      1, a, b,
      Eigen::Vector3d(0.0, q, 0.0),
      Eigen::Vector3d(0.0, q, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto response = fem::ElementRecorder{}.record(model, 1, result);

  const double total = q * (b - a);
  const double moment_i = q * (b * b - a * a) / 2.0;

  near(result.reactionAt(1, fem::Dof::UY), -total, 1.0e-10,
       "partial uniform root shear");
  near(result.reactionAt(1, fem::Dof::RZ), -moment_i, 1.0e-10,
       "partial uniform root moment");
  near(response.local_end_force[1], -total, 1.0e-10,
       "partial uniform recovered root shear");
  near(response.local_end_force[5], -moment_i, 1.0e-10,
       "partial uniform recovered root moment");
  near(response.local_end_force[7], 0.0, 1.0e-8,
       "partial uniform free-end shear");
  near(response.local_end_force[11], 0.0, 1.0e-8,
       "partial uniform free-end moment");
}

void testPartialTrapezoidalLoad() {
  constexpr double length = 5.0;
  constexpr double a = 1.0;
  constexpr double b = 4.0;
  constexpr double qa = -800.0;
  constexpr double qb = -2600.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);

  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(200.0e9, 0.3, 7850.0),
      fem::BeamSection(0.01, 6.0e-6, 8.0e-6, 1.0e-5));
  model.addElementLoad<fem::BeamPartialLinearLoad3D>(
      1, a, b,
      Eigen::Vector3d(0.0, qa, 0.0),
      Eigen::Vector3d(0.0, qb, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto response = fem::ElementRecorder{}.record(model, 1, result);

  const double span = b - a;
  const double total = span * (qa + qb) / 2.0;
  const double moment_i =
      a * total + span * span * (qa + 2.0 * qb) / 6.0;

  near(result.reactionAt(1, fem::Dof::UY), -total, 1.0e-10,
       "partial trapezoidal root shear");
  near(result.reactionAt(1, fem::Dof::RZ), -moment_i, 1.0e-10,
       "partial trapezoidal root moment");
  near(response.local_end_force[1], -total, 1.0e-10,
       "partial trapezoidal recovered root shear");
  near(response.local_end_force[5], -moment_i, 1.0e-10,
       "partial trapezoidal recovered root moment");
}

void testInvalidPartialInterval() {
  fem::Model model;
  model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {4.0, 0.0, 0.0});
  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(200.0e9, 0.3, 7850.0),
      fem::BeamSection(0.01, 6.0e-6, 8.0e-6, 1.0e-5));
  model.addElementLoad<fem::BeamPartialLinearLoad3D>(
      1, 3.0, 2.0,
      Eigen::Vector3d(0.0, -1000.0, 0.0),
      Eigen::Vector3d(0.0, -1000.0, 0.0));

  bool threw = false;
  try {
    (void)model.assemble();
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  if (!threw) {
    std::cerr << "FAIL: invalid partial-load interval was accepted\n";
    std::exit(EXIT_FAILURE);
  }
}

}  // namespace

int main() {
  testPartialUniformLoad();
  testPartialTrapezoidalLoad();
  testInvalidPartialInterval();
  std::cout << "All partial Beam3D load tests passed.\n";
  return EXIT_SUCCESS;
}
