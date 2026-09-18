#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/loads/BeamUniformLoad3D.hpp"
#include "fem/recorder/ElementRecorder.hpp"
#include "fem/solver/LinearStaticSolver.hpp"
#include "fem/solver/NewmarkBetaSolver.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

void near(double actual,
          double expected,
          double relative_tolerance,
          const std::string& message) {
  const double scale = std::max(1.0, std::abs(expected));
  if (std::abs(actual - expected) > relative_tolerance * scale) {
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

void testNodalLoadResponse() {
  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {3.0, 0.0, 0.0});
  fixAll(root);

  const double e = 200.0e9;
  const double iz = 8.0e-6;
  const double length = 3.0;
  const double p = -12000.0;

  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(e, 0.3, 7850.0),
      fem::BeamSection(0.01, 6.0e-6, iz, 1.0e-5));
  tip.addLoad(fem::Dof::UY, p);

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto response = fem::ElementRecorder{}.record(model, 1, result);

  near(response.local_end_force[1], -p, 1.0e-10,
       "nodal-load root shear");
  near(response.local_end_force[5], -p * length, 1.0e-10,
       "nodal-load root bending moment");
  near(response.local_end_force[7], p, 1.0e-10,
       "nodal-load tip shear");
  near(response.local_end_force[11], 0.0, 1.0e-8,
       "nodal-load tip moment");
}

void testUniformLoadAndFixedEndForceRecovery() {
  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {4.0, 0.0, 0.0});
  fixAll(root);

  const double e = 200.0e9;
  const double iz = 8.0e-6;
  const double length = 4.0;
  const double qy = -2000.0;

  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(e, 0.3, 7850.0),
      fem::BeamSection(0.01, 6.0e-6, iz, 1.0e-5));
  model.addElementLoad<fem::BeamUniformLoad3D>(
      1, Eigen::Vector3d(0.0, qy, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto response = fem::ElementRecorder{}.record(model, 1, result);

  near(result.displacementAt(2, fem::Dof::UY),
       qy * std::pow(length, 4) / (8.0 * e * iz),
       1.0e-10,
       "uniform-load cantilever tip deflection");
  near(result.displacementAt(2, fem::Dof::RZ),
       qy * std::pow(length, 3) / (6.0 * e * iz),
       1.0e-10,
       "uniform-load cantilever tip rotation");

  near(result.reactionAt(1, fem::Dof::UY), -qy * length, 1.0e-10,
       "uniform-load support shear");
  near(result.reactionAt(1, fem::Dof::RZ),
       -qy * length * length / 2.0,
       1.0e-10,
       "uniform-load support moment");

  // A free beam end under only a distributed load has zero section resultants.
  near(response.local_end_force[1], -qy * length, 1.0e-10,
       "recovered root shear includes fixed-end contribution");
  near(response.local_end_force[5],
       -qy * length * length / 2.0,
       1.0e-10,
       "recovered root moment includes fixed-end contribution");
  near(response.local_end_force[7], 0.0, 1.0e-8,
       "recovered free-end shear");
  near(response.local_end_force[11], 0.0, 1.0e-8,
       "recovered free-end moment");
}

void testDynamicElementRecorder() {
  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {2.0, 0.0, 0.0});
  fixAll(root);

  tip.fix(fem::Dof::UY);
  tip.fix(fem::Dof::UZ);
  tip.fix(fem::Dof::RX);
  tip.fix(fem::Dof::RY);
  tip.fix(fem::Dof::RZ);

  const double e = 210.0e9;
  const double area = 0.01;
  const double length = 2.0;
  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(e, 0.3, 7850.0),
      fem::BeamSection(area, 2.0e-5, 3.0e-5, 4.0e-5));

  const auto assembled = model.assemble();
  fem::DynamicInitialState initial;
  initial.displacement =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(assembled.dofs.size()));
  initial.velocity =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(assembled.dofs.size()));

  const double initial_u = 1.0e-4;
  initial.displacement[
      static_cast<Eigen::Index>(assembled.dofs.equation(2, fem::Dof::UX))] =
      initial_u;

  fem::NewmarkSettings settings;
  settings.time_step = 1.0e-4;
  settings.step_count = 10;

  const auto dynamic =
      fem::NewmarkBetaSolver{}.solve(model, settings, fem::RayleighDamping{}, {}, initial);
  const auto history = fem::ElementRecorder{}.record(model, 1, dynamic);

  if (history.local_end_force.rows() != 12 ||
      history.local_end_force.cols() != 11) {
    std::cerr << "FAIL: unexpected element force history dimensions\n";
    std::exit(EXIT_FAILURE);
  }

  const double axial_force = e * area * initial_u / length;
  near(history.local_end_force(0, 0), -axial_force, 1.0e-10,
       "dynamic recorder initial root axial force");
  near(history.local_end_force(6, 0), axial_force, 1.0e-10,
       "dynamic recorder initial tip axial force");
}

}  // namespace

int main() {
  testNodalLoadResponse();
  testUniformLoadAndFixedEndForceRecovery();
  testDynamicElementRecorder();
  std::cout << "All element-response tests passed.\n";
  return EXIT_SUCCESS;
}
