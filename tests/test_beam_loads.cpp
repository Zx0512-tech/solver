#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/loads/BeamLinearLoad3D.hpp"
#include "fem/loads/BeamPointLoad3D.hpp"
#include "fem/loads/BeamUniformLoad3D.hpp"
#include "fem/recorder/ElementRecorder.hpp"
#include "fem/solver/LinearStaticSolver.hpp"

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

fem::BeamSection testSection(double iy, double iz) {
  return fem::BeamSection(0.01, iy, iz, 1.0e-5);
}

void testInteriorPointForceLocalY() {
  const double e = 200.0e9;
  const double iz = 8.0e-6;
  const double length = 5.0;
  const double a = 2.0;
  const double p = -15000.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);

  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(e, 0.3, 7850.0),
      testSection(6.0e-6, iz));
  model.addElementLoad<fem::BeamPointLoad3D>(
      1, a, Eigen::Vector3d(0.0, p, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto response = fem::ElementRecorder{}.record(model, 1, result);

  near(result.displacementAt(2, fem::Dof::UY),
       p * a * a * (3.0 * length - a) / (6.0 * e * iz),
       1.0e-10,
       "interior point-load tip y displacement");
  near(result.displacementAt(2, fem::Dof::RZ),
       p * a * a / (2.0 * e * iz),
       1.0e-10,
       "interior point-load tip z rotation");
  near(result.reactionAt(1, fem::Dof::UY), -p, 1.0e-10,
       "interior point-load root shear");
  near(result.reactionAt(1, fem::Dof::RZ), -p * a, 1.0e-10,
       "interior point-load root moment");

  near(response.local_end_force[1], -p, 1.0e-10,
       "point-load recovered root shear");
  near(response.local_end_force[5], -p * a, 1.0e-10,
       "point-load recovered root moment");
  near(response.local_end_force[7], 0.0, 1.0e-8,
       "point-load recovered free-end shear");
  near(response.local_end_force[11], 0.0, 1.0e-8,
       "point-load recovered free-end moment");
}

void testInteriorPointForceLocalZ() {
  const double e = 200.0e9;
  const double iy = 6.0e-6;
  const double length = 4.0;
  const double a = 1.5;
  const double p = 10000.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);

  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(e, 0.3, 7850.0),
      testSection(iy, 8.0e-6));
  model.addElementLoad<fem::BeamPointLoad3D>(
      1, a, Eigen::Vector3d(0.0, 0.0, p));

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto response = fem::ElementRecorder{}.record(model, 1, result);

  near(result.displacementAt(2, fem::Dof::UZ),
       p * a * a * (3.0 * length - a) / (6.0 * e * iy),
       1.0e-10,
       "interior local-z point-load tip displacement");
  near(result.displacementAt(2, fem::Dof::RY),
       -p * a * a / (2.0 * e * iy),
       1.0e-10,
       "interior local-z point-load tip rotation");
  near(result.reactionAt(1, fem::Dof::UZ), -p, 1.0e-10,
       "interior local-z point-load root shear");
  near(result.reactionAt(1, fem::Dof::RY), p * a, 1.0e-10,
       "interior local-z point-load root moment");

  near(response.local_end_force[2], -p, 1.0e-10,
       "local-z point-load recovered root shear");
  near(response.local_end_force[4], p * a, 1.0e-10,
       "local-z point-load recovered root moment");
  near(response.local_end_force[8], 0.0, 1.0e-8,
       "local-z point-load recovered free-end shear");
  near(response.local_end_force[10], 0.0, 1.0e-8,
       "local-z point-load recovered free-end moment");
}

void testInteriorPointMoment() {
  const double e = 200.0e9;
  const double iz = 8.0e-6;
  const double length = 5.0;
  const double a = 2.0;
  const double moment_z = 6000.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);

  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(e, 0.3, 7850.0),
      testSection(6.0e-6, iz));
  model.addElementLoad<fem::BeamPointLoad3D>(
      1, a,
      Eigen::Vector3d::Zero(),
      Eigen::Vector3d(0.0, 0.0, moment_z));

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto response = fem::ElementRecorder{}.record(model, 1, result);

  near(result.displacementAt(2, fem::Dof::UY),
       moment_z * a * (2.0 * length - a) / (2.0 * e * iz),
       1.0e-10,
       "interior point-moment tip displacement");
  near(result.displacementAt(2, fem::Dof::RZ),
       moment_z * a / (e * iz),
       1.0e-10,
       "interior point-moment tip rotation");
  near(result.reactionAt(1, fem::Dof::UY), 0.0, 1.0e-8,
       "interior point-moment root shear");
  near(result.reactionAt(1, fem::Dof::RZ), -moment_z, 1.0e-10,
       "interior point-moment root moment");

  near(response.local_end_force[1], 0.0, 1.0e-8,
       "point-moment recovered root shear");
  near(response.local_end_force[5], -moment_z, 1.0e-10,
       "point-moment recovered root moment");
}

void testLinearTrapezoidalLoad() {
  const double e = 200.0e9;
  const double length = 4.0;
  const double qi = -1000.0;
  const double qj = -3000.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);

  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(e, 0.3, 7850.0),
      testSection(6.0e-6, 8.0e-6));
  model.addElementLoad<fem::BeamLinearLoad3D>(
      1,
      Eigen::Vector3d(0.0, qi, 0.0),
      Eigen::Vector3d(0.0, qj, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto response = fem::ElementRecorder{}.record(model, 1, result);

  const double total_load = length * (qi + qj) / 2.0;
  const double load_moment_about_i =
      length * length * (qi + 2.0 * qj) / 6.0;

  near(result.reactionAt(1, fem::Dof::UY), -total_load, 1.0e-10,
       "trapezoidal-load root shear");
  near(result.reactionAt(1, fem::Dof::RZ),
       -load_moment_about_i,
       1.0e-10,
       "trapezoidal-load root moment");

  near(response.local_end_force[1], -total_load, 1.0e-10,
       "trapezoidal-load recovered root shear");
  near(response.local_end_force[5],
       -load_moment_about_i,
       1.0e-10,
       "trapezoidal-load recovered root moment");
  near(response.local_end_force[7], 0.0, 1.0e-8,
       "trapezoidal-load recovered free-end shear");
  near(response.local_end_force[11], 0.0, 1.0e-8,
       "trapezoidal-load recovered free-end moment");
}

void testUniformLoadInGlobalAxes() {
  const Eigen::Vector3d end(3.0, 4.0, 0.0);
  const double length = end.norm();
  const Eigen::Vector3d q_global(0.0, 0.0, -1200.0);

  fem::Model model;
  model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, end);
  model.addElement<fem::Beam3D>(
      1, 1, 2,
      fem::LinearElasticMaterial(200.0e9, 0.3, 7850.0),
      testSection(6.0e-6, 8.0e-6),
      Eigen::Vector3d::UnitZ());
  model.addElementLoad<fem::BeamUniformLoad3D>(
      1,
      q_global,
      fem::BeamLoadCoordinateSystem::Global);

  const auto system = model.assemble();
  Eigen::Vector3d resultant = Eigen::Vector3d::Zero();
  for (const fem::NodeId node_id : {1, 2}) {
    resultant.x() += system.load[
        static_cast<Eigen::Index>(system.dofs.equation(node_id, fem::Dof::UX))];
    resultant.y() += system.load[
        static_cast<Eigen::Index>(system.dofs.equation(node_id, fem::Dof::UY))];
    resultant.z() += system.load[
        static_cast<Eigen::Index>(system.dofs.equation(node_id, fem::Dof::UZ))];
  }

  near(resultant.x(), q_global.x() * length, 1.0e-10,
       "global line-load resultant x");
  near(resultant.y(), q_global.y() * length, 1.0e-10,
       "global line-load resultant y");
  near(resultant.z(), q_global.z() * length, 1.0e-10,
       "global line-load resultant z");
}

}  // namespace

int main() {
  testInteriorPointForceLocalY();
  testInteriorPointForceLocalZ();
  testInteriorPointMoment();
  testLinearTrapezoidalLoad();
  testUniformLoadInGlobalAxes();
  std::cout << "All Beam3D load tests passed.\n";
  return EXIT_SUCCESS;
}
