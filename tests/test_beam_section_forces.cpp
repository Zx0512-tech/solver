#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/loads/BeamPartialLinearLoad3D.hpp"
#include "fem/loads/BeamPointLoad3D.hpp"
#include "fem/loads/BeamUniformLoad3D.hpp"
#include "fem/response/BeamSectionForces.hpp"
#include "fem/solver/LinearStaticSolver.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
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

fem::BeamSection section() {
  return fem::BeamSection(0.01, 6.0e-6, 8.0e-6, 1.0e-5);
}

fem::LinearElasticMaterial material() {
  return fem::LinearElasticMaterial(200.0e9, 0.3, 7850.0);
}

void testCombinedTipLoads() {
  constexpr double l = 4.0;
  constexpr double x = 1.5;
  constexpr double fx = 5000.0;
  constexpr double fy = -12000.0;
  constexpr double fz = 7000.0;
  constexpr double mx = 3000.0;
  constexpr double my = 4000.0;
  constexpr double mz = -2000.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {l, 0.0, 0.0});
  fixAll(root);
  model.addElement<fem::Beam3D>(1, 1, 2, material(), section());

  tip.addLoad(fem::Dof::UX, fx);
  tip.addLoad(fem::Dof::UY, fy);
  tip.addLoad(fem::Dof::UZ, fz);
  tip.addLoad(fem::Dof::RX, mx);
  tip.addLoad(fem::Dof::RY, my);
  tip.addLoad(fem::Dof::RZ, mz);

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto forces = model.beamSectionForces(1, x, result.displacement);

  near(forces.x, x, 1.0e-12, "section coordinate");
  near(forces.N, fx, 1.0e-10, "tip-load axial force");
  near(forces.Vy, fy, 1.0e-10, "tip-load local-y shear");
  near(forces.Vz, fz, 1.0e-10, "tip-load local-z shear");
  near(forces.T, mx, 1.0e-10, "tip-load torsion");
  near(forces.My, my - fz * (l - x), 1.0e-10, "tip-load My");
  near(forces.Mz, mz + fy * (l - x), 1.0e-10, "tip-load Mz");
}

void testPointLoadJumpHasLeftAndRightLimits() {
  constexpr double l = 4.0;
  constexpr double a = 1.5;
  constexpr double p = -10000.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {l, 0.0, 0.0});
  fixAll(root);
  model.addElement<fem::Beam3D>(1, 1, 2, material(), section());
  model.addElementLoad<fem::BeamPointLoad3D>(
      1, a, Eigen::Vector3d(0.0, p, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);

  const auto at_root = model.beamSectionForces(1, 0.0, result.displacement);
  near(at_root.Vy, p, 1.0e-10, "point-load root shear");
  near(at_root.Mz, p * a, 1.0e-10, "point-load root moment");

  const auto left = model.beamSectionForces(
      1, a, result.displacement, 0.0, fem::BeamSectionSide::Left);
  const auto right = model.beamSectionForces(
      1, a, result.displacement, 0.0, fem::BeamSectionSide::Right);

  near(left.Vy, p, 1.0e-10, "point-load left-limit shear");
  near(right.Vy, 0.0, 1.0e-8, "point-load right-limit shear");
  near(left.Mz, 0.0, 1.0e-8, "point-load left-limit moment");
  near(right.Mz, 0.0, 1.0e-8, "point-load right-limit moment");

  const auto after = model.beamSectionForces(1, 3.0, result.displacement);
  near(after.Vy, 0.0, 1.0e-8, "point-load shear after load");
  near(after.Mz, 0.0, 1.0e-8, "point-load moment after load");
}

void testUniformLoadDiagram() {
  constexpr double l = 4.0;
  constexpr double q = -2000.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {l, 0.0, 0.0});
  fixAll(root);
  model.addElement<fem::Beam3D>(1, 1, 2, material(), section());
  model.addElementLoad<fem::BeamUniformLoad3D>(
      1, Eigen::Vector3d(0.0, q, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);

  for (const double x : {0.0, 1.0, 2.5, l}) {
    const auto forces = model.beamSectionForces(1, x, result.displacement);
    near(forces.Vy, q * (l - x), 1.0e-10, "uniform-load shear");
    near(forces.Mz, q * (l - x) * (l - x) / 2.0, 1.0e-10,
         "uniform-load bending moment");
  }
}

void testPartialUniformLoadDiagram() {
  constexpr double l = 5.0;
  constexpr double a = 1.0;
  constexpr double b = 4.0;
  constexpr double q = -1000.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {l, 0.0, 0.0});
  fixAll(root);
  model.addElement<fem::Beam3D>(1, 1, 2, material(), section());
  model.addElementLoad<fem::BeamPartialLinearLoad3D>(
      1, a, b,
      Eigen::Vector3d(0.0, q, 0.0),
      Eigen::Vector3d(0.0, q, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);

  {
    constexpr double x = 0.5;
    const auto forces = model.beamSectionForces(1, x, result.displacement);
    near(forces.Vy, q * (b - a), 1.0e-10,
         "partial-load shear before loaded interval");
    near(forces.Mz,
         q * ((b * b - a * a) / 2.0 - x * (b - a)),
         1.0e-10,
         "partial-load moment before loaded interval");
  }

  {
    constexpr double x = 2.5;
    const auto forces = model.beamSectionForces(1, x, result.displacement);
    near(forces.Vy, q * (b - x), 1.0e-10,
         "partial-load shear inside loaded interval");
    near(forces.Mz, q * (b - x) * (b - x) / 2.0, 1.0e-10,
         "partial-load moment inside loaded interval");
  }

  {
    const auto forces = model.beamSectionForces(1, 4.5, result.displacement);
    near(forces.Vy, 0.0, 1.0e-8, "partial-load shear after interval");
    near(forces.Mz, 0.0, 1.0e-8, "partial-load moment after interval");
  }
}

void testTimeDependentElementLoadScalesSectionForces() {
  constexpr double l = 4.0;
  constexpr double a = 2.0;
  constexpr double p = -5000.0;

  fem::Model model;
  model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {l, 0.0, 0.0});
  model.addElement<fem::Beam3D>(1, 1, 2, material(), section());
  model.addTimeDependentElementLoad<fem::BeamPointLoad3D>(
      [](double t) { return t; },
      1, a, Eigen::Vector3d(0.0, p, 0.0));

  const auto assembled = model.assemble();
  const Eigen::VectorXd zero =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(assembled.dofs.size()));

  const auto f1 = model.beamSectionForces(1, 1.0, zero, 1.0);
  const auto f2 = model.beamSectionForces(1, 1.0, zero, 2.0);

  near(f2.Vy, 2.0 * f1.Vy, 1.0e-10,
       "time-dependent point-load shear scaling");
  near(f2.Mz, 2.0 * f1.Mz, 1.0e-10,
       "time-dependent point-load moment scaling");
}

void testSectionCoordinateValidation() {
  fem::Model model;
  model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {3.0, 0.0, 0.0});
  model.addElement<fem::Beam3D>(1, 1, 2, material(), section());
  const auto assembled = model.assemble();
  const Eigen::VectorXd zero =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(assembled.dofs.size()));

  for (const double x : {-0.1, 3.1}) {
    bool threw = false;
    try {
      (void)model.beamSectionForces(1, x, zero);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    if (!threw) {
      std::cerr << "FAIL: invalid section coordinate was accepted\n";
      std::exit(EXIT_FAILURE);
    }
  }
}

}  // namespace

int main() {
  testCombinedTipLoads();
  testPointLoadJumpHasLeftAndRightLimits();
  testUniformLoadDiagram();
  testPartialUniformLoadDiagram();
  testTimeDependentElementLoadScalesSectionForces();
  testSectionCoordinateValidation();
  std::cout << "All Beam3D section-force tests passed.\n";
  return EXIT_SUCCESS;
}
