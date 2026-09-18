#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/response/BeamSectionStressRecovery.hpp"
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

void testRectangleSectionProperties() {
  constexpr double size_y = 0.20;
  constexpr double size_z = 0.40;
  fem::RectangleSection section(size_y, size_z);

  near(section.area, size_y * size_z, 1.0e-12, "rectangle area");
  near(section.iy, size_y * std::pow(size_z, 3) / 12.0, 1.0e-12,
       "rectangle Iy");
  near(section.iz, size_z * std::pow(size_y, 3) / 12.0, 1.0e-12,
       "rectangle Iz");
  require(section.torsion_constant > 0.0, "rectangle J should be positive");
  require(section.kind == fem::BeamSectionKind::Rectangle,
          "rectangle section kind");
  near(section.size_y, size_y, 1.0e-12, "rectangle size_y");
  near(section.size_z, size_z, 1.0e-12, "rectangle size_z");
}

void testCircularSectionProperties() {
  constexpr double radius = 0.15;
  const double pi = std::acos(-1.0);
  fem::CircularSection section(radius);

  near(section.area, pi * radius * radius, 1.0e-12, "circle area");
  near(section.iy, pi * std::pow(radius, 4) / 4.0, 1.0e-12, "circle Iy");
  near(section.iz, section.iy, 1.0e-12, "circle Iz");
  near(section.torsion_constant, pi * std::pow(radius, 4) / 2.0, 1.0e-12,
       "circle polar J");
  require(section.kind == fem::BeamSectionKind::SolidCircle,
          "circle section kind");
  near(section.radius, radius, 1.0e-12, "circle radius");
}

void testInvalidGeometryRejected() {
  for (const double bad : {0.0, -0.1}) {
    bool rectangle_threw = false;
    try {
      (void)fem::RectangleSection(bad, 0.2);
    } catch (const std::invalid_argument&) {
      rectangle_threw = true;
    }
    require(rectangle_threw, "invalid rectangle dimension should throw");

    bool circle_threw = false;
    try {
      (void)fem::CircularSection(bad);
    } catch (const std::invalid_argument&) {
      circle_threw = true;
    }
    require(circle_threw, "invalid circle radius should throw");
  }
}

void testGeneralSectionAxialAndBiaxialNormalStress() {
  const fem::GeneralSection section(0.02, 8.0e-5, 1.2e-4, 6.0e-5);
  fem::BeamSectionForces forces;
  forces.N = 100000.0;
  forces.My = 12000.0;
  forces.Mz = -8000.0;

  constexpr double y = 0.06;
  constexpr double z = -0.04;
  const double expected =
      forces.N / section.area -
      forces.Mz * y / section.iz +
      forces.My * z / section.iy;

  const double sigma =
      fem::BeamSectionStressRecovery{}.normalStressAt(section, forces, y, z);
  near(sigma, expected, 1.0e-12, "general section normal stress");
}

void testRectangleFullStress() {
  const fem::RectangleSection section(0.20, 0.40);
  fem::BeamSectionForces forces;
  forces.N = 80000.0;
  forces.My = 6000.0;
  forces.Mz = -4000.0;
  forces.Vy = 12000.0;
  forces.Vz = -6000.0;

  constexpr double y = 0.0;
  constexpr double z = 0.0;
  const auto stress =
      fem::BeamSectionStressRecovery{}.stressAt(section, forces, y, z);

  near(stress.sigma_x, forces.N / section.area, 1.0e-12,
       "rectangle center normal stress");
  near(stress.tau_xy, 1.5 * forces.Vy / section.area, 1.0e-12,
       "rectangle center tau_xy");
  near(stress.tau_xz, 1.5 * forces.Vz / section.area, 1.0e-12,
       "rectangle center tau_xz");

  const auto edge_y =
      fem::BeamSectionStressRecovery{}.stressAt(
          section, forces, section.size_y / 2.0, 0.0);
  near(edge_y.tau_xy, 0.0, 1.0e-10,
       "rectangle y-edge transverse shear should vanish");

  const auto edge_z =
      fem::BeamSectionStressRecovery{}.stressAt(
          section, forces, 0.0, section.size_z / 2.0);
  near(edge_z.tau_xz, 0.0, 1.0e-10,
       "rectangle z-edge transverse shear should vanish");
}

void testCircleTorsionStress() {
  const fem::CircularSection section(0.10);
  fem::BeamSectionForces forces;
  forces.T = 5000.0;

  constexpr double y = 0.06;
  constexpr double z = 0.08;
  const auto stress =
      fem::BeamSectionStressRecovery{}.stressAt(section, forces, y, z);

  near(stress.sigma_x, 0.0, 1.0e-12, "circle torsion sigma");
  near(stress.tau_xy, -forces.T * z / section.torsion_constant, 1.0e-12,
       "circle torsion tau_xy");
  near(stress.tau_xz, forces.T * y / section.torsion_constant, 1.0e-12,
       "circle torsion tau_xz");
}

void testPointOutsideKnownSectionRejected() {
  const fem::RectangleSection rectangle(0.20, 0.40);
  const fem::CircularSection circle(0.10);
  fem::BeamSectionForces forces;

  bool rectangle_threw = false;
  try {
    (void)fem::BeamSectionStressRecovery{}.normalStressAt(
        rectangle, forces, 0.11, 0.0);
  } catch (const std::invalid_argument&) {
    rectangle_threw = true;
  }
  require(rectangle_threw, "rectangle point outside section should throw");

  bool circle_threw = false;
  try {
    (void)fem::BeamSectionStressRecovery{}.normalStressAt(
        circle, forces, 0.08, 0.08);
  } catch (const std::invalid_argument&) {
    circle_threw = true;
  }
  require(circle_threw, "circle point outside section should throw");
}

void testUnsupportedShearTorsionCombinationsFailExplicitly() {
  fem::BeamSectionStressRecovery recovery;

  {
    const fem::GeneralSection section(0.02, 8.0e-5, 1.2e-4, 6.0e-5);
    fem::BeamSectionForces forces;
    forces.Vy = 1.0;
    bool threw = false;
    try {
      (void)recovery.stressAt(section, forces, 0.0, 0.0);
    } catch (const std::logic_error&) {
      threw = true;
    }
    require(threw, "general-section shear stress should fail explicitly");
  }

  {
    const fem::RectangleSection section(0.20, 0.40);
    fem::BeamSectionForces forces;
    forces.T = 1.0;
    bool threw = false;
    try {
      (void)recovery.stressAt(section, forces, 0.0, 0.0);
    } catch (const std::logic_error&) {
      threw = true;
    }
    require(threw, "rectangle torsion stress should fail explicitly");
  }

  {
    const fem::CircularSection section(0.10);
    fem::BeamSectionForces forces;
    forces.Vz = 1.0;
    bool threw = false;
    try {
      (void)recovery.stressAt(section, forces, 0.0, 0.0);
    } catch (const std::logic_error&) {
      threw = true;
    }
    require(threw, "circle transverse shear stress should fail explicitly");
  }
}

void testRectangleNormalStressExtrema() {
  const fem::RectangleSection section(0.20, 0.40);
  fem::BeamSectionForces forces;
  forces.N = 60000.0;
  forces.My = 5000.0;
  forces.Mz = -7000.0;

  const auto extrema =
      fem::BeamSectionStressRecovery{}.normalStressExtrema(section, forces);

  double expected_min = 0.0;
  double expected_max = 0.0;
  bool first = true;
  for (const double y : {-section.size_y / 2.0, section.size_y / 2.0}) {
    for (const double z : {-section.size_z / 2.0, section.size_z / 2.0}) {
      const double sigma =
          forces.N / section.area -
          forces.Mz * y / section.iz +
          forces.My * z / section.iy;
      if (first) {
        expected_min = sigma;
        expected_max = sigma;
        first = false;
      } else {
        expected_min = std::min(expected_min, sigma);
        expected_max = std::max(expected_max, sigma);
      }
    }
  }

  near(extrema.sigma_min, expected_min, 1.0e-12,
       "rectangle minimum normal stress");
  near(extrema.sigma_max, expected_max, 1.0e-12,
       "rectangle maximum normal stress");
}

void testCircleNormalStressExtrema() {
  const fem::CircularSection section(0.10);
  fem::BeamSectionForces forces;
  forces.N = 30000.0;
  forces.My = 4000.0;
  forces.Mz = -3000.0;

  const double base = forces.N / section.area;
  const double gradient_y = -forces.Mz / section.iz;
  const double gradient_z = forces.My / section.iy;
  const double delta =
      section.radius *
      std::sqrt(gradient_y * gradient_y + gradient_z * gradient_z);

  const auto extrema =
      fem::BeamSectionStressRecovery{}.normalStressExtrema(section, forces);

  near(extrema.sigma_min, base - delta, 1.0e-12,
       "circle minimum normal stress");
  near(extrema.sigma_max, base + delta, 1.0e-12,
       "circle maximum normal stress");
}

void testModelNormalStressRemainsAvailableWhenFullStressIsUnsupported() {
  constexpr double length = 2.0;
  constexpr double axial = 40000.0;
  constexpr double shear_y = -5000.0;
  const fem::GeneralSection section(0.02, 8.0e-5, 1.2e-4, 6.0e-5);

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);
  model.addElement<fem::Beam3D>(
      1,
      1,
      2,
      fem::LinearElasticMaterial(210.0e9, 0.3, 7850.0),
      section);

  tip.addLoad(fem::Dof::UX, axial);
  tip.addLoad(fem::Dof::UY, shear_y);

  const auto result = fem::LinearStaticSolver{}.solve(model);

  const double sigma =
      model.beamNormalStressAt(
          1,
          0.0,
          0.0,
          0.0,
          result.displacement);

  near(sigma, axial / section.area, 1.0e-10,
       "model normal stress should remain available for GeneralSection");

  bool full_stress_threw = false;
  try {
    (void)model.beamStressAt(
        1, 0.0, 0.0, 0.0, result.displacement);
  } catch (const std::logic_error&) {
    full_stress_threw = true;
  }
  require(full_stress_threw,
          "full GeneralSection stress should still reject transverse shear");
}

void testModelBeamStressRecoveryUsesElementSection() {
  constexpr double length = 3.0;
  constexpr double axial = 50000.0;
  constexpr double moment_z = 6000.0;
  const fem::RectangleSection section(0.20, 0.40);

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);
  model.addElement<fem::Beam3D>(
      1,
      1,
      2,
      fem::LinearElasticMaterial(210.0e9, 0.3, 7850.0),
      section);

  tip.addLoad(fem::Dof::UX, axial);
  tip.addLoad(fem::Dof::RZ, moment_z);

  const auto result = fem::LinearStaticSolver{}.solve(model);

  const double y = section.size_y / 2.0;
  const double expected =
      axial / section.area - moment_z * y / section.iz;

  const auto stress =
      model.beamStressAt(1, 0.0, y, 0.0, result.displacement);
  near(stress.sigma_x, expected, 1.0e-10,
       "model-integrated beam normal stress");

  const auto extrema =
      model.beamNormalStressExtrema(1, 0.0, result.displacement);
  require(extrema.sigma_min <= stress.sigma_x + 1.0e-10,
          "model stress should lie inside extrema");
  require(extrema.sigma_max >= stress.sigma_x - 1.0e-10,
          "model stress should lie inside extrema");
}

}  // namespace

int main() {
  testRectangleSectionProperties();
  testCircularSectionProperties();
  testInvalidGeometryRejected();
  testGeneralSectionAxialAndBiaxialNormalStress();
  testRectangleFullStress();
  testCircleTorsionStress();
  testPointOutsideKnownSectionRejected();
  testUnsupportedShearTorsionCombinationsFailExplicitly();
  testRectangleNormalStressExtrema();
  testCircleNormalStressExtrema();
  testModelNormalStressRemainsAvailableWhenFullStressIsUnsupported();
  testModelBeamStressRecoveryUsesElementSection();
  std::cout << "All section/stress recovery tests passed.\n";
  return EXIT_SUCCESS;
}
