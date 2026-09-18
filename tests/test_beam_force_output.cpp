#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/io/BeamForceCsvWriter.hpp"
#include "fem/io/BeamForceDiagramSvgWriter.hpp"
#include "fem/loads/BeamPartialLinearLoad3D.hpp"
#include "fem/loads/BeamPointLoad3D.hpp"
#include "fem/response/BeamResponseSampler.hpp"
#include "fem/solver/LinearStaticSolver.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
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

fem::BeamSection section() {
  return fem::BeamSection(0.01, 6.0e-6, 8.0e-6, 1.0e-5);
}

fem::LinearElasticMaterial material() {
  return fem::LinearElasticMaterial(200.0e9, 0.3, 7850.0);
}

void testUniformSamplingAndForces() {
  constexpr double l = 4.0;
  constexpr double p = -10000.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {l, 0.0, 0.0});
  fixAll(root);
  model.addElement<fem::Beam3D>(1, 1, 2, material(), section());
  tip.addLoad(fem::Dof::UY, p);

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto samples =
      fem::BeamResponseSampler{}.sample(model, 1, result.displacement, 5);

  require(samples.size() == 5, "uniform sampler should return five points");
  for (std::size_t i = 0; i < samples.size(); ++i) {
    const double x = static_cast<double>(i);
    near(samples[i].forces.x, x, 1.0e-12, "uniform sample x");
    near(samples[i].forces.Vy, p, 1.0e-10, "uniform sample Vy");
    near(samples[i].forces.Mz, p * (l - x), 1.0e-10, "uniform sample Mz");
    require(samples[i].side == fem::BeamSectionSide::Right,
            "regular samples use right-side convention");
  }
}

void testPointLoadAddsLeftAndRightSamples() {
  constexpr double l = 4.0;
  constexpr double a = 1.3;
  constexpr double p = -9000.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {l, 0.0, 0.0});
  fixAll(root);
  model.addElement<fem::Beam3D>(1, 1, 2, material(), section());
  model.addElementLoad<fem::BeamPointLoad3D>(
      1, a, Eigen::Vector3d(0.0, p, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto samples =
      fem::BeamResponseSampler{}.sample(model, 1, result.displacement, 5);

  require(samples.size() == 7,
          "point-load location should add left and right samples");

  std::vector<std::size_t> at_load;
  for (std::size_t i = 0; i < samples.size(); ++i) {
    if (std::abs(samples[i].forces.x - a) < 1.0e-12) {
      at_load.push_back(i);
    }
  }
  require(at_load.size() == 2, "point load should have two samples at same x");

  const auto& left = samples[at_load[0]];
  const auto& right = samples[at_load[1]];
  require(left.side == fem::BeamSectionSide::Left,
          "first discontinuity sample should be left limit");
  require(right.side == fem::BeamSectionSide::Right,
          "second discontinuity sample should be right limit");
  near(left.forces.Vy, p, 1.0e-10, "left-limit shear");
  near(right.forces.Vy, 0.0, 1.0e-8, "right-limit shear");
}

void testTipPointLoadPreservesLeftLimitAtBeamEnd() {
  constexpr double l = 4.0;
  constexpr double p = -9000.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {l, 0.0, 0.0});
  fixAll(root);
  model.addElement<fem::Beam3D>(1, 1, 2, material(), section());
  model.addElementLoad<fem::BeamPointLoad3D>(
      1, l, Eigen::Vector3d(0.0, p, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto samples =
      fem::BeamResponseSampler{}.sample(model, 1, result.displacement, 3);

  std::vector<const fem::BeamForceSample*> at_tip;
  for (const auto& sample : samples) {
    if (std::abs(sample.forces.x - l) < 1.0e-12) {
      at_tip.push_back(&sample);
    }
  }

  require(at_tip.size() == 2,
          "tip point load should preserve left and right limits");
  require(at_tip[0]->side == fem::BeamSectionSide::Left,
          "tip point load first sample should be left limit");
  require(at_tip[1]->side == fem::BeamSectionSide::Right,
          "tip point load second sample should be right limit");
  near(at_tip[0]->forces.Vy, p, 1.0e-10,
       "tip point-load left-limit shear");
  near(at_tip[1]->forces.Vy, 0.0, 1.0e-8,
       "tip point-load right-limit shear");
}

void testPartialLoadBoundariesAreInsertedOnce() {
  constexpr double l = 4.0;
  constexpr double a = 0.7;
  constexpr double b = 3.2;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {l, 0.0, 0.0});
  fixAll(root);
  model.addElement<fem::Beam3D>(1, 1, 2, material(), section());
  model.addElementLoad<fem::BeamPartialLinearLoad3D>(
      1, a, b,
      Eigen::Vector3d(0.0, -1000.0, 0.0),
      Eigen::Vector3d(0.0, -2000.0, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto samples =
      fem::BeamResponseSampler{}.sample(model, 1, result.displacement, 3);

  // Uniform points: 0, 2, 4. Partial-load breakpoints: 0.7, 3.2.
  require(samples.size() == 5,
          "partial-load boundaries should be inserted without duplication");

  const auto count_x = [&samples](double x) {
    return std::count_if(
        samples.begin(), samples.end(),
        [x](const fem::BeamForceSample& sample) {
          return std::abs(sample.forces.x - x) < 1.0e-12;
        });
  };
  require(count_x(a) == 1, "partial-load start should appear once");
  require(count_x(b) == 1, "partial-load end should appear once");
}

void testSamplerRejectsTooFewPoints() {
  fem::Model model;
  model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {2.0, 0.0, 0.0});
  model.addElement<fem::Beam3D>(1, 1, 2, material(), section());
  const auto assembled = model.assemble();
  const Eigen::VectorXd zero =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(assembled.dofs.size()));

  bool threw = false;
  try {
    (void)fem::BeamResponseSampler{}.sample(model, 1, zero, 1);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  require(threw, "sampler should reject point_count < 2");
}

void testCsvWriterIncludesSideAndAllForceComponents() {
  std::vector<fem::BeamForceSample> samples = {
      {{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, fem::BeamSectionSide::Right},
      {{1.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0}, fem::BeamSectionSide::Left},
  };

  std::ostringstream stream;
  fem::BeamForceCsvWriter{}.write(stream, samples);
  const std::string csv = stream.str();

  require(csv.find("x,side,N,Vy,Vz,T,My,Mz\n") == 0,
          "CSV header should be stable");
  require(csv.find("0,right,1,2,3,4,5,6") != std::string::npos,
          "CSV should include right-side row");
  require(csv.find("1,left,7,8,9,10,11,12") != std::string::npos,
          "CSV should include left-side row");
}

void testSvgWritersProduceFourDiagramTypes() {
  std::vector<fem::BeamForceSample> samples = {
      {{0.0, 0.0, -10.0, 5.0, 3.0, 2.0, -20.0}, fem::BeamSectionSide::Right},
      {{1.0, 0.0, -5.0, 2.0, 2.0, 1.0, -10.0}, fem::BeamSectionSide::Right},
      {{2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, fem::BeamSectionSide::Right},
  };

  fem::BeamForceDiagramSvgWriter writer;

  std::ostringstream axial;
  writer.writeAxial(axial, samples);
  require(axial.str().find("<svg") != std::string::npos,
          "axial diagram should be SVG");
  require(axial.str().find("Axial force N") != std::string::npos,
          "axial diagram title");

  std::ostringstream shear;
  writer.writeShear(shear, samples);
  require(shear.str().find("Shear forces Vy / Vz") != std::string::npos,
          "shear diagram title");
  require(shear.str().find("Vy") != std::string::npos &&
              shear.str().find("Vz") != std::string::npos,
          "shear diagram should label both components");

  std::ostringstream bending;
  writer.writeBending(bending, samples);
  require(bending.str().find("Bending moments My / Mz") != std::string::npos,
          "bending diagram title");

  std::ostringstream torsion;
  writer.writeTorsion(torsion, samples);
  require(torsion.str().find("Torsion T") != std::string::npos,
          "torsion diagram title");
}

void testFileOutputCreatesCsvAndFourSvgFiles() {
  const auto temp =
      std::filesystem::temp_directory_path() / "solver_beam_force_output_test";
  std::filesystem::remove_all(temp);
  std::filesystem::create_directories(temp);

  const std::vector<fem::BeamForceSample> samples = {
      {{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, fem::BeamSectionSide::Right},
      {{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, fem::BeamSectionSide::Right},
  };

  fem::BeamForceCsvWriter{}.writeFile(temp / "beam_1_forces.csv", samples);
  fem::BeamForceDiagramSvgWriter{}.writeSet(temp, "beam_1", samples);

  for (const auto& filename : {
           "beam_1_forces.csv",
           "beam_1_axial.svg",
           "beam_1_shear.svg",
           "beam_1_bending.svg",
           "beam_1_torsion.svg"}) {
    const auto path = temp / filename;
    require(std::filesystem::exists(path), "expected output file missing");
    require(std::filesystem::file_size(path) > 0, "output file should not be empty");
  }

  std::filesystem::remove_all(temp);
}

}  // namespace

int main() {
  testUniformSamplingAndForces();
  testPointLoadAddsLeftAndRightSamples();
  testTipPointLoadPreservesLeftLimitAtBeamEnd();
  testPartialLoadBoundariesAreInsertedOnce();
  testSamplerRejectsTooFewPoints();
  testCsvWriterIncludesSideAndAllForceComponents();
  testSvgWritersProduceFourDiagramTypes();
  testFileOutputCreatesCsvAndFourSvgFiles();
  std::cout << "All Beam3D force sampling/output tests passed.\n";
  return EXIT_SUCCESS;
}
