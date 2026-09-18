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

#include <filesystem>
#include <iostream>

namespace {

void fixAll(fem::Node& node) {
  for (std::size_t i = 0; i < fem::kDofsPerFrameNode; ++i) {
    node.fix(fem::dofFromOffset(i));
  }
}

}  // namespace

int main() {
  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {6.0, 0.0, 0.0});
  fixAll(root);

  model.addElement<fem::Beam3D>(
      1,
      1,
      2,
      fem::LinearElasticMaterial(210.0e9, 0.3, 7850.0),
      fem::BeamSection(0.02, 8.0e-5, 1.2e-4, 6.0e-5));

  model.addElementLoad<fem::BeamPointLoad3D>(
      1,
      2.25,
      Eigen::Vector3d(0.0, -18000.0, 0.0));

  model.addElementLoad<fem::BeamPartialLinearLoad3D>(
      1,
      3.0,
      5.5,
      Eigen::Vector3d(0.0, -1500.0, 0.0),
      Eigen::Vector3d(0.0, -3500.0, 0.0));

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const auto samples =
      fem::BeamResponseSampler{}.sample(model, 1, result.displacement, 31);

  const std::filesystem::path output_dir = "beam_force_output";
  fem::BeamForceCsvWriter{}.writeFile(
      output_dir / "beam_1_forces.csv", samples);
  fem::BeamForceDiagramSvgWriter{}.writeSet(
      output_dir, "beam_1", samples);

  std::cout << "Wrote Beam3D force results to "
            << output_dir.string() << "\n";
  return 0;
}
