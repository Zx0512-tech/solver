#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/solver/LinearStaticSolver.hpp"

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
  auto& tip = model.addNode(2, {3.0, 0.0, 0.0});
  fixAll(root);

  const fem::RectangleSection section(0.20, 0.40);

  model.addElement<fem::Beam3D>(
      1,
      1,
      2,
      fem::LinearElasticMaterial(210.0e9, 0.3, 7850.0),
      section);

  tip.addLoad(fem::Dof::UX, 50000.0);
  tip.addLoad(fem::Dof::UY, -12000.0);
  tip.addLoad(fem::Dof::RZ, 6000.0);

  const auto result = fem::LinearStaticSolver{}.solve(model);

  // Full point stress is available for RectangleSection when torsion is zero.
  const auto center =
      model.beamStressAt(
          1,
          0.0,
          0.0,
          0.0,
          result.displacement);

  // Normal-stress extrema are recovered from the actual section boundary.
  const auto extrema =
      model.beamNormalStressExtrema(
          1,
          0.0,
          result.displacement);

  std::cout << "Root center sigma_x = " << center.sigma_x << " Pa\n";
  std::cout << "Root center tau_xy = " << center.tau_xy << " Pa\n";
  std::cout << "Root center tau_xz = " << center.tau_xz << " Pa\n";
  std::cout << "Root sigma_min = " << extrema.sigma_min << " Pa"
            << " at (y,z)=(" << extrema.min_point.y
            << "," << extrema.min_point.z << ")\n";
  std::cout << "Root sigma_max = " << extrema.sigma_max << " Pa"
            << " at (y,z)=(" << extrema.max_point.y
            << "," << extrema.max_point.z << ")\n";

  return 0;
}
