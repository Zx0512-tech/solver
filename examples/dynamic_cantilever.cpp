#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/dynamics/RayleighDamping.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/solver/ModalSolver.hpp"
#include "fem/solver/NewmarkBetaSolver.hpp"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>

namespace {
void fixAll(fem::Node& node) {
  for (std::size_t i = 0; i < fem::kDofsPerFrameNode; ++i) {
    node.fix(fem::dofFromOffset(i));
  }
}
}

int main() {
  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {2.0, 0.0, 0.0});
  fixAll(root);

  const fem::LinearElasticMaterial steel(210.0e9, 0.30, 7850.0);
  const fem::BeamSection section(8.0e-3, 6.0e-6, 8.0e-6, 1.0e-5);
  model.addElement<fem::Beam3D>(1, 1, 2, steel, section);

  const auto modes = fem::ModalSolver{}.solve(model, 2);
  const auto damping = fem::RayleighDamping::fromModalTargets(
      modes.angular_frequencies[0], 0.02,
      modes.angular_frequencies[1], 0.02);

  fem::NewmarkSettings settings;
  settings.time_step = 0.001;
  settings.step_count = 2000;

  const double pulse_duration = 0.05;
  const std::vector<fem::NodalTimeLoad> loads = {
      {2, fem::Dof::UY, [pulse_duration](double t) {
         if (t < 0.0 || t > pulse_duration) {
           return 0.0;
         }
         constexpr double pi = 3.14159265358979323846;
         return -10000.0 * std::sin(pi * t / pulse_duration);
       }}
  };

  const auto response = fem::NewmarkBetaSolver{}.solve(model, settings, damping, loads);

  std::cout << std::fixed << std::setprecision(4)
            << "f1 = " << modes.frequencies_hz[0] << " Hz\n"
            << "f2 = " << modes.frequencies_hz[1] << " Hz\n"
            << std::scientific << std::setprecision(6)
            << "final Uy = "
            << response.displacementAt(settings.step_count, 2, fem::Dof::UY)
            << " m\n";
}
