#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/dynamics/RayleighDamping.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/recorder/EnvelopeRecorder.hpp"
#include "fem/recorder/NodeRecorder.hpp"
#include "fem/recorder/SectionRecorder.hpp"
#include "fem/solver/ModalSolver.hpp"
#include "fem/solver/NewmarkBetaSolver.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace {
void fixAll(fem::Node& node) {
  for (std::size_t i = 0; i < fem::kDofsPerFrameNode; ++i) {
    node.fix(fem::dofFromOffset(i));
  }
}

void printEnvelope(const char* name, const fem::EnvelopeResult& e) {
  std::cout << name
            << ": min=" << e.minimum.value
            << " at t=" << e.minimum.time
            << ", max=" << e.maximum.value
            << " at t=" << e.maximum.time
            << ", max|.|=" << e.maximum_absolute.magnitude
            << " (value=" << e.maximum_absolute.value
            << ") at t=" << e.maximum_absolute.time
            << "\n";
}
}  // namespace

int main() {
  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {2.0, 0.0, 0.0});
  fixAll(root);

  const fem::LinearElasticMaterial steel(210.0e9, 0.30, 7850.0);
  const fem::RectangleSection section(0.10, 0.20);
  model.addElement<fem::Beam3D>(1, 1, 2, steel, section);

  const auto modes = fem::ModalSolver{}.solve(model, 2);
  const auto damping = fem::RayleighDamping::fromModalTargets(
      modes.angular_frequencies[0], 0.02,
      modes.angular_frequencies[1], 0.02);

  fem::NewmarkSettings settings;
  settings.time_step = 0.001;
  settings.step_count = 1000;

  const double pulse_duration = 0.05;
  const std::vector<fem::NodalTimeLoad> loads = {
      {2, fem::Dof::UY, [pulse_duration](double t) {
         if (t < 0.0 || t > pulse_duration) {
           return 0.0;
         }
         const double pi = std::acos(-1.0);
         return -10000.0 * std::sin(pi * t / pulse_duration);
       }}
  };

  const auto result =
      fem::NewmarkBetaSolver{}.solve(model, settings, damping, loads);

  const auto node_history =
      fem::NodeRecorder{}.record(2, result);
  const auto root_section =
      fem::SectionRecorder{}.record(model, 1, 0.0, result);
  const auto root_stress =
      fem::SectionRecorder{}.recordNormalStress(
          model,
          1,
          0.0,
          section.size_y / 2.0,
          0.0,
          result);

  const fem::EnvelopeRecorder envelope;
  printEnvelope(
      "tip Uy",
      envelope.record(
          node_history,
          fem::NodeResponseQuantity::Displacement,
          fem::Dof::UY));
  printEnvelope(
      "root Mz",
      envelope.record(
          root_section,
          fem::BeamSectionForceComponent::Mz));
  printEnvelope(
      "root edge sigma_x",
      envelope.record(root_stress));

  return 0;
}
