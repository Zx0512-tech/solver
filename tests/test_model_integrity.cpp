#include "fem/core/DofManager.hpp"
#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/solver/LinearStaticSolver.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

void near(double actual, double expected, double tol, const std::string& message) {
  const double scale = std::max(1.0, std::abs(expected));
  if (std::abs(actual - expected) > tol * scale) {
    std::cerr << "FAIL: " << message << " actual=" << actual
              << " expected=" << expected << '\n';
    std::exit(EXIT_FAILURE);
  }
}

bool contains(
    const std::vector<Eigen::Index>& values,
    Eigen::Index target) {
  return std::find(values.begin(), values.end(), target) != values.end();
}

void fixAll(fem::Node& node) {
  for (std::size_t i = 0; i < fem::kDofsPerFrameNode; ++i) {
    node.fix(fem::dofFromOffset(i));
  }
}

fem::LinearElasticMaterial material() {
  return fem::LinearElasticMaterial(210.0e9, 0.3, 7850.0);
}

fem::BeamSection section() {
  return fem::BeamSection(0.01, 2.0e-5, 3.0e-5, 4.0e-5);
}

void testDeterministicNumberingAndConstraintPartition() {
  fem::Model model;
  model.addNode(100, {2.0, 0.0, 0.0});
  model.addNode(7, {0.0, 0.0, 0.0});
  model.addNode(5000, {3.0, 0.0, 0.0});
  model.addNode(21, {1.0, 0.0, 0.0});

  model.node(21).fix(fem::Dof::UX);
  model.node(100).fix(fem::Dof::RZ);

  const fem::DofManager dm(model.nodes());
  const std::vector<fem::NodeId> expected_order = {7, 21, 100, 5000};
  require(dm.nodeOrder() == expected_order,
          "DofManager must number nodes deterministically by NodeId");

  require(dm.equation(7, fem::Dof::UX) == 0U,
          "node 7 base equation");
  require(dm.equation(21, fem::Dof::UX) == 6U,
          "node 21 base equation");
  require(dm.equation(100, fem::Dof::UX) == 12U,
          "node 100 base equation");
  require(dm.equation(5000, fem::Dof::RZ) == 23U,
          "node 5000 final equation");

  const auto& free = dm.freeEquations();
  const auto& constrained = dm.constrainedEquations();
  require(free.size() == 22U && constrained.size() == 2U,
          "free/constrained equation counts");
  require(contains(constrained, 6),
          "node 21 UX must be constrained");
  require(contains(constrained, 17),
          "node 100 RZ must be constrained");

  std::vector<unsigned char> seen(dm.size(), 0U);
  for (const Eigen::Index eq : free) {
    require(seen[static_cast<std::size_t>(eq)] == 0U,
            "free equations must be unique");
    seen[static_cast<std::size_t>(eq)] = 1U;
  }
  for (const Eigen::Index eq : constrained) {
    require(seen[static_cast<std::size_t>(eq)] == 0U,
            "free and constrained equations must be disjoint");
    seen[static_cast<std::size_t>(eq)] = 1U;
  }
  require(
      std::all_of(seen.begin(), seen.end(),
                  [](unsigned char value) { return value == 1U; }),
      "free and constrained equations must cover every global DOF");
}

void testElementEquationMappingPreservesInputOrder() {
  fem::Model model;
  model.addNode(100, {2.0, 0.0, 0.0});
  model.addNode(7, {0.0, 0.0, 0.0});
  model.addNode(5000, {3.0, 0.0, 0.0});
  model.addNode(21, {1.0, 0.0, 0.0});

  const fem::DofManager dm(model.nodes());
  const std::vector<fem::ElementDof> element_dofs = {
      {100, fem::Dof::RZ},
      {7, fem::Dof::UX},
      {100, fem::Dof::RZ},
      {21, fem::Dof::UY},
  };

  const auto equations = dm.equations(element_dofs);
  const std::vector<Eigen::Index> expected = {17, 0, 17, 7};
  require(equations == expected,
          "element equation mapping must preserve order and duplicates");
}

void testConstraintSnapshotSemantics() {
  fem::Model model;
  model.addNode(2, {1.0, 0.0, 0.0});
  model.addNode(1, {0.0, 0.0, 0.0});

  const auto before = model.assemble();
  const Eigen::Index ux1 = static_cast<Eigen::Index>(
      before.dofs.equation(1, fem::Dof::UX));
  require(contains(before.dofs.freeEquations(), ux1),
          "initial snapshot must see node 1 UX as free");

  model.node(1).fix(fem::Dof::UX);

  require(contains(before.dofs.freeEquations(), ux1),
          "existing DofManager snapshot must not mutate retroactively");

  const auto after = model.assemble();
  require(contains(after.dofs.constrainedEquations(), ux1),
          "fresh assembly must rebuild constraint partition");
}

void testBeamConnectivityMatchesDofDeclaration() {
  const fem::Beam3D beam(
      5, 30, 10, material(), section());

  const auto nodes = beam.nodeIds();
  require(nodes.size() == 2U && nodes[0] == 30 && nodes[1] == 10,
          "Beam3D nodeIds must preserve connectivity order");

  const auto dofs = beam.dofs();
  require(dofs.size() == 12U, "Beam3D must expose 12 element DOFs");
  for (std::size_t i = 0; i < 6U; ++i) {
    require(dofs[i].first == 30,
            "Beam3D first six DOFs must belong to node i");
  }
  for (std::size_t i = 6U; i < 12U; ++i) {
    require(dofs[i].first == 10,
            "Beam3D final six DOFs must belong to node j");
  }
}

void testElementInsertionValidatesConnectivityWithoutPartialMutation() {
  fem::Model model;
  model.addNode(1, {0.0, 0.0, 0.0});

  bool threw = false;
  std::string message;
  try {
    model.addElement<fem::Beam3D>(
        10, 1, 999, material(), section());
  } catch (const std::invalid_argument& error) {
    threw = true;
    message = error.what();
  }

  require(threw, "element insertion must reject missing connectivity");
  require(message.find("Element 10") != std::string::npos &&
              message.find("999") != std::string::npos,
          "connectivity error must include element and node ids");

  bool element_exists = true;
  try {
    (void)model.element(10);
  } catch (const std::out_of_range&) {
    element_exists = false;
  }
  require(!element_exists,
          "failed insertion must not leave a partial element in the model");

  model.addNode(999, {2.0, 0.0, 0.0});
  auto& beam = model.addElement<fem::Beam3D>(
      10, 1, 999, material(), section());
  require(beam.id() == 10,
          "element id must remain reusable after failed insertion");
}

void testOutOfOrderNodeIdsSolveCorrectly() {
  constexpr double e = 210.0e9;
  constexpr double iz = 3.0e-5;
  constexpr double length = 2.0;
  constexpr double load = -10000.0;

  fem::Model model;
  auto& root = model.addNode(100, {0.0, 0.0, 0.0});
  model.addNode(7, {length, 0.0, 0.0});
  fixAll(root);

  model.addElement<fem::Beam3D>(
      1,
      100,
      7,
      fem::LinearElasticMaterial(e, 0.3, 7850.0),
      fem::BeamSection(0.01, 2.0e-5, iz, 4.0e-5));
  model.node(7).addLoad(fem::Dof::UY, load);

  const auto result = fem::LinearStaticSolver{}.solve(model);
  const double expected =
      load * length * length * length / (3.0 * e * iz);

  near(
      result.displacementAt(7, fem::Dof::UY),
      expected,
      1.0e-10,
      "out-of-order/non-contiguous NodeIds must not change the solution");
  near(
      result.reactionAt(100, fem::Dof::UY),
      -load,
      1.0e-10,
      "out-of-order/non-contiguous NodeIds must preserve reactions");
}

}  // namespace

int main() {
  testDeterministicNumberingAndConstraintPartition();
  testElementEquationMappingPreservesInputOrder();
  testConstraintSnapshotSemantics();
  testBeamConnectivityMatchesDofDeclaration();
  testElementInsertionValidatesConnectivityWithoutPartialMutation();
  testOutOfOrderNodeIdsSolveCorrectly();
  std::cout << "All model-integrity tests passed.\n";
  return EXIT_SUCCESS;
}
