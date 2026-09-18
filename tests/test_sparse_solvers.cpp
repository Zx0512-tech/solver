#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/dynamics/RayleighDamping.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/solver/LinearStaticSolver.hpp"
#include "fem/solver/NewmarkBetaSolver.hpp"
#include "fem/solver/SparseDofReducer.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <type_traits>
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

void constrainToAxial(fem::Node& node) {
  node.fix(fem::Dof::UY);
  node.fix(fem::Dof::UZ);
  node.fix(fem::Dof::RX);
  node.fix(fem::Dof::RY);
  node.fix(fem::Dof::RZ);
}

void testSparseDofReducerKeepsSparseTypeAndValues() {
  using ReducedType = decltype(
      fem::SparseDofReducer::matrix(
          std::declval<const Eigen::SparseMatrix<double>&>(),
          std::declval<const std::vector<Eigen::Index>&>()));

  static_assert(
      std::is_same_v<ReducedType, Eigen::SparseMatrix<double>>,
      "SparseDofReducer::matrix must return Eigen::SparseMatrix<double>");

  Eigen::SparseMatrix<double> matrix(5, 5);
  std::vector<Eigen::Triplet<double>> triplets = {
      {0, 0, 10.0},
      {0, 2, -2.0},
      {1, 1, 20.0},
      {2, 0, -2.0},
      {2, 2, 30.0},
      {2, 4, -4.0},
      {4, 2, -4.0},
      {4, 4, 50.0},
  };
  matrix.setFromTriplets(triplets.begin(), triplets.end());

  const std::vector<Eigen::Index> keep = {0, 2, 4};
  const auto reduced = fem::SparseDofReducer::matrix(matrix, keep);

  require(reduced.rows() == 3 && reduced.cols() == 3,
          "reduced sparse matrix shape");
  near(reduced.coeff(0, 0), 10.0, 1e-12, "reduced (0,0)");
  near(reduced.coeff(0, 1), -2.0, 1e-12, "reduced (0,1)");
  near(reduced.coeff(1, 2), -4.0, 1e-12, "reduced (1,2)");
  near(reduced.coeff(2, 2), 50.0, 1e-12, "reduced (2,2)");
  require(reduced.nonZeros() == 7, "reduced matrix should preserve sparsity");

  const Eigen::VectorXd full =
      (Eigen::VectorXd(5) << 1.0, 2.0, 3.0, 4.0, 5.0).finished();
  const Eigen::VectorXd reduced_vector =
      fem::SparseDofReducer::vector(full, keep);
  near(reduced_vector[0], 1.0, 1e-12, "reduced vector first");
  near(reduced_vector[1], 3.0, 1e-12, "reduced vector second");
  near(reduced_vector[2], 5.0, 1e-12, "reduced vector third");
}

void testSparseRayleighDamping() {
  Eigen::SparseMatrix<double> mass(2, 2);
  Eigen::SparseMatrix<double> stiffness(2, 2);
  mass.insert(0, 0) = 2.0;
  mass.insert(1, 1) = 3.0;
  stiffness.insert(0, 0) = 100.0;
  stiffness.insert(0, 1) = -20.0;
  stiffness.insert(1, 0) = -20.0;
  stiffness.insert(1, 1) = 80.0;
  mass.makeCompressed();
  stiffness.makeCompressed();

  const fem::RayleighDamping damping(0.5, 0.02);
  const auto c = damping.matrix(mass, stiffness);

  static_assert(
      std::is_same_v<std::decay_t<decltype(c)>, Eigen::SparseMatrix<double>>,
      "Sparse Rayleigh damping must remain sparse");

  near(c.coeff(0, 0), 0.5 * 2.0 + 0.02 * 100.0, 1e-12,
       "sparse Rayleigh C00");
  near(c.coeff(0, 1), 0.02 * -20.0, 1e-12,
       "sparse Rayleigh C01");
  near(c.coeff(1, 1), 0.5 * 3.0 + 0.02 * 80.0, 1e-12,
       "sparse Rayleigh C11");
}

void testSparseStaticSolverOnAxialChain() {
  constexpr double e = 210.0e9;
  constexpr double area = 0.01;
  constexpr double element_length = 1.25;
  constexpr int element_count = 20;
  constexpr double load = 50000.0;

  fem::Model model;
  for (int i = 0; i <= element_count; ++i) {
    auto& node = model.addNode(
        static_cast<fem::NodeId>(i + 1),
        {static_cast<double>(i) * element_length, 0.0, 0.0});
    constrainToAxial(node);
  }
  model.node(1).fix(fem::Dof::UX);

  const fem::LinearElasticMaterial material(e, 0.3, 7850.0);
  const fem::BeamSection section(area, 2e-5, 3e-5, 4e-5);
  for (int i = 0; i < element_count; ++i) {
    model.addElement<fem::Beam3D>(
        static_cast<fem::ElementId>(i + 1),
        static_cast<fem::NodeId>(i + 1),
        static_cast<fem::NodeId>(i + 2),
        material,
        section);
  }

  model.node(static_cast<fem::NodeId>(element_count + 1))
      .addLoad(fem::Dof::UX, load);

  const auto result = fem::LinearStaticSolver{}.solve(model);

  const double total_length =
      static_cast<double>(element_count) * element_length;
  const double expected_tip = load * total_length / (e * area);

  near(
      result.displacementAt(
          static_cast<fem::NodeId>(element_count + 1),
          fem::Dof::UX),
      expected_tip,
      1e-10,
      "sparse static axial-chain tip displacement");
  near(
      result.reactionAt(1, fem::Dof::UX),
      -load,
      1e-10,
      "sparse static axial-chain root reaction");
}

void testSparseStaticSolverRejectsRigidBodyMechanism() {
  fem::Model model;
  model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {2.0, 0.0, 0.0});

  model.addElement<fem::Beam3D>(
      1,
      1,
      2,
      fem::LinearElasticMaterial(200.0e9, 0.3, 7850.0),
      fem::BeamSection(0.01, 2e-5, 3e-5, 4e-5));

  bool threw = false;
  try {
    (void)fem::LinearStaticSolver{}.solve(model);
  } catch (const std::runtime_error&) {
    threw = true;
  }

  require(
      threw,
      "sparse static solver must reject an unconstrained rigid-body mechanism");
}

void testSparseNewmarkMatchesClosedFormFirstStepForAxialSdoF() {
  constexpr double length = 2.0;
  constexpr double e = 200.0e9;
  constexpr double area = 0.01;
  constexpr double rho = 7800.0;

  fem::Model model;
  auto& root = model.addNode(1, {0.0, 0.0, 0.0});
  auto& tip = model.addNode(2, {length, 0.0, 0.0});
  fixAll(root);
  constrainToAxial(tip);

  model.addElement<fem::Beam3D>(
      1,
      1,
      2,
      fem::LinearElasticMaterial(e, 0.3, rho),
      fem::BeamSection(area, 2e-5, 3e-5, 4e-5));

  fem::NewmarkSettings settings;
  settings.time_step = 0.001;
  settings.step_count = 2;
  settings.beta = 0.25;
  settings.gamma = 0.5;

  constexpr double force = 1000.0;
  const std::vector<fem::NodalTimeLoad> loads = {
      {2, fem::Dof::UX,
       [](double t) { return t > 0.0 ? force : 0.0; }}};

  const auto result =
      fem::NewmarkBetaSolver{}.solve(model, settings, {}, loads);

  const double k = e * area / length;
  const double m = rho * area * length / 3.0;
  const double a0 =
      1.0 / (settings.beta * settings.time_step * settings.time_step);
  const double expected_u1 = force / (k + a0 * m);

  near(
      result.displacementAt(1, 2, fem::Dof::UX),
      expected_u1,
      1e-10,
      "sparse Newmark first-step displacement");
  require(
      std::isfinite(result.displacementAt(2, 2, fem::Dof::UX)),
      "sparse Newmark second-step displacement must be finite");
}

}  // namespace

int main() {
  testSparseDofReducerKeepsSparseTypeAndValues();
  testSparseRayleighDamping();
  testSparseStaticSolverOnAxialChain();
  testSparseStaticSolverRejectsRigidBodyMechanism();
  testSparseNewmarkMatchesClosedFormFirstStepForAxialSdoF();
  std::cout << "All sparse-solver tests passed.\n";
  return EXIT_SUCCESS;
}
