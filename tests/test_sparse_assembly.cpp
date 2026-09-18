#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <type_traits>

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

void testAssemblyUsesSparseMatrices() {
  using StiffnessType =
      std::decay_t<decltype(std::declval<fem::AssembledSystem>().stiffness)>;
  using MassType =
      std::decay_t<decltype(std::declval<fem::AssembledSystem>().mass)>;

  static_assert(
      std::is_same_v<StiffnessType, Eigen::SparseMatrix<double>>,
      "AssembledSystem::stiffness must be Eigen::SparseMatrix<double>");
  static_assert(
      std::is_same_v<MassType, Eigen::SparseMatrix<double>>,
      "AssembledSystem::mass must be Eigen::SparseMatrix<double>");
}

void testTwoElementAssemblyAccumulatesSharedDofs() {
  constexpr double e = 210.0e9;
  constexpr double rho = 7850.0;
  constexpr double area = 0.01;
  constexpr double length = 2.0;

  fem::Model model;
  model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {length, 0.0, 0.0});
  model.addNode(3, {2.0 * length, 0.0, 0.0});

  const fem::LinearElasticMaterial material(e, 0.3, rho);
  const fem::BeamSection section(area, 2.0e-5, 3.0e-5, 4.0e-5);
  model.addElement<fem::Beam3D>(1, 1, 2, material, section);
  model.addElement<fem::Beam3D>(2, 2, 3, material, section);

  const auto system = model.assemble();

  require(system.stiffness.rows() == 18 && system.stiffness.cols() == 18,
          "sparse global stiffness size");
  require(system.mass.rows() == 18 && system.mass.cols() == 18,
          "sparse global mass size");

  const Eigen::Index ux1 =
      static_cast<Eigen::Index>(system.dofs.equation(1, fem::Dof::UX));
  const Eigen::Index ux2 =
      static_cast<Eigen::Index>(system.dofs.equation(2, fem::Dof::UX));
  const Eigen::Index ux3 =
      static_cast<Eigen::Index>(system.dofs.equation(3, fem::Dof::UX));

  const double axial_k = e * area / length;
  near(system.stiffness.coeff(ux1, ux1), axial_k, 1e-12,
       "node 1 axial stiffness diagonal");
  near(system.stiffness.coeff(ux2, ux2), 2.0 * axial_k, 1e-12,
       "shared node axial stiffness accumulates two elements");
  near(system.stiffness.coeff(ux3, ux3), axial_k, 1e-12,
       "node 3 axial stiffness diagonal");
  near(system.stiffness.coeff(ux1, ux2), -axial_k, 1e-12,
       "element 1 axial coupling");
  near(system.stiffness.coeff(ux2, ux3), -axial_k, 1e-12,
       "element 2 axial coupling");

  const double axial_mass_diag = rho * area * length / 3.0;
  const double axial_mass_offdiag = rho * area * length / 6.0;
  near(system.mass.coeff(ux1, ux1), axial_mass_diag, 1e-12,
       "node 1 axial consistent mass diagonal");
  near(system.mass.coeff(ux2, ux2), 2.0 * axial_mass_diag, 1e-12,
       "shared node consistent mass accumulates");
  near(system.mass.coeff(ux1, ux2), axial_mass_offdiag, 1e-12,
       "element 1 axial mass coupling");
  near(system.mass.coeff(ux2, ux3), axial_mass_offdiag, 1e-12,
       "element 2 axial mass coupling");
}

void testSparseAssemblyPreservesSymmetryAndSparsity() {
  fem::Model model;
  model.addNode(1, {0.0, 0.0, 0.0});
  model.addNode(2, {2.0, 0.0, 0.0});
  model.addNode(3, {4.0, 0.0, 0.0});
  model.addNode(4, {6.0, 0.0, 0.0});

  const fem::LinearElasticMaterial material(200.0e9, 0.3, 7850.0);
  const fem::BeamSection section(0.02, 3.0e-5, 5.0e-5, 7.0e-5);
  model.addElement<fem::Beam3D>(1, 1, 2, material, section);
  model.addElement<fem::Beam3D>(2, 2, 3, material, section);
  model.addElement<fem::Beam3D>(3, 3, 4, material, section);

  const auto system = model.assemble();
  const Eigen::MatrixXd k = Eigen::MatrixXd(system.stiffness);
  const Eigen::MatrixXd m = Eigen::MatrixXd(system.mass);

  near((k - k.transpose()).cwiseAbs().maxCoeff(), 0.0, 1e-10,
       "sparse stiffness symmetry");
  near((m - m.transpose()).cwiseAbs().maxCoeff(), 0.0, 1e-10,
       "sparse mass symmetry");

  const Eigen::Index dense_entries =
      system.stiffness.rows() * system.stiffness.cols();
  require(system.stiffness.nonZeros() < dense_entries,
          "global stiffness should not store a dense matrix");
  require(system.mass.nonZeros() < dense_entries,
          "global mass should not store a dense matrix");
}

}  // namespace

int main() {
  testAssemblyUsesSparseMatrices();
  testTwoElementAssemblyAccumulatesSharedDofs();
  testSparseAssemblyPreservesSymmetryAndSparsity();
  std::cout << "All sparse-assembly tests passed.\n";
  return EXIT_SUCCESS;
}
