#include "fem/solver/LinearStaticSolver.hpp"
#include "fem/solver/SparseDofReducer.hpp"

#include <Eigen/SparseCholesky>

#include <stdexcept>

namespace fem {

StaticResult LinearStaticSolver::solve(const Model& model) const {
  AssembledSystem system = model.assemble();
  const Eigen::Index global_size = system.load.size();

  const auto& free = system.dofs.freeEquations();

  if (free.empty()) {
    throw std::runtime_error(
        "Model has no free degrees of freedom");
  }

  const Eigen::SparseMatrix<double> kff =
      SparseDofReducer::matrix(system.stiffness, free);
  const Eigen::VectorXd ff =
      SparseDofReducer::vector(system.load, free);

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> factor;
  factor.compute(kff);
  if (factor.info() != Eigen::Success) {
    throw std::runtime_error(
        "Reduced stiffness matrix factorization failed; "
        "check constraints, connectivity and section properties");
  }

  const Eigen::VectorXd uf = factor.solve(ff);
  if (factor.info() != Eigen::Success || !uf.allFinite()) {
    throw std::runtime_error(
        "Reduced stiffness solve failed; "
        "check constraints, connectivity and section properties");
  }

  Eigen::VectorXd displacement =
      Eigen::VectorXd::Zero(global_size);
  for (std::size_t i = 0; i < free.size(); ++i) {
    displacement[free[i]] =
        uf[static_cast<Eigen::Index>(i)];
  }

  const Eigen::VectorXd reaction =
      system.stiffness * displacement - system.load;

  return {
      std::move(displacement),
      reaction,
      std::move(system.dofs)};
}

}  // namespace fem
