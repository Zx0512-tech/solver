#include "fem/solver/LinearStaticSolver.hpp"
#include "fem/solver/SparseDofReducer.hpp"

#include <Eigen/SparseCholesky>

#include <stdexcept>
#include <vector>

namespace fem {

StaticResult LinearStaticSolver::solve(const Model& model) const {
  AssembledSystem system = model.assemble();
  const Eigen::Index global_size = system.load.size();

  std::vector<Eigen::Index> free;
  free.reserve(static_cast<std::size_t>(global_size));

  for (const NodeId node_id : system.dofs.nodeOrder()) {
    const Node& node = model.node(node_id);
    for (std::size_t offset = 0; offset < kDofsPerFrameNode; ++offset) {
      const Dof dof = dofFromOffset(offset);
      if (!node.isFixed(dof)) {
        free.push_back(
            static_cast<Eigen::Index>(
                system.dofs.equation(node_id, dof)));
      }
    }
  }

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
