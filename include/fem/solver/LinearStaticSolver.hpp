#pragma once
#include "fem/core/Model.hpp"
#include <Eigen/Core>
namespace fem {
struct StaticResult {
  Eigen::VectorXd displacement;
  Eigen::VectorXd reaction;
  DofManager dofs;
  double displacementAt(NodeId node_id,Dof dof) const { return displacement[static_cast<Eigen::Index>(dofs.equation(node_id,dof))]; }
  double reactionAt(NodeId node_id,Dof dof) const { return reaction[static_cast<Eigen::Index>(dofs.equation(node_id,dof))]; }
};
class LinearStaticSolver { public: StaticResult solve(const Model& model) const; };
}  // namespace fem
