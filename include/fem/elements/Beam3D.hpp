#pragma once
#include "fem/core/Element.hpp"
#include "fem/core/Material.hpp"
#include "fem/core/Section.hpp"
#include <Eigen/Core>
namespace fem {
class Beam3D final : public Element {
 public:
  using Matrix12d=Eigen::Matrix<double,12,12>;
  Beam3D(ElementId id,NodeId node_i,NodeId node_j,LinearElasticMaterial material,BeamSection section,
         const Eigen::Vector3d& local_y_hint=Eigen::Vector3d::UnitY());
  NodeId nodeI() const noexcept { return node_i_; }
  NodeId nodeJ() const noexcept { return node_j_; }
  std::vector<ElementDof> dofs() const override;
  Eigen::MatrixXd stiffness(const NodeResolver& node) const override;
  double length(const NodeResolver& node) const;
  Matrix12d localStiffness(double length) const;
  Matrix12d transformation(const NodeResolver& node) const;
  Matrix12d globalStiffness(const NodeResolver& node) const;
 private:
  Eigen::Matrix3d rotationToLocal(const NodeResolver& node) const;
  NodeId node_i_, node_j_;
  LinearElasticMaterial material_;
  BeamSection section_;
  Eigen::Vector3d local_y_hint_;
};
}  // namespace fem
