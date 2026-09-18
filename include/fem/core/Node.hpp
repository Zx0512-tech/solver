#pragma once
#include "fem/core/Types.hpp"
#include <Eigen/Core>
#include <array>
namespace fem {
class Node {
 public:
  Node(NodeId id, const Eigen::Vector3d& coordinates) : id_(id), coordinates_(coordinates) { fixed_.fill(false); }
  NodeId id() const noexcept { return id_; }
  const Eigen::Vector3d& coordinates() const noexcept { return coordinates_; }
  void fix(Dof dof, bool value=true) noexcept { fixed_[dofOffset(dof)] = value; }
  bool isFixed(Dof dof) const noexcept { return fixed_[dofOffset(dof)]; }
  void setLoad(Dof dof, double value) noexcept { loads_[dofOffset(dof)] = value; }
  void addLoad(Dof dof, double value) noexcept { loads_[dofOffset(dof)] += value; }
  double load(Dof dof) const noexcept { return loads_[dofOffset(dof)]; }
  const Eigen::Matrix<double,6,1>& loads() const noexcept { return loads_; }
 private:
  NodeId id_;
  Eigen::Vector3d coordinates_;
  std::array<bool,kDofsPerFrameNode> fixed_{};
  Eigen::Matrix<double,6,1> loads_ = Eigen::Matrix<double,6,1>::Zero();
};
}  // namespace fem
