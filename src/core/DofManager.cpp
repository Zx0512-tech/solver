#include "fem/core/DofManager.hpp"

#include <algorithm>
#include <stdexcept>

namespace fem {

DofManager::DofManager(const std::unordered_map<NodeId, Node>& nodes) {
  node_order_.reserve(nodes.size());
  for (const auto& [id, unused] : nodes) {
    (void)unused;
    node_order_.push_back(id);
  }
  std::sort(node_order_.begin(), node_order_.end());

  const std::size_t total_dofs = nodes.size() * kDofsPerFrameNode;
  free_equations_.reserve(total_dofs);
  constrained_equations_.reserve(total_dofs);

  for (std::size_t i = 0; i < node_order_.size(); ++i) {
    const NodeId node_id = node_order_[i];
    const std::size_t base = i * kDofsPerFrameNode;
    node_base_.emplace(node_id, base);

    const Node& node = nodes.at(node_id);
    for (std::size_t offset = 0; offset < kDofsPerFrameNode; ++offset) {
      const Eigen::Index eq =
          static_cast<Eigen::Index>(base + offset);
      if (node.isFixed(dofFromOffset(offset))) {
        constrained_equations_.push_back(eq);
      } else {
        free_equations_.push_back(eq);
      }
    }
  }
}

std::size_t DofManager::equation(NodeId node_id, Dof dof) const {
  const auto it = node_base_.find(node_id);
  if (it == node_base_.end()) {
    throw std::out_of_range("Unknown node id in DOF manager");
  }
  return it->second + dofOffset(dof);
}

std::vector<Eigen::Index> DofManager::equations(
    const std::vector<ElementDof>& dofs) const {
  std::vector<Eigen::Index> result;
  result.reserve(dofs.size());

  for (const auto& [node_id, dof] : dofs) {
    result.push_back(
        static_cast<Eigen::Index>(equation(node_id, dof)));
  }
  return result;
}

}  // namespace fem
