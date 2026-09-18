#include "fem/core/DofManager.hpp"
#include <algorithm>
#include <stdexcept>
namespace fem {
DofManager::DofManager(const std::unordered_map<NodeId,Node>& nodes) {
  node_order_.reserve(nodes.size());
  for (const auto& [id,unused]:nodes) { (void)unused; node_order_.push_back(id); }
  std::sort(node_order_.begin(),node_order_.end());
  for (std::size_t i=0;i<node_order_.size();++i) node_base_.emplace(node_order_[i],i*kDofsPerFrameNode);
}
std::size_t DofManager::equation(NodeId node_id,Dof dof) const {
  const auto it=node_base_.find(node_id);
  if (it==node_base_.end()) throw std::out_of_range("Unknown node id in DOF manager");
  return it->second+dofOffset(dof);
}
}  // namespace fem
