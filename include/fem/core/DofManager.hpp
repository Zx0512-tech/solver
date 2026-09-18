#pragma once
#include "fem/core/Node.hpp"
#include "fem/core/Types.hpp"
#include <cstddef>
#include <unordered_map>
#include <vector>
namespace fem {
class DofManager {
 public:
  DofManager() = default;
  explicit DofManager(const std::unordered_map<NodeId,Node>& nodes);
  std::size_t equation(NodeId node_id,Dof dof) const;
  std::size_t size() const noexcept { return node_base_.size()*kDofsPerFrameNode; }
  const std::vector<NodeId>& nodeOrder() const noexcept { return node_order_; }
 private:
  std::unordered_map<NodeId,std::size_t> node_base_;
  std::vector<NodeId> node_order_;
};
}  // namespace fem
