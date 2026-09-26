#pragma once

#include "fem/core/Element.hpp"
#include "fem/core/Node.hpp"
#include "fem/core/Types.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <unordered_map>
#include <vector>

namespace fem {

class DofManager {
 public:
  DofManager() = default;
  explicit DofManager(const std::unordered_map<NodeId, Node>& nodes);

  std::size_t equation(NodeId node_id, Dof dof) const;

  // Maps element DOFs to global equations without sorting or deduplication.
  // Output order always matches the input order because element matrix
  // rows/columns are defined by Element::dofs().
  std::vector<Eigen::Index> equations(
      const std::vector<ElementDof>& dofs) const;

  std::size_t size() const noexcept {
    return node_base_.size() * kDofsPerFrameNode;
  }

  const std::vector<NodeId>& nodeOrder() const noexcept {
    return node_order_;
  }

  // These sets are a snapshot of the nodal constraint state at construction
  // time. Model::assemble() constructs a fresh DofManager for each analysis.
  const std::vector<Eigen::Index>& freeEquations() const noexcept {
    return free_equations_;
  }

  const std::vector<Eigen::Index>& constrainedEquations() const noexcept {
    return constrained_equations_;
  }

 private:
  std::unordered_map<NodeId, std::size_t> node_base_;
  std::vector<NodeId> node_order_;
  std::vector<Eigen::Index> free_equations_;
  std::vector<Eigen::Index> constrained_equations_;
};

}  // namespace fem
