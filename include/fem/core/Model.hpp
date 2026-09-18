#pragma once

#include "fem/core/DofManager.hpp"
#include "fem/core/Element.hpp"
#include "fem/core/Node.hpp"
#include "fem/core/Types.hpp"

#include <Eigen/Core>

#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace fem {

struct AssembledSystem {
  Eigen::MatrixXd stiffness;
  Eigen::MatrixXd mass;
  Eigen::VectorXd load;
  DofManager dofs;
};

class Model {
 public:
  Node& addNode(NodeId id, const Eigen::Vector3d& coordinates);

  template <typename ElementType, typename... Args>
  ElementType& addElement(Args&&... args) {
    auto element = std::make_unique<ElementType>(std::forward<Args>(args)...);
    const ElementId id = element->id();
    if (elements_.count(id) != 0U) {
      throw std::invalid_argument("Duplicate element id");
    }
    ElementType& ref = *element;
    elements_.emplace(id, std::move(element));
    return ref;
  }

  Node& node(NodeId id);
  const Node& node(NodeId id) const;

  AssembledSystem assemble() const;

  const std::unordered_map<NodeId, Node>& nodes() const noexcept { return nodes_; }

 private:
  std::unordered_map<NodeId, Node> nodes_;
  std::unordered_map<ElementId, std::unique_ptr<Element>> elements_;
};

}  // namespace fem
