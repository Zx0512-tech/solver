#pragma once

#include "fem/core/Element.hpp"

#include <Eigen/Core>

namespace fem {

class ElementLoad {
 public:
  explicit ElementLoad(ElementId element_id) : element_id_(element_id) {}
  virtual ~ElementLoad() = default;

  ElementId elementId() const noexcept { return element_id_; }

  // Returns the load vector in the element's GLOBAL DOF ordering.
  virtual Eigen::VectorXd equivalentNodalLoad(
      const Element& element,
      const NodeResolver& node) const = 0;

  // Static loads use the same vector at every time. Time-dependent wrappers
  // override this method while preserving the spatial load implementation.
  virtual Eigen::VectorXd equivalentNodalLoadAt(
      const Element& element,
      const NodeResolver& node,
      double time) const {
    (void)time;
    return equivalentNodalLoad(element, node);
  }

 private:
  ElementId element_id_;
};

}  // namespace fem
