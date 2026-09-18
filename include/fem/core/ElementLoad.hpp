#pragma once

#include "fem/core/Element.hpp"

#include <Eigen/Core>

#include <vector>

namespace fem {

struct ElementLoadSamplingLocation {
  double x{};
  bool has_jump{};
};

class ElementLoad {
 public:
  explicit ElementLoad(ElementId element_id) : element_id_(element_id) {}
  virtual ~ElementLoad() = default;

  ElementId elementId() const noexcept { return element_id_; }

  virtual Eigen::VectorXd equivalentNodalLoad(
      const Element& element,
      const NodeResolver& node) const = 0;

  virtual Eigen::VectorXd equivalentNodalLoadAt(
      const Element& element,
      const NodeResolver& node,
      double time) const {
    (void)time;
    return equivalentNodalLoad(element, node);
  }

  // Additional x-locations that improve response diagrams.
  // has_jump=true requests left/right limits at a discontinuity.
  virtual std::vector<ElementLoadSamplingLocation> responseSampleLocations(
      const Element& element,
      const NodeResolver& node) const {
    (void)element;
    (void)node;
    return {};
  }

 private:
  ElementId element_id_;
};

}  // namespace fem
