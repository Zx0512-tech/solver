#include "fem/loads/TimeDependentElementLoad.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fem {

ElementId TimeDependentElementLoad::requireElementId(
    const std::unique_ptr<ElementLoad>& load) {
  if (!load) {
    throw std::invalid_argument("TimeDependentElementLoad requires a spatial load");
  }
  return load->elementId();
}

TimeDependentElementLoad::TimeDependentElementLoad(
    std::unique_ptr<ElementLoad> spatial_load,
    std::function<double(double)> scale)
    : ElementLoad(requireElementId(spatial_load)),
      spatial_load_(std::move(spatial_load)),
      scale_(std::move(scale)) {
  if (!scale_) {
    throw std::invalid_argument("TimeDependentElementLoad requires a scale function");
  }
}

Eigen::VectorXd TimeDependentElementLoad::equivalentNodalLoad(
    const Element& element,
    const NodeResolver& node) const {
  return equivalentNodalLoadAt(element, node, 0.0);
}

Eigen::VectorXd TimeDependentElementLoad::equivalentNodalLoadAt(
    const Element& element,
    const NodeResolver& node,
    double time) const {
  const double factor = scale_(time);
  if (!std::isfinite(factor)) {
    throw std::runtime_error("Time-dependent element-load scale must be finite");
  }
  return factor * spatial_load_->equivalentNodalLoadAt(element, node, time);
}

}  // namespace fem
