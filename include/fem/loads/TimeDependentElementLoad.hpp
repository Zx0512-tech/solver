#pragma once

#include "fem/core/ElementLoad.hpp"

#include <Eigen/Core>

#include <functional>
#include <memory>

namespace fem {

// Applies a scalar time history to any existing spatial ElementLoad.
// Example: P(t) = P0 * scale(t), q(x,t) = q0(x) * scale(t).
class TimeDependentElementLoad final : public ElementLoad {
 public:
  TimeDependentElementLoad(
      std::unique_ptr<ElementLoad> spatial_load,
      std::function<double(double)> scale);

  const ElementLoad& spatialLoad() const noexcept { return *spatial_load_; }

  Eigen::VectorXd equivalentNodalLoad(
      const Element& element,
      const NodeResolver& node) const override;

  Eigen::VectorXd equivalentNodalLoadAt(
      const Element& element,
      const NodeResolver& node,
      double time) const override;

 private:
  static ElementId requireElementId(const std::unique_ptr<ElementLoad>& load);

  std::unique_ptr<ElementLoad> spatial_load_;
  std::function<double(double)> scale_;
};

}  // namespace fem
