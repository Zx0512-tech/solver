#pragma once

#include "fem/core/ElementLoad.hpp"

#include <Eigen/Core>

#include <functional>
#include <memory>

namespace fem {

// Applies a scalar time history to any existing spatial ElementLoad.
class TimeDependentElementLoad final : public ElementLoad {
 public:
  TimeDependentElementLoad(
      std::unique_ptr<ElementLoad> spatial_load,
      std::function<double(double)> scale);

  const ElementLoad& spatialLoad() const noexcept { return *spatial_load_; }
  double scaleAt(double time) const;

  Eigen::VectorXd equivalentNodalLoad(
      const Element& element,
      const NodeResolver& node) const override;

  Eigen::VectorXd equivalentNodalLoadAt(
      const Element& element,
      const NodeResolver& node,
      double time) const override;

  std::vector<ElementLoadSamplingLocation> responseSampleLocations(
      const Element& element,
      const NodeResolver& node) const override {
    return spatial_load_->responseSampleLocations(element, node);
  }

 private:
  static ElementId requireElementId(const std::unique_ptr<ElementLoad>& load);

  std::unique_ptr<ElementLoad> spatial_load_;
  std::function<double(double)> scale_;
};

}  // namespace fem
