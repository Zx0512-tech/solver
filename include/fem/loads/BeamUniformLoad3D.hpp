#pragma once

#include "fem/core/ElementLoad.hpp"

#include <Eigen/Core>

namespace fem {

// Uniform line load expressed in the Beam3D LOCAL axes.
// q = [qx, qy, qz] in force/length.
class BeamUniformLoad3D final : public ElementLoad {
 public:
  BeamUniformLoad3D(ElementId element_id, const Eigen::Vector3d& load_per_length)
      : ElementLoad(element_id), load_per_length_(load_per_length) {}

  const Eigen::Vector3d& loadPerLength() const noexcept { return load_per_length_; }

  Eigen::VectorXd equivalentNodalLoad(
      const Element& element,
      const NodeResolver& node) const override;

 private:
  Eigen::Vector3d load_per_length_;
};

}  // namespace fem
