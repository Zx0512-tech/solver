#pragma once

#include "fem/core/ElementLoad.hpp"
#include "fem/loads/BeamLoadCoordinateSystem.hpp"

#include <Eigen/Core>

namespace fem {

// Uniform line load q = [qx, qy, qz] in force/length.
// Local is the default; Global is useful for gravity-like loads on inclined beams.
class BeamUniformLoad3D final : public ElementLoad {
 public:
  BeamUniformLoad3D(
      ElementId element_id,
      const Eigen::Vector3d& load_per_length,
      BeamLoadCoordinateSystem coordinate_system = BeamLoadCoordinateSystem::Local)
      : ElementLoad(element_id),
        load_per_length_(load_per_length),
        coordinate_system_(coordinate_system) {}

  const Eigen::Vector3d& loadPerLength() const noexcept { return load_per_length_; }
  BeamLoadCoordinateSystem coordinateSystem() const noexcept { return coordinate_system_; }

  Eigen::VectorXd equivalentNodalLoad(
      const Element& element,
      const NodeResolver& node) const override;

 private:
  Eigen::Vector3d load_per_length_;
  BeamLoadCoordinateSystem coordinate_system_;
};

}  // namespace fem
