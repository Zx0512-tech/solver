#pragma once

#include "fem/core/ElementLoad.hpp"
#include "fem/loads/BeamLoadCoordinateSystem.hpp"

#include <Eigen/Core>

namespace fem {

// Linearly varying line load over the complete Beam3D.
// q_i and q_j are the force/length vectors at node i and node j.
// Equal endpoints reduce exactly to a uniform load; unequal endpoints represent
// triangular or trapezoidal loading.
class BeamLinearLoad3D final : public ElementLoad {
 public:
  BeamLinearLoad3D(
      ElementId element_id,
      const Eigen::Vector3d& load_per_length_i,
      const Eigen::Vector3d& load_per_length_j,
      BeamLoadCoordinateSystem coordinate_system = BeamLoadCoordinateSystem::Local)
      : ElementLoad(element_id),
        load_per_length_i_(load_per_length_i),
        load_per_length_j_(load_per_length_j),
        coordinate_system_(coordinate_system) {}

  const Eigen::Vector3d& loadPerLengthI() const noexcept { return load_per_length_i_; }
  const Eigen::Vector3d& loadPerLengthJ() const noexcept { return load_per_length_j_; }
  BeamLoadCoordinateSystem coordinateSystem() const noexcept { return coordinate_system_; }

  Eigen::VectorXd equivalentNodalLoad(
      const Element& element,
      const NodeResolver& node) const override;

 private:
  Eigen::Vector3d load_per_length_i_;
  Eigen::Vector3d load_per_length_j_;
  BeamLoadCoordinateSystem coordinate_system_;
};

using BeamTrapezoidalLoad3D = BeamLinearLoad3D;

}  // namespace fem
