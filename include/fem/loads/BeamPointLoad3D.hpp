#pragma once

#include "fem/core/ElementLoad.hpp"
#include "fem/loads/BeamLoadCoordinateSystem.hpp"

#include <Eigen/Core>

namespace fem {

// Concentrated force and/or moment at an arbitrary distance from Beam3D node i.
// force and moment may be supplied in local or global axes.
class BeamPointLoad3D final : public ElementLoad {
 public:
  BeamPointLoad3D(
      ElementId element_id,
      double distance_from_i,
      const Eigen::Vector3d& force,
      const Eigen::Vector3d& moment = Eigen::Vector3d::Zero(),
      BeamLoadCoordinateSystem coordinate_system = BeamLoadCoordinateSystem::Local)
      : ElementLoad(element_id),
        distance_from_i_(distance_from_i),
        force_(force),
        moment_(moment),
        coordinate_system_(coordinate_system) {}

  double distanceFromI() const noexcept { return distance_from_i_; }
  const Eigen::Vector3d& force() const noexcept { return force_; }
  const Eigen::Vector3d& moment() const noexcept { return moment_; }
  BeamLoadCoordinateSystem coordinateSystem() const noexcept { return coordinate_system_; }

  Eigen::VectorXd equivalentNodalLoad(
      const Element& element,
      const NodeResolver& node) const override;

 private:
  double distance_from_i_;
  Eigen::Vector3d force_;
  Eigen::Vector3d moment_;
  BeamLoadCoordinateSystem coordinate_system_;
};

}  // namespace fem
