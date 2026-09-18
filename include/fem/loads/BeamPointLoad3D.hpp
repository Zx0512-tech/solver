#pragma once

#include "fem/loads/BeamElementLoad3D.hpp"
#include "fem/loads/BeamLoadCoordinateSystem.hpp"

#include <Eigen/Core>

namespace fem {

// Concentrated force and/or moment at an arbitrary distance from Beam3D node i.
class BeamPointLoad3D final : public BeamElementLoad3D {
 public:
  BeamPointLoad3D(
      ElementId element_id,
      double distance_from_i,
      const Eigen::Vector3d& force,
      const Eigen::Vector3d& moment = Eigen::Vector3d::Zero(),
      BeamLoadCoordinateSystem coordinate_system = BeamLoadCoordinateSystem::Local)
      : BeamElementLoad3D(element_id),
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

  BeamLoadResultant3D localResultantTo(
      const Beam3D& beam,
      const NodeResolver& node,
      double x,
      BeamSectionSide side) const override;

 private:
  double distance_from_i_;
  Eigen::Vector3d force_;
  Eigen::Vector3d moment_;
  BeamLoadCoordinateSystem coordinate_system_;
};

}  // namespace fem
