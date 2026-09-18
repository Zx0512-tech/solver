#pragma once

#include "fem/loads/BeamElementLoad3D.hpp"
#include "fem/loads/BeamLoadCoordinateSystem.hpp"

#include <Eigen/Core>

namespace fem {

// Linearly varying distributed load applied only on [start_from_i, end_from_i].
class BeamPartialLinearLoad3D final : public BeamElementLoad3D {
 public:
  BeamPartialLinearLoad3D(
      ElementId element_id,
      double start_from_i,
      double end_from_i,
      const Eigen::Vector3d& load_per_length_start,
      const Eigen::Vector3d& load_per_length_end,
      BeamLoadCoordinateSystem coordinate_system = BeamLoadCoordinateSystem::Local)
      : BeamElementLoad3D(element_id),
        start_from_i_(start_from_i),
        end_from_i_(end_from_i),
        load_per_length_start_(load_per_length_start),
        load_per_length_end_(load_per_length_end),
        coordinate_system_(coordinate_system) {}

  double startFromI() const noexcept { return start_from_i_; }
  double endFromI() const noexcept { return end_from_i_; }

  Eigen::VectorXd equivalentNodalLoad(
      const Element& element,
      const NodeResolver& node) const override;

  BeamLoadResultant3D localResultantTo(
      const Beam3D& beam,
      const NodeResolver& node,
      double x,
      BeamSectionSide side) const override;

 private:
  double start_from_i_;
  double end_from_i_;
  Eigen::Vector3d load_per_length_start_;
  Eigen::Vector3d load_per_length_end_;
  BeamLoadCoordinateSystem coordinate_system_;
};

}  // namespace fem
