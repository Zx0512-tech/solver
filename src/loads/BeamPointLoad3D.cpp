#include "fem/loads/BeamPointLoad3D.hpp"

#include "fem/elements/Beam3D.hpp"

#include <stdexcept>

namespace fem {

Eigen::VectorXd BeamPointLoad3D::equivalentNodalLoad(
    const Element& element,
    const NodeResolver& node) const {
  const auto* beam = dynamic_cast<const Beam3D*>(&element);
  if (beam == nullptr) {
    throw std::invalid_argument("BeamPointLoad3D can only target Beam3D elements");
  }

  const double l = beam->length(node);
  if (distance_from_i_ < 0.0 || distance_from_i_ > l) {
    throw std::invalid_argument(
        "BeamPointLoad3D distance must lie between node i and node j");
  }

  const Beam3D::Matrix12d transform = beam->transformation(node);
  const Eigen::Matrix3d rotation = transform.block<3, 3>(0, 0);

  Eigen::Vector3d force = force_;
  Eigen::Vector3d moment = moment_;
  if (coordinate_system_ == BeamLoadCoordinateSystem::Global) {
    force = rotation * force_;
    moment = rotation * moment_;
  }

  const double s = distance_from_i_ / l;
  const double s2 = s * s;
  const double s3 = s2 * s;

  const double n1 = 1.0 - 3.0 * s2 + 2.0 * s3;
  const double n2 = l * (s - 2.0 * s2 + s3);
  const double n3 = 3.0 * s2 - 2.0 * s3;
  const double n4 = l * (-s2 + s3);

  const double dn1 = 6.0 * (s2 - s) / l;
  const double dn2 = 1.0 - 4.0 * s + 3.0 * s2;
  const double dn3 = 6.0 * (s - s2) / l;
  const double dn4 = -2.0 * s + 3.0 * s2;

  Eigen::Matrix<double, 12, 1> local =
      Eigen::Matrix<double, 12, 1>::Zero();

  local[0] += (1.0 - s) * force.x();
  local[6] += s * force.x();

  local[1] += n1 * force.y();
  local[5] += n2 * force.y();
  local[7] += n3 * force.y();
  local[11] += n4 * force.y();

  local[2] += n1 * force.z();
  local[4] += -n2 * force.z();
  local[8] += n3 * force.z();
  local[10] += -n4 * force.z();

  local[3] += (1.0 - s) * moment.x();
  local[9] += s * moment.x();

  local[1] += dn1 * moment.z();
  local[5] += dn2 * moment.z();
  local[7] += dn3 * moment.z();
  local[11] += dn4 * moment.z();

  local[2] += -dn1 * moment.y();
  local[4] += dn2 * moment.y();
  local[8] += -dn3 * moment.y();
  local[10] += dn4 * moment.y();

  return transform.transpose() * local;
}

BeamLoadResultant3D BeamPointLoad3D::localResultantTo(
    const Beam3D& beam,
    const NodeResolver& node,
    double x,
    BeamSectionSide side) const {
  const double l = beam.length(node);
  if (distance_from_i_ < 0.0 || distance_from_i_ > l) {
    throw std::invalid_argument(
        "BeamPointLoad3D distance must lie between node i and node j");
  }
  if (x < 0.0 || x > l) {
    throw std::invalid_argument("Beam section coordinate must lie in [0, L]");
  }

  const bool before_section = distance_from_i_ < x;
  const bool exactly_at_section = distance_from_i_ == x;
  if (!before_section &&
      !(exactly_at_section && side == BeamSectionSide::Right)) {
    return {};
  }

  Eigen::Vector3d force = force_;
  Eigen::Vector3d moment = moment_;
  if (coordinate_system_ == BeamLoadCoordinateSystem::Global) {
    const Eigen::Matrix3d rotation =
        beam.transformation(node).block<3, 3>(0, 0);
    force = rotation * force;
    moment = rotation * moment;
  }

  BeamLoadResultant3D result;
  result.force = force;
  result.moment =
      moment +
      (distance_from_i_ - x) * Eigen::Vector3d::UnitX().cross(force);
  return result;
}

}  // namespace fem
