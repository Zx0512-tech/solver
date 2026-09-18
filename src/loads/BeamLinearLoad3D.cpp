#include "fem/loads/BeamLinearLoad3D.hpp"

#include "fem/elements/Beam3D.hpp"

#include <stdexcept>

namespace fem {

Eigen::VectorXd BeamLinearLoad3D::equivalentNodalLoad(
    const Element& element,
    const NodeResolver& node) const {
  const auto* beam = dynamic_cast<const Beam3D*>(&element);
  if (beam == nullptr) {
    throw std::invalid_argument("BeamLinearLoad3D can only target Beam3D elements");
  }

  const double l = beam->length(node);
  const double l2 = l * l;
  const Beam3D::Matrix12d transform = beam->transformation(node);
  const Eigen::Matrix3d rotation = transform.block<3, 3>(0, 0);

  Eigen::Vector3d qi = load_per_length_i_;
  Eigen::Vector3d qj = load_per_length_j_;
  if (coordinate_system_ == BeamLoadCoordinateSystem::Global) {
    qi = rotation * load_per_length_i_;
    qj = rotation * load_per_length_j_;
  }

  Eigen::Matrix<double, 12, 1> local =
      Eigen::Matrix<double, 12, 1>::Zero();

  // Linear axial interpolation.
  local[0] = l * (2.0 * qi.x() + qj.x()) / 6.0;
  local[6] = l * (qi.x() + 2.0 * qj.x()) / 6.0;

  // Hermite-consistent local-y loading.
  local[1] = l * (7.0 * qi.y() + 3.0 * qj.y()) / 20.0;
  local[5] = l2 * (qi.y() / 20.0 + qj.y() / 30.0);
  local[7] = l * (3.0 * qi.y() + 7.0 * qj.y()) / 20.0;
  local[11] = -l2 * (qi.y() / 30.0 + qj.y() / 20.0);

  // Same interpolation for local-z, with Beam3D's ry sign convention.
  local[2] = l * (7.0 * qi.z() + 3.0 * qj.z()) / 20.0;
  local[4] = -l2 * (qi.z() / 20.0 + qj.z() / 30.0);
  local[8] = l * (3.0 * qi.z() + 7.0 * qj.z()) / 20.0;
  local[10] = l2 * (qi.z() / 30.0 + qj.z() / 20.0);

  return transform.transpose() * local;
}

}  // namespace fem
