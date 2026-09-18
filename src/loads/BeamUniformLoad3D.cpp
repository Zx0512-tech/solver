#include "fem/loads/BeamUniformLoad3D.hpp"

#include "fem/elements/Beam3D.hpp"

#include <stdexcept>

namespace fem {

Eigen::VectorXd BeamUniformLoad3D::equivalentNodalLoad(
    const Element& element,
    const NodeResolver& node) const {
  const auto* beam = dynamic_cast<const Beam3D*>(&element);
  if (beam == nullptr) {
    throw std::invalid_argument("BeamUniformLoad3D can only target Beam3D elements");
  }

  const double l = beam->length(node);
  const double l2 = l * l;
  const Beam3D::Matrix12d transform = beam->transformation(node);

  Eigen::Vector3d local_q = load_per_length_;
  if (coordinate_system_ == BeamLoadCoordinateSystem::Global) {
    local_q = transform.block<3, 3>(0, 0) * load_per_length_;
  }

  const double qx = local_q.x();
  const double qy = local_q.y();
  const double qz = local_q.z();

  Eigen::Matrix<double, 12, 1> local =
      Eigen::Matrix<double, 12, 1>::Zero();

  local[0] = qx * l / 2.0;
  local[6] = qx * l / 2.0;

  local[1] = qy * l / 2.0;
  local[5] = qy * l2 / 12.0;
  local[7] = qy * l / 2.0;
  local[11] = -qy * l2 / 12.0;

  // Positive local uz corresponds to negative local ry slope.
  local[2] = qz * l / 2.0;
  local[4] = -qz * l2 / 12.0;
  local[8] = qz * l / 2.0;
  local[10] = qz * l2 / 12.0;

  return transform.transpose() * local;
}

}  // namespace fem
