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
  const double qx = load_per_length_.x();
  const double qy = load_per_length_.y();
  const double qz = load_per_length_.z();

  Eigen::Matrix<double, 12, 1> local =
      Eigen::Matrix<double, 12, 1>::Zero();

  // Axial consistent load.
  local[0] = qx * l / 2.0;
  local[6] = qx * l / 2.0;

  // Uniform local-y load -> bending about local z.
  local[1] = qy * l / 2.0;
  local[5] = qy * l2 / 12.0;
  local[7] = qy * l / 2.0;
  local[11] = -qy * l2 / 12.0;

  // Uniform local-z load -> bending about local y.
  // Beam3D uses the sign convention where positive uz bending gives negative ry.
  local[2] = qz * l / 2.0;
  local[4] = -qz * l2 / 12.0;
  local[8] = qz * l / 2.0;
  local[10] = qz * l2 / 12.0;

  const Beam3D::Matrix12d transform = beam->transformation(node);
  return transform.transpose() * local;
}

}  // namespace fem
