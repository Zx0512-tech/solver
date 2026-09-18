#include "fem/loads/BeamPartialLinearLoad3D.hpp"

#include "fem/elements/Beam3D.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace fem {

Eigen::VectorXd BeamPartialLinearLoad3D::equivalentNodalLoad(
    const Element& element,
    const NodeResolver& node) const {
  const auto* beam = dynamic_cast<const Beam3D*>(&element);
  if (beam == nullptr) {
    throw std::invalid_argument(
        "BeamPartialLinearLoad3D can only target Beam3D elements");
  }

  const double l = beam->length(node);
  if (start_from_i_ < 0.0 || end_from_i_ > l || start_from_i_ >= end_from_i_) {
    throw std::invalid_argument(
        "BeamPartialLinearLoad3D interval must satisfy 0 <= start < end <= L");
  }

  const Beam3D::Matrix12d transform = beam->transformation(node);
  const Eigen::Matrix3d rotation = transform.block<3, 3>(0, 0);

  Eigen::Vector3d qa = load_per_length_start_;
  Eigen::Vector3d qb = load_per_length_end_;
  if (coordinate_system_ == BeamLoadCoordinateSystem::Global) {
    qa = rotation * qa;
    qb = rotation * qb;
  }

  Eigen::Matrix<double, 12, 1> local =
      Eigen::Matrix<double, 12, 1>::Zero();

  constexpr double root_three_fifths = 0.77459666924148337704;
  constexpr std::array<double, 3> xi = {
      -root_three_fifths, 0.0, root_three_fifths};
  constexpr std::array<double, 3> weight = {
      5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

  const double interval = end_from_i_ - start_from_i_;
  const double midpoint = 0.5 * (start_from_i_ + end_from_i_);
  const double jacobian = 0.5 * interval;

  for (std::size_t gp = 0; gp < xi.size(); ++gp) {
    const double x = midpoint + jacobian * xi[gp];
    const double eta = (x - start_from_i_) / interval;
    const Eigen::Vector3d q = (1.0 - eta) * qa + eta * qb;

    const double s = x / l;
    const double s2 = s * s;
    const double s3 = s2 * s;

    const double n_axial_i = 1.0 - s;
    const double n_axial_j = s;

    const double n1 = 1.0 - 3.0 * s2 + 2.0 * s3;
    const double n2 = l * (s - 2.0 * s2 + s3);
    const double n3 = 3.0 * s2 - 2.0 * s3;
    const double n4 = l * (-s2 + s3);

    const double dvol = jacobian * weight[gp];

    local[0] += n_axial_i * q.x() * dvol;
    local[6] += n_axial_j * q.x() * dvol;

    local[1] += n1 * q.y() * dvol;
    local[5] += n2 * q.y() * dvol;
    local[7] += n3 * q.y() * dvol;
    local[11] += n4 * q.y() * dvol;

    // Beam3D local-z interpolation uses theta_y = -duz/dx.
    local[2] += n1 * q.z() * dvol;
    local[4] += -n2 * q.z() * dvol;
    local[8] += n3 * q.z() * dvol;
    local[10] += -n4 * q.z() * dvol;
  }

  return transform.transpose() * local;
}

}  // namespace fem
