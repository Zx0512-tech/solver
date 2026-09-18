#include "fem/loads/BeamPartialLinearLoad3D.hpp"

#include "fem/elements/Beam3D.hpp"

#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <stdexcept>

namespace fem {
namespace {

BeamLoadResultant3D integratePartialLinearLoad(
    double start,
    double end,
    double reference_x,
    double load_start,
    double load_end,
    const Eigen::Vector3d& qa,
    const Eigen::Vector3d& qb) {
  BeamLoadResultant3D result;
  if (end <= start) {
    return result;
  }

  constexpr double root_three_fifths = 0.77459666924148337704;
  constexpr std::array<double, 3> xi = {
      -root_three_fifths, 0.0, root_three_fifths};
  constexpr std::array<double, 3> weight = {
      5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

  const double midpoint = 0.5 * (start + end);
  const double jacobian = 0.5 * (end - start);
  const double load_span = load_end - load_start;
  const Eigen::Vector3d ex = Eigen::Vector3d::UnitX();

  for (std::size_t gp = 0; gp < xi.size(); ++gp) {
    const double s = midpoint + jacobian * xi[gp];
    const double eta = (s - load_start) / load_span;
    const Eigen::Vector3d q = (1.0 - eta) * qa + eta * qb;
    const double dv = jacobian * weight[gp];

    result.force += q * dv;
    result.moment += (s - reference_x) * ex.cross(q) * dv;
  }
  return result;
}

}  // namespace

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

    local[2] += n1 * q.z() * dvol;
    local[4] += -n2 * q.z() * dvol;
    local[8] += n3 * q.z() * dvol;
    local[10] += -n4 * q.z() * dvol;
  }

  return transform.transpose() * local;
}

BeamLoadResultant3D BeamPartialLinearLoad3D::localResultantTo(
    const Beam3D& beam,
    const NodeResolver& node,
    double x,
    BeamSectionSide side) const {
  (void)side;
  const double l = beam.length(node);
  if (start_from_i_ < 0.0 || end_from_i_ > l || start_from_i_ >= end_from_i_) {
    throw std::invalid_argument(
        "BeamPartialLinearLoad3D interval must satisfy 0 <= start < end <= L");
  }
  if (x < 0.0 || x > l) {
    throw std::invalid_argument("Beam section coordinate must lie in [0, L]");
  }

  const double integration_end = std::min(x, end_from_i_);
  if (integration_end <= start_from_i_) {
    return {};
  }

  Eigen::Vector3d qa = load_per_length_start_;
  Eigen::Vector3d qb = load_per_length_end_;
  if (coordinate_system_ == BeamLoadCoordinateSystem::Global) {
    const Eigen::Matrix3d rotation =
        beam.transformation(node).block<3, 3>(0, 0);
    qa = rotation * qa;
    qb = rotation * qb;
  }

  return integratePartialLinearLoad(
      start_from_i_,
      integration_end,
      x,
      start_from_i_,
      end_from_i_,
      qa,
      qb);
}

}  // namespace fem
