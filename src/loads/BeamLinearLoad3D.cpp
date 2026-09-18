#include "fem/loads/BeamLinearLoad3D.hpp"

#include "fem/elements/Beam3D.hpp"

#include <array>
#include <stdexcept>

namespace fem {
namespace {

BeamLoadResultant3D integrateLinearLoad(
    double start,
    double end,
    double reference_x,
    const Eigen::Vector3d& qi,
    const Eigen::Vector3d& qj,
    double full_length) {
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
  const Eigen::Vector3d ex = Eigen::Vector3d::UnitX();

  for (std::size_t gp = 0; gp < xi.size(); ++gp) {
    const double s = midpoint + jacobian * xi[gp];
    const double eta = s / full_length;
    const Eigen::Vector3d q = (1.0 - eta) * qi + eta * qj;
    const double dv = jacobian * weight[gp];

    result.force += q * dv;
    result.moment += (s - reference_x) * ex.cross(q) * dv;
  }
  return result;
}

}  // namespace

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

  local[0] = l * (2.0 * qi.x() + qj.x()) / 6.0;
  local[6] = l * (qi.x() + 2.0 * qj.x()) / 6.0;

  local[1] = l * (7.0 * qi.y() + 3.0 * qj.y()) / 20.0;
  local[5] = l2 * (qi.y() / 20.0 + qj.y() / 30.0);
  local[7] = l * (3.0 * qi.y() + 7.0 * qj.y()) / 20.0;
  local[11] = -l2 * (qi.y() / 30.0 + qj.y() / 20.0);

  local[2] = l * (7.0 * qi.z() + 3.0 * qj.z()) / 20.0;
  local[4] = -l2 * (qi.z() / 20.0 + qj.z() / 30.0);
  local[8] = l * (3.0 * qi.z() + 7.0 * qj.z()) / 20.0;
  local[10] = l2 * (qi.z() / 30.0 + qj.z() / 20.0);

  return transform.transpose() * local;
}

BeamLoadResultant3D BeamLinearLoad3D::localResultantTo(
    const Beam3D& beam,
    const NodeResolver& node,
    double x,
    BeamSectionSide side) const {
  (void)side;
  const double l = beam.length(node);
  if (x < 0.0 || x > l) {
    throw std::invalid_argument("Beam section coordinate must lie in [0, L]");
  }

  Eigen::Vector3d qi = load_per_length_i_;
  Eigen::Vector3d qj = load_per_length_j_;
  if (coordinate_system_ == BeamLoadCoordinateSystem::Global) {
    const Eigen::Matrix3d rotation =
        beam.transformation(node).block<3, 3>(0, 0);
    qi = rotation * qi;
    qj = rotation * qj;
  }

  return integrateLinearLoad(0.0, x, x, qi, qj, l);
}

}  // namespace fem
