#include "fem/elements/Beam3D.hpp"

#include <Eigen/Geometry>

#include <cmath>
#include <stdexcept>

namespace fem {
namespace {
constexpr double kTol = 1.0e-12;
}

Beam3D::Beam3D(ElementId id,
               NodeId node_i,
               NodeId node_j,
               LinearElasticMaterial material,
               BeamSection section,
               const Eigen::Vector3d& local_y_hint)
    : Element(id),
      node_i_(node_i),
      node_j_(node_j),
      material_(material),
      section_(section),
      local_y_hint_(local_y_hint) {
  if (node_i_ == node_j_) {
    throw std::invalid_argument("Beam3D requires two distinct node ids");
  }
  if (local_y_hint_.norm() <= kTol) {
    throw std::invalid_argument("Beam3D local y-axis hint cannot be zero");
  }
}

std::vector<ElementDof> Beam3D::dofs() const {
  std::vector<ElementDof> result;
  result.reserve(12);
  for (const NodeId node_id : {node_i_, node_j_}) {
    for (std::size_t i = 0; i < kDofsPerFrameNode; ++i) {
      result.emplace_back(node_id, dofFromOffset(i));
    }
  }
  return result;
}

Eigen::MatrixXd Beam3D::stiffness(const NodeResolver& node) const {
  return globalStiffness(node);
}

Eigen::MatrixXd Beam3D::mass(const NodeResolver& node) const {
  return globalMass(node);
}

double Beam3D::length(const NodeResolver& node) const {
  const double l = (node(node_j_).coordinates() - node(node_i_).coordinates()).norm();
  if (l <= kTol) {
    throw std::runtime_error("Beam3D has zero or near-zero length");
  }
  return l;
}

Beam3D::Matrix12d Beam3D::localStiffness(double l) const {
  if (l <= kTol) {
    throw std::invalid_argument("Beam3D length must be positive");
  }

  const double e = material_.youngs_modulus;
  const double g = material_.shearModulus();
  const double a = section_.area;
  const double iy = section_.iy;
  const double iz = section_.iz;
  const double j = section_.torsion_constant;
  const double l2 = l * l;
  const double l3 = l2 * l;

  Matrix12d k = Matrix12d::Zero();

  const double ea_l = e * a / l;
  k(0, 0) = ea_l;
  k(0, 6) = -ea_l;
  k(6, 0) = -ea_l;
  k(6, 6) = ea_l;

  const double gj_l = g * j / l;
  k(3, 3) = gj_l;
  k(3, 9) = -gj_l;
  k(9, 3) = -gj_l;
  k(9, 9) = gj_l;

  const double z1 = 12.0 * e * iz / l3;
  const double z2 = 6.0 * e * iz / l2;
  const double z3 = 4.0 * e * iz / l;
  const double z4 = 2.0 * e * iz / l;
  const int z_dofs[4] = {1, 5, 7, 11};
  const double z_values[4][4] = {
      {z1, z2, -z1, z2},
      {z2, z3, -z2, z4},
      {-z1, -z2, z1, -z2},
      {z2, z4, -z2, z3},
  };
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) {
      k(z_dofs[r], z_dofs[c]) = z_values[r][c];
    }
  }

  const double y1 = 12.0 * e * iy / l3;
  const double y2 = 6.0 * e * iy / l2;
  const double y3 = 4.0 * e * iy / l;
  const double y4 = 2.0 * e * iy / l;
  const int y_dofs[4] = {2, 4, 8, 10};
  const double y_values[4][4] = {
      {y1, -y2, -y1, -y2},
      {-y2, y3, y2, y4},
      {-y1, y2, y1, y2},
      {-y2, y4, y2, y3},
  };
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) {
      k(y_dofs[r], y_dofs[c]) = y_values[r][c];
    }
  }

  return k;
}

Beam3D::Matrix12d Beam3D::localMass(double l) const {
  if (l <= kTol) {
    throw std::invalid_argument("Beam3D length must be positive");
  }

  const double rho = material_.density;
  const double area = section_.area;
  const double total_mass = rho * area * l;
  const double l2 = l * l;

  Matrix12d m = Matrix12d::Zero();

  const double axial = total_mass / 6.0;
  m(0, 0) = 2.0 * axial;
  m(0, 6) = axial;
  m(6, 0) = axial;
  m(6, 6) = 2.0 * axial;

  const double torsion = rho * (section_.iy + section_.iz) * l / 6.0;
  m(3, 3) = 2.0 * torsion;
  m(3, 9) = torsion;
  m(9, 3) = torsion;
  m(9, 9) = 2.0 * torsion;

  const double bending = total_mass / 420.0;

  const int z_dofs[4] = {1, 5, 7, 11};
  const double z_values[4][4] = {
      {156.0, 22.0 * l, 54.0, -13.0 * l},
      {22.0 * l, 4.0 * l2, 13.0 * l, -3.0 * l2},
      {54.0, 13.0 * l, 156.0, -22.0 * l},
      {-13.0 * l, -3.0 * l2, -22.0 * l, 4.0 * l2},
  };
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) {
      m(z_dofs[r], z_dofs[c]) = bending * z_values[r][c];
    }
  }

  const int y_dofs[4] = {2, 4, 8, 10};
  const double y_values[4][4] = {
      {156.0, -22.0 * l, 54.0, 13.0 * l},
      {-22.0 * l, 4.0 * l2, -13.0 * l, -3.0 * l2},
      {54.0, -13.0 * l, 156.0, 22.0 * l},
      {13.0 * l, -3.0 * l2, 22.0 * l, 4.0 * l2},
  };
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) {
      m(y_dofs[r], y_dofs[c]) = bending * y_values[r][c];
    }
  }

  return m;
}

Eigen::Matrix3d Beam3D::rotationToLocal(const NodeResolver& node) const {
  const Eigen::Vector3d x =
      (node(node_j_).coordinates() - node(node_i_).coordinates()).normalized();

  Eigen::Vector3d y_projected = local_y_hint_ - local_y_hint_.dot(x) * x;
  if (y_projected.norm() <= kTol) {
    const Eigen::Vector3d axes[3] = {
        Eigen::Vector3d::UnitX(), Eigen::Vector3d::UnitY(), Eigen::Vector3d::UnitZ()};
    int best = 0;
    double min_alignment = std::abs(x.dot(axes[0]));
    for (int i = 1; i < 3; ++i) {
      const double alignment = std::abs(x.dot(axes[i]));
      if (alignment < min_alignment) {
        min_alignment = alignment;
        best = i;
      }
    }
    y_projected = axes[best] - axes[best].dot(x) * x;
  }

  const Eigen::Vector3d y = y_projected.normalized();
  const Eigen::Vector3d z = x.cross(y).normalized();
  const Eigen::Vector3d corrected_y = z.cross(x).normalized();

  Eigen::Matrix3d r;
  r.row(0) = x.transpose();
  r.row(1) = corrected_y.transpose();
  r.row(2) = z.transpose();
  return r;
}

Beam3D::Matrix12d Beam3D::transformation(const NodeResolver& node) const {
  const Eigen::Matrix3d r = rotationToLocal(node);
  Matrix12d t = Matrix12d::Zero();
  for (int block = 0; block < 4; ++block) {
    t.block<3, 3>(3 * block, 3 * block) = r;
  }
  return t;
}

Beam3D::Matrix12d Beam3D::globalStiffness(const NodeResolver& node) const {
  const Matrix12d t = transformation(node);
  return t.transpose() * localStiffness(length(node)) * t;
}

Beam3D::Matrix12d Beam3D::globalMass(const NodeResolver& node) const {
  const Matrix12d t = transformation(node);
  return t.transpose() * localMass(length(node)) * t;
}

}  // namespace fem
