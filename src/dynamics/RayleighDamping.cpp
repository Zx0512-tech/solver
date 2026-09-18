#include "fem/dynamics/RayleighDamping.hpp"

#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fem {

RayleighDamping::RayleighDamping(double alpha_mass, double beta_stiffness)
    : alpha_mass_(alpha_mass), beta_stiffness_(beta_stiffness) {
  if (!std::isfinite(alpha_mass_) || !std::isfinite(beta_stiffness_)) {
    throw std::invalid_argument("Rayleigh damping coefficients must be finite");
  }
}

RayleighDamping RayleighDamping::fromModalTargets(double omega_1,
                                                  double damping_ratio_1,
                                                  double omega_2,
                                                  double damping_ratio_2) {
  if (omega_1 <= 0.0 || omega_2 <= 0.0) {
    throw std::invalid_argument("Target circular frequencies must be positive");
  }
  if (damping_ratio_1 < 0.0 || damping_ratio_2 < 0.0) {
    throw std::invalid_argument("Target damping ratios cannot be negative");
  }
  if (std::abs(omega_1 - omega_2) <= 1.0e-12 * std::max(omega_1, omega_2)) {
    throw std::invalid_argument("Rayleigh target frequencies must be distinct");
  }

  Eigen::Matrix2d a;
  a << 1.0 / (2.0 * omega_1), omega_1 / 2.0,
       1.0 / (2.0 * omega_2), omega_2 / 2.0;
  const Eigen::Vector2d zeta(damping_ratio_1, damping_ratio_2);
  const Eigen::Vector2d coefficients = a.fullPivLu().solve(zeta);

  return {coefficients[0], coefficients[1]};
}

double RayleighDamping::dampingRatio(double omega) const {
  if (omega <= 0.0) {
    throw std::invalid_argument("Circular frequency must be positive");
  }
  return alpha_mass_ / (2.0 * omega) + beta_stiffness_ * omega / 2.0;
}

Eigen::MatrixXd RayleighDamping::matrix(const Eigen::MatrixXd& mass,
                                        const Eigen::MatrixXd& stiffness) const {
  if (mass.rows() != mass.cols() || stiffness.rows() != stiffness.cols() ||
      mass.rows() != stiffness.rows()) {
    throw std::invalid_argument("Mass and stiffness matrices must be square and the same size");
  }
  return alpha_mass_ * mass + beta_stiffness_ * stiffness;
}

}  // namespace fem
