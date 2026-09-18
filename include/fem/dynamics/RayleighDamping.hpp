#pragma once

#include <Eigen/Core>

namespace fem {

class RayleighDamping {
 public:
  RayleighDamping(double alpha_mass = 0.0, double beta_stiffness = 0.0);

  static RayleighDamping fromModalTargets(double omega_1,
                                          double damping_ratio_1,
                                          double omega_2,
                                          double damping_ratio_2);

  double alphaMass() const noexcept { return alpha_mass_; }
  double betaStiffness() const noexcept { return beta_stiffness_; }

  double dampingRatio(double omega) const;

  Eigen::MatrixXd matrix(const Eigen::MatrixXd& mass,
                         const Eigen::MatrixXd& stiffness) const;

 private:
  double alpha_mass_;
  double beta_stiffness_;
};

}  // namespace fem
