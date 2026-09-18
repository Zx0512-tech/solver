#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <vector>

namespace fem {

class SparseDofReducer {
 public:
  static Eigen::SparseMatrix<double> matrix(
      const Eigen::SparseMatrix<double>& full,
      const std::vector<Eigen::Index>& indices);

  static Eigen::VectorXd vector(
      const Eigen::VectorXd& full,
      const std::vector<Eigen::Index>& indices);
};

}  // namespace fem
