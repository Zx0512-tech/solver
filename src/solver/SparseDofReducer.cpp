#include "fem/solver/SparseDofReducer.hpp"

#include <stdexcept>
#include <vector>

namespace fem {
namespace {

void validateIndices(
    Eigen::Index full_size,
    const std::vector<Eigen::Index>& indices) {
  std::vector<unsigned char> seen(
      static_cast<std::size_t>(full_size), 0U);

  for (const Eigen::Index index : indices) {
    if (index < 0 || index >= full_size) {
      throw std::out_of_range(
          "Reduced-system DOF index is out of range");
    }
    auto& flag = seen[static_cast<std::size_t>(index)];
    if (flag != 0U) {
      throw std::invalid_argument(
          "Reduced-system DOF indices must be unique");
    }
    flag = 1U;
  }
}

}  // namespace

Eigen::SparseMatrix<double> SparseDofReducer::matrix(
    const Eigen::SparseMatrix<double>& full,
    const std::vector<Eigen::Index>& indices) {
  if (full.rows() != full.cols()) {
    throw std::invalid_argument(
        "SparseDofReducer requires a square matrix");
  }

  validateIndices(full.rows(), indices);

  const Eigen::Index reduced_size =
      static_cast<Eigen::Index>(indices.size());

  std::vector<Eigen::Index> global_to_local(
      static_cast<std::size_t>(full.rows()),
      Eigen::Index{-1});
  for (Eigen::Index local = 0; local < reduced_size; ++local) {
    global_to_local[
        static_cast<std::size_t>(
            indices[static_cast<std::size_t>(local)])] = local;
  }

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(full.nonZeros()));

  for (int outer = 0; outer < full.outerSize(); ++outer) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(full, outer);
         it;
         ++it) {
      const Eigen::Index local_row =
          global_to_local[static_cast<std::size_t>(it.row())];
      const Eigen::Index local_col =
          global_to_local[static_cast<std::size_t>(it.col())];

      if (local_row >= 0 && local_col >= 0) {
        triplets.emplace_back(
            local_row, local_col, it.value());
      }
    }
  }

  Eigen::SparseMatrix<double> reduced(
      reduced_size, reduced_size);
  reduced.setFromTriplets(
      triplets.begin(), triplets.end(),
      [](double a, double b) { return a + b; });
  reduced.prune(0.0);
  reduced.makeCompressed();
  return reduced;
}

Eigen::VectorXd SparseDofReducer::vector(
    const Eigen::VectorXd& full,
    const std::vector<Eigen::Index>& indices) {
  validateIndices(full.size(), indices);

  Eigen::VectorXd reduced(
      static_cast<Eigen::Index>(indices.size()));
  for (std::size_t i = 0; i < indices.size(); ++i) {
    reduced[static_cast<Eigen::Index>(i)] =
        full[indices[i]];
  }
  return reduced;
}

}  // namespace fem
