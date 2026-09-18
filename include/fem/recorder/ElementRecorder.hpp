#pragma once

#include "fem/core/ElementResponse.hpp"
#include "fem/core/Model.hpp"
#include "fem/solver/LinearStaticSolver.hpp"
#include "fem/solver/NewmarkBetaSolver.hpp"

#include <Eigen/Core>

#include <vector>

namespace fem {

struct ElementResponseHistory {
  ElementId element_id{};
  std::vector<double> time;
  Eigen::MatrixXd local_end_force;
  Eigen::MatrixXd global_end_force;
};

class ElementRecorder {
 public:
  ElementResponse record(
      const Model& model,
      ElementId element_id,
      const StaticResult& result) const;

  ElementResponseHistory record(
      const Model& model,
      ElementId element_id,
      const NewmarkResult& result) const;
};

}  // namespace fem
