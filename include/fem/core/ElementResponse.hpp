#pragma once

#include "fem/core/Types.hpp"

#include <Eigen/Core>

namespace fem {

struct ElementResponse {
  ElementId element_id{};
  Eigen::VectorXd local_deformation;
  Eigen::VectorXd local_end_force;
  Eigen::VectorXd global_end_force;
};

}  // namespace fem
