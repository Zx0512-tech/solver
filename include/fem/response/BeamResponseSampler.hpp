#pragma once

#include "fem/core/Types.hpp"
#include "fem/response/BeamSectionForces.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace fem {

class Model;

struct BeamForceSample {
  BeamSectionForces forces;
  BeamSectionSide side{BeamSectionSide::Right};
};

using BeamForceSeries = std::vector<BeamForceSample>;

class BeamResponseSampler {
 public:
  BeamForceSeries sample(
      const Model& model,
      ElementId element_id,
      const Eigen::VectorXd& global_displacement,
      std::size_t point_count,
      double time = 0.0) const;
};

}  // namespace fem
