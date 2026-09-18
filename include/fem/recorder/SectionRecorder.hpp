#pragma once

#include "fem/core/Model.hpp"
#include "fem/recorder/ResultHistory.hpp"
#include "fem/response/BeamSectionForces.hpp"
#include "fem/solver/LinearStaticSolver.hpp"
#include "fem/solver/NewmarkBetaSolver.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace fem {

enum class BeamSectionForceComponent : std::size_t {
  N = 0,
  Vy = 1,
  Vz = 2,
  T = 3,
  My = 4,
  Mz = 5,
};

struct SectionResponse {
  ElementId element_id{};
  double x{};
  BeamSectionSide side{BeamSectionSide::Right};
  BeamSectionForces forces;
};

struct SectionResponseHistory {
  ElementId element_id{};
  double x{};
  BeamSectionSide side{BeamSectionSide::Right};
  std::vector<double> time;
  Eigen::MatrixXd forces;

  double forceAt(
      std::size_t step,
      BeamSectionForceComponent component) const;

  ScalarHistory series(
      BeamSectionForceComponent component) const;
};

class SectionRecorder {
 public:
  SectionResponse record(
      const Model& model,
      ElementId element_id,
      double x,
      const StaticResult& result,
      BeamSectionSide side = BeamSectionSide::Right) const;

  SectionResponseHistory record(
      const Model& model,
      ElementId element_id,
      double x,
      const NewmarkResult& result,
      BeamSectionSide side = BeamSectionSide::Right) const;

  ScalarHistory recordNormalStress(
      const Model& model,
      ElementId element_id,
      double x,
      double y,
      double z,
      const NewmarkResult& result,
      BeamSectionSide side = BeamSectionSide::Right) const;
};

}  // namespace fem
