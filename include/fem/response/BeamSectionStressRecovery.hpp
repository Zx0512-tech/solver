#pragma once

#include "fem/core/Section.hpp"
#include "fem/response/BeamSectionForces.hpp"
#include "fem/response/BeamSectionStress.hpp"

namespace fem {

class BeamSectionStressRecovery {
 public:
  double normalStressAt(
      const BeamSection& section,
      const BeamSectionForces& forces,
      double y,
      double z) const;

  BeamSectionStress stressAt(
      const BeamSection& section,
      const BeamSectionForces& forces,
      double y,
      double z) const;

  BeamNormalStressExtrema normalStressExtrema(
      const BeamSection& section,
      const BeamSectionForces& forces) const;
};

}  // namespace fem
