#pragma once

#include "fem/core/ElementLoad.hpp"
#include "fem/loads/BeamLoadResultant3D.hpp"
#include "fem/response/BeamSectionForces.hpp"

namespace fem {

class Beam3D;

// Common interface for spatial Beam3D loads that can participate in section
// equilibrium recovery in addition to providing equivalent nodal loads.
class BeamElementLoad3D : public ElementLoad {
 public:
  using ElementLoad::ElementLoad;
  ~BeamElementLoad3D() override = default;

  virtual BeamLoadResultant3D localResultantTo(
      const Beam3D& beam,
      const NodeResolver& node,
      double x,
      BeamSectionSide side) const = 0;
};

}  // namespace fem
