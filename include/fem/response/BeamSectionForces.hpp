#pragma once

namespace fem {

// Section forces act on the positive-local-x cut face.
// N > 0 is tension. Point loads are discontinuities, so callers can request
// the limit immediately to the left or right of the load location.
enum class BeamSectionSide {
  Left,
  Right,
};

struct BeamSectionForces {
  double x{};
  double N{};
  double Vy{};
  double Vz{};
  double T{};
  double My{};
  double Mz{};
};

}  // namespace fem
