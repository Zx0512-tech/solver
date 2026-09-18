#pragma once

namespace fem {

struct BeamSectionPoint {
  double y{};
  double z{};
};

struct BeamSectionStress {
  double sigma_x{};
  double tau_xy{};
  double tau_xz{};
};

struct BeamNormalStressExtrema {
  double sigma_min{};
  double sigma_max{};
  BeamSectionPoint min_point{};
  BeamSectionPoint max_point{};
};

}  // namespace fem
