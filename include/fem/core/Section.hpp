#pragma once
#include <stdexcept>
namespace fem {
struct BeamSection {
  double area{0.0}, iy{0.0}, iz{0.0}, torsion_constant{0.0};
  BeamSection(double a,double iy_value,double iz_value,double j) : area(a),iy(iy_value),iz(iz_value),torsion_constant(j) {
    if (area<=0.0 || iy<=0.0 || iz<=0.0 || torsion_constant<=0.0)
      throw std::invalid_argument("Beam section properties A, Iy, Iz and J must be positive");
  }
};
}  // namespace fem
