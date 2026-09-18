#pragma once
#include <stdexcept>
namespace fem {
struct LinearElasticMaterial {
  double youngs_modulus{0.0};
  double poisson_ratio{0.0};
  double density{0.0};
  LinearElasticMaterial(double e, double nu, double rho=0.0) : youngs_modulus(e), poisson_ratio(nu), density(rho) {
    if (youngs_modulus <= 0.0) throw std::invalid_argument("Young's modulus must be positive");
    if (poisson_ratio <= -1.0 || poisson_ratio >= 0.5) throw std::invalid_argument("Poisson ratio must be in (-1, 0.5)");
    if (density < 0.0) throw std::invalid_argument("Density cannot be negative");
  }
  double shearModulus() const noexcept { return youngs_modulus / (2.0 * (1.0 + poisson_ratio)); }
};
}  // namespace fem
