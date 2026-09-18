#include "fem/response/BeamSectionStressRecovery.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace fem {
namespace {

double forceTolerance(const BeamSectionForces& f) {
  const double scale = std::max({
      1.0,
      std::abs(f.N),
      std::abs(f.Vy),
      std::abs(f.Vz),
      std::abs(f.T),
      std::abs(f.My),
      std::abs(f.Mz)});
  return 1.0e-12 * scale;
}

void validatePoint(
    const BeamSection& section,
    double y,
    double z) {
  switch (section.kind) {
    case BeamSectionKind::General:
      return;

    case BeamSectionKind::Rectangle: {
      const double tolerance =
          1.0e-12 * std::max({1.0, section.size_y, section.size_z});
      if (std::abs(y) > section.size_y / 2.0 + tolerance ||
          std::abs(z) > section.size_z / 2.0 + tolerance) {
        throw std::invalid_argument(
            "Stress recovery point lies outside RectangleSection");
      }
      return;
    }

    case BeamSectionKind::SolidCircle: {
      const double tolerance =
          1.0e-12 * std::max(1.0, section.radius);
      const double limit = section.radius + tolerance;
      if (y * y + z * z > limit * limit) {
        throw std::invalid_argument(
            "Stress recovery point lies outside CircularSection");
      }
      return;
    }
  }

  throw std::logic_error("Unknown BeamSectionKind");
}

BeamSectionPoint rectanglePoint(
    const BeamSection& section,
    bool positive_y,
    bool positive_z) {
  return {
      (positive_y ? 0.5 : -0.5) * section.size_y,
      (positive_z ? 0.5 : -0.5) * section.size_z};
}

}  // namespace

double BeamSectionStressRecovery::normalStressAt(
    const BeamSection& section,
    const BeamSectionForces& forces,
    double y,
    double z) const {
  validatePoint(section, y, z);

  return forces.N / section.area -
         forces.Mz * y / section.iz +
         forces.My * z / section.iy;
}

BeamSectionStress BeamSectionStressRecovery::stressAt(
    const BeamSection& section,
    const BeamSectionForces& forces,
    double y,
    double z) const {
  BeamSectionStress stress;
  stress.sigma_x = normalStressAt(section, forces, y, z);

  const double tolerance = forceTolerance(forces);

  switch (section.kind) {
    case BeamSectionKind::General:
      if (std::abs(forces.Vy) > tolerance ||
          std::abs(forces.Vz) > tolerance ||
          std::abs(forces.T) > tolerance) {
        throw std::logic_error(
            "GeneralSection full shear/torsion stress recovery is not implemented; "
            "normalStressAt remains available");
      }
      return stress;

    case BeamSectionKind::Rectangle: {
      if (std::abs(forces.T) > tolerance) {
        throw std::logic_error(
            "RectangleSection torsional point stress is not implemented");
      }

      const double eta_y = 2.0 * y / section.size_y;
      const double eta_z = 2.0 * z / section.size_z;

      stress.tau_xy =
          1.5 * forces.Vy / section.area *
          std::max(0.0, 1.0 - eta_y * eta_y);
      stress.tau_xz =
          1.5 * forces.Vz / section.area *
          std::max(0.0, 1.0 - eta_z * eta_z);
      return stress;
    }

    case BeamSectionKind::SolidCircle:
      if (std::abs(forces.Vy) > tolerance ||
          std::abs(forces.Vz) > tolerance) {
        throw std::logic_error(
            "CircularSection transverse-shear point stress is not implemented");
      }

      stress.tau_xy =
          -forces.T * z / section.torsion_constant;
      stress.tau_xz =
          forces.T * y / section.torsion_constant;
      return stress;
  }

  throw std::logic_error("Unknown BeamSectionKind");
}

BeamNormalStressExtrema BeamSectionStressRecovery::normalStressExtrema(
    const BeamSection& section,
    const BeamSectionForces& forces) const {
  switch (section.kind) {
    case BeamSectionKind::General:
      throw std::logic_error(
          "GeneralSection has no boundary geometry for normal-stress extrema");

    case BeamSectionKind::Rectangle: {
      const std::array<BeamSectionPoint, 4> points = {
          rectanglePoint(section, false, false),
          rectanglePoint(section, false, true),
          rectanglePoint(section, true, false),
          rectanglePoint(section, true, true)};

      BeamNormalStressExtrema result;
      result.min_point = points[0];
      result.max_point = points[0];
      result.sigma_min =
          normalStressAt(section, forces, points[0].y, points[0].z);
      result.sigma_max = result.sigma_min;

      for (std::size_t i = 1; i < points.size(); ++i) {
        const double sigma =
            normalStressAt(section, forces, points[i].y, points[i].z);
        if (sigma < result.sigma_min) {
          result.sigma_min = sigma;
          result.min_point = points[i];
        }
        if (sigma > result.sigma_max) {
          result.sigma_max = sigma;
          result.max_point = points[i];
        }
      }
      return result;
    }

    case BeamSectionKind::SolidCircle: {
      BeamNormalStressExtrema result;
      const double base = forces.N / section.area;
      const double gradient_y = -forces.Mz / section.iz;
      const double gradient_z = forces.My / section.iy;
      const double norm =
          std::sqrt(
              gradient_y * gradient_y +
              gradient_z * gradient_z);

      if (norm <= 1.0e-30) {
        result.sigma_min = base;
        result.sigma_max = base;
        return result;
      }

      const double y =
          section.radius * gradient_y / norm;
      const double z =
          section.radius * gradient_z / norm;

      result.max_point = {y, z};
      result.min_point = {-y, -z};
      result.sigma_max = base + section.radius * norm;
      result.sigma_min = base - section.radius * norm;
      return result;
    }
  }

  throw std::logic_error("Unknown BeamSectionKind");
}

}  // namespace fem
