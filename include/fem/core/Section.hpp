#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fem {

enum class BeamSectionKind {
  General,
  Rectangle,
  SolidCircle,
};

namespace detail {

inline void requirePositive(double value, const char* message) {
  if (value <= 0.0) {
    throw std::invalid_argument(message);
  }
}

inline double rectangleIy(double size_y, double size_z) {
  requirePositive(size_y, "RectangleSection size_y must be positive");
  requirePositive(size_z, "RectangleSection size_z must be positive");
  return size_y * size_z * size_z * size_z / 12.0;
}

inline double rectangleIz(double size_y, double size_z) {
  requirePositive(size_y, "RectangleSection size_y must be positive");
  requirePositive(size_z, "RectangleSection size_z must be positive");
  return size_z * size_y * size_y * size_y / 12.0;
}

inline double rectangleTorsionConstant(double size_y, double size_z) {
  requirePositive(size_y, "RectangleSection size_y must be positive");
  requirePositive(size_z, "RectangleSection size_z must be positive");

  const double a = std::max(size_y, size_z);
  const double b = std::min(size_y, size_z);
  const double ratio = b / a;
  const double ratio4 = ratio * ratio * ratio * ratio;

  // Standard engineering approximation for the Saint-Venant torsion
  // constant of a solid rectangle. a is the long side and b the short side.
  return a * b * b * b *
         (1.0 / 3.0 - 0.21 * ratio * (1.0 - ratio4 / 12.0));
}

inline double pi() {
  return std::acos(-1.0);
}

}  // namespace detail

struct BeamSection {
  double area{0.0};
  double iy{0.0};
  double iz{0.0};
  double torsion_constant{0.0};

  BeamSectionKind kind{BeamSectionKind::General};

  // Geometry metadata used for stress recovery when known.
  // Rectangle: full dimensions along local y and z.
  double size_y{0.0};
  double size_z{0.0};

  // SolidCircle: radius.
  double radius{0.0};

  BeamSection(double a, double iy_value, double iz_value, double j)
      : area(a), iy(iy_value), iz(iz_value), torsion_constant(j) {
    validateProperties();
  }

 protected:
  BeamSection(double a,
              double iy_value,
              double iz_value,
              double j,
              BeamSectionKind section_kind,
              double section_size_y,
              double section_size_z,
              double section_radius)
      : area(a),
        iy(iy_value),
        iz(iz_value),
        torsion_constant(j),
        kind(section_kind),
        size_y(section_size_y),
        size_z(section_size_z),
        radius(section_radius) {
    validateProperties();
  }

 private:
  void validateProperties() const {
    if (area <= 0.0 || iy <= 0.0 || iz <= 0.0 || torsion_constant <= 0.0) {
      throw std::invalid_argument(
          "Beam section properties A, Iy, Iz and J must be positive");
    }
  }
};

struct GeneralSection final : BeamSection {
  GeneralSection(double a, double iy_value, double iz_value, double j)
      : BeamSection(a, iy_value, iz_value, j) {}
};

struct RectangleSection final : BeamSection {
  RectangleSection(double section_size_y, double section_size_z)
      : BeamSection(
            checkedArea(section_size_y, section_size_z),
            detail::rectangleIy(section_size_y, section_size_z),
            detail::rectangleIz(section_size_y, section_size_z),
            detail::rectangleTorsionConstant(section_size_y, section_size_z),
            BeamSectionKind::Rectangle,
            section_size_y,
            section_size_z,
            0.0) {}

 private:
  static double checkedArea(double section_size_y, double section_size_z) {
    detail::requirePositive(
        section_size_y, "RectangleSection size_y must be positive");
    detail::requirePositive(
        section_size_z, "RectangleSection size_z must be positive");
    return section_size_y * section_size_z;
  }
};

struct CircularSection final : BeamSection {
  explicit CircularSection(double section_radius)
      : BeamSection(
            areaFromRadius(section_radius),
            bendingInertia(section_radius),
            bendingInertia(section_radius),
            polarInertia(section_radius),
            BeamSectionKind::SolidCircle,
            0.0,
            0.0,
            section_radius) {}

 private:
  static double checkedRadius(double section_radius) {
    detail::requirePositive(
        section_radius, "CircularSection radius must be positive");
    return section_radius;
  }

  static double areaFromRadius(double section_radius) {
    const double r = checkedRadius(section_radius);
    return detail::pi() * r * r;
  }

  static double bendingInertia(double section_radius) {
    const double r = checkedRadius(section_radius);
    const double r2 = r * r;
    return detail::pi() * r2 * r2 / 4.0;
  }

  static double polarInertia(double section_radius) {
    const double r = checkedRadius(section_radius);
    const double r2 = r * r;
    return detail::pi() * r2 * r2 / 2.0;
  }
};

}  // namespace fem
