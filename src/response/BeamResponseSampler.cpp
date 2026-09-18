#include "fem/response/BeamResponseSampler.hpp"

#include "fem/core/Model.hpp"
#include "fem/elements/Beam3D.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace fem {
namespace {

struct SampleLocation {
  double x{};
  bool has_jump{};
};

void addOrMerge(
    std::vector<SampleLocation>& locations,
    double x,
    bool has_jump,
    double tolerance) {
  for (auto& location : locations) {
    if (std::abs(location.x - x) <= tolerance) {
      location.has_jump = location.has_jump || has_jump;
      return;
    }
  }
  locations.push_back({x, has_jump});
}

}  // namespace

BeamForceSeries BeamResponseSampler::sample(
    const Model& model,
    ElementId element_id,
    const Eigen::VectorXd& global_displacement,
    std::size_t point_count,
    double time) const {
  if (point_count < 2U) {
    throw std::invalid_argument(
        "BeamResponseSampler requires at least two uniform points");
  }

  const Element& target = model.element(element_id);
  const auto* beam = dynamic_cast<const Beam3D*>(&target);
  if (beam == nullptr) {
    throw std::invalid_argument(
        "BeamResponseSampler requires a Beam3D element");
  }

  const NodeResolver resolver = [&model](NodeId id) -> const Node& {
    return model.node(id);
  };
  const double length = beam->length(resolver);
  const double tolerance = 1.0e-10 * std::max(1.0, length);

  std::vector<SampleLocation> locations;
  locations.reserve(point_count + 8U);

  for (std::size_t i = 0; i < point_count; ++i) {
    const double ratio =
        static_cast<double>(i) / static_cast<double>(point_count - 1U);
    addOrMerge(locations, ratio * length, false, tolerance);
  }

  for (const auto& extra : model.elementLoadSampleLocations(element_id)) {
    if (extra.x < -tolerance || extra.x > length + tolerance) {
      throw std::runtime_error(
          "Element load supplied a sampling location outside the beam");
    }
    const double x = std::clamp(extra.x, 0.0, length);
    addOrMerge(locations, x, extra.has_jump, tolerance);
  }

  std::sort(
      locations.begin(), locations.end(),
      [](const SampleLocation& a, const SampleLocation& b) {
        return a.x < b.x;
      });

  BeamForceSeries result;
  result.reserve(locations.size() + 4U);
  for (const auto& location : locations) {
    const bool has_left_domain = location.x > tolerance;
    if (location.has_jump && has_left_domain) {
      result.push_back({
          model.beamSectionForces(
              element_id,
              location.x,
              global_displacement,
              time,
              BeamSectionSide::Left),
          BeamSectionSide::Left});
    }

    result.push_back({
        model.beamSectionForces(
            element_id,
            location.x,
            global_displacement,
            time,
            BeamSectionSide::Right),
        BeamSectionSide::Right});
  }

  return result;
}

}  // namespace fem
