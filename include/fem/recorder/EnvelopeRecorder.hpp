#pragma once

#include "fem/recorder/NodeRecorder.hpp"
#include "fem/recorder/ResultHistory.hpp"
#include "fem/recorder/SectionRecorder.hpp"

#include <cstddef>

namespace fem {

struct EnvelopePoint {
  double value{};
  double time{};
  std::size_t step{};
};

struct AbsoluteEnvelopePoint {
  double magnitude{};
  double value{};
  double time{};
  std::size_t step{};
};

struct EnvelopeResult {
  EnvelopePoint minimum;
  EnvelopePoint maximum;
  AbsoluteEnvelopePoint maximum_absolute;
};

class EnvelopeRecorder {
 public:
  EnvelopeResult record(
      const ScalarHistory& history) const;

  EnvelopeResult record(
      const NodeResponseHistory& history,
      NodeResponseQuantity quantity,
      Dof dof) const;

  EnvelopeResult record(
      const SectionResponseHistory& history,
      BeamSectionForceComponent component) const;
};

}  // namespace fem
