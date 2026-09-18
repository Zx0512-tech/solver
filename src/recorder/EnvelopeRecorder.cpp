#include "fem/recorder/EnvelopeRecorder.hpp"

#include <cmath>
#include <stdexcept>

namespace fem {

EnvelopeResult EnvelopeRecorder::record(
    const ScalarHistory& history) const {
  if (history.time.empty()) {
    throw std::invalid_argument(
        "Cannot compute an envelope of an empty history");
  }
  if (history.time.size() != history.value.size()) {
    throw std::invalid_argument(
        "Envelope time and value history sizes must match");
  }

  EnvelopeResult result;
  result.minimum = {history.value[0], history.time[0], 0U};
  result.maximum = result.minimum;
  result.maximum_absolute = {
      std::abs(history.value[0]),
      history.value[0],
      history.time[0],
      0U};

  for (std::size_t step = 1; step < history.value.size(); ++step) {
    const double value = history.value[step];

    if (value < result.minimum.value) {
      result.minimum = {value, history.time[step], step};
    }
    if (value > result.maximum.value) {
      result.maximum = {value, history.time[step], step};
    }

    const double magnitude = std::abs(value);
    if (magnitude > result.maximum_absolute.magnitude) {
      result.maximum_absolute = {
          magnitude, value, history.time[step], step};
    }
  }

  return result;
}

EnvelopeResult EnvelopeRecorder::record(
    const NodeResponseHistory& history,
    NodeResponseQuantity quantity,
    Dof dof) const {
  return record(history.series(quantity, dof));
}

EnvelopeResult EnvelopeRecorder::record(
    const SectionResponseHistory& history,
    BeamSectionForceComponent component) const {
  return record(history.series(component));
}

}  // namespace fem
