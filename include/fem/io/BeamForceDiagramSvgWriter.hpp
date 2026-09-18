#pragma once

#include "fem/response/BeamResponseSampler.hpp"

#include <filesystem>
#include <iosfwd>
#include <string>

namespace fem {

class BeamForceDiagramSvgWriter {
 public:
  void writeAxial(
      std::ostream& stream,
      const BeamForceSeries& samples) const;
  void writeShear(
      std::ostream& stream,
      const BeamForceSeries& samples) const;
  void writeBending(
      std::ostream& stream,
      const BeamForceSeries& samples) const;
  void writeTorsion(
      std::ostream& stream,
      const BeamForceSeries& samples) const;

  void writeSet(
      const std::filesystem::path& directory,
      const std::string& base_name,
      const BeamForceSeries& samples) const;
};

}  // namespace fem
