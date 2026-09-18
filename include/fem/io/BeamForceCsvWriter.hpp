#pragma once

#include "fem/response/BeamResponseSampler.hpp"

#include <filesystem>
#include <iosfwd>

namespace fem {

class BeamForceCsvWriter {
 public:
  void write(
      std::ostream& stream,
      const BeamForceSeries& samples) const;

  void writeFile(
      const std::filesystem::path& path,
      const BeamForceSeries& samples) const;
};

}  // namespace fem
