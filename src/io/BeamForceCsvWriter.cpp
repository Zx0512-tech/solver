#include "fem/io/BeamForceCsvWriter.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>

namespace fem {
namespace {

const char* sideName(BeamSectionSide side) {
  return side == BeamSectionSide::Left ? "left" : "right";
}

}  // namespace

void BeamForceCsvWriter::write(
    std::ostream& stream,
    const BeamForceSeries& samples) const {
  stream << "x,side,N,Vy,Vz,T,My,Mz\n";
  stream << std::setprecision(17);

  for (const auto& sample : samples) {
    const auto& f = sample.forces;
    stream << f.x << ',' << sideName(sample.side) << ','
           << f.N << ',' << f.Vy << ',' << f.Vz << ','
           << f.T << ',' << f.My << ',' << f.Mz << '\n';
  }

  if (!stream) {
    throw std::runtime_error("Failed while writing Beam3D force CSV");
  }
}

void BeamForceCsvWriter::writeFile(
    const std::filesystem::path& path,
    const BeamForceSeries& samples) const {
  if (path.has_parent_path()) {
    std::filesystem::create_directories(path.parent_path());
  }
  std::ofstream stream(path);
  if (!stream) {
    throw std::runtime_error("Could not open Beam3D force CSV output file");
  }
  write(stream, samples);
}

}  // namespace fem
