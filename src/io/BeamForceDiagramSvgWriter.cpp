#include "fem/io/BeamForceDiagramSvgWriter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem {
namespace {

using ValueGetter = std::function<double(const BeamSectionForces&)>;

struct SeriesDefinition {
  std::string label;
  ValueGetter value;
  bool dashed{};
};

void writeChart(
    std::ostream& stream,
    const BeamForceSeries& samples,
    const std::string& title,
    const std::vector<SeriesDefinition>& series) {
  if (samples.empty()) {
    throw std::invalid_argument("Cannot draw an empty Beam3D force series");
  }

  constexpr double width = 900.0;
  constexpr double height = 420.0;
  constexpr double left = 80.0;
  constexpr double right = 30.0;
  constexpr double top = 55.0;
  constexpr double bottom = 65.0;

  double x_min = samples.front().forces.x;
  double x_max = samples.front().forces.x;
  double y_min = 0.0;
  double y_max = 0.0;

  for (const auto& sample : samples) {
    x_min = std::min(x_min, sample.forces.x);
    x_max = std::max(x_max, sample.forces.x);
    for (const auto& definition : series) {
      const double value = definition.value(sample.forces);
      y_min = std::min(y_min, value);
      y_max = std::max(y_max, value);
    }
  }

  if (std::abs(x_max - x_min) < 1.0e-14) {
    x_max = x_min + 1.0;
  }
  if (std::abs(y_max - y_min) < 1.0e-14) {
    y_min -= 1.0;
    y_max += 1.0;
  } else {
    const double pad = 0.08 * (y_max - y_min);
    y_min -= pad;
    y_max += pad;
  }

  const double plot_width = width - left - right;
  const double plot_height = height - top - bottom;
  const auto px = [&](double x) {
    return left + (x - x_min) / (x_max - x_min) * plot_width;
  };
  const auto py = [&](double y) {
    return top + (y_max - y) / (y_max - y_min) * plot_height;
  };

  stream << std::setprecision(10);
  stream << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 "
         << width << ' ' << height << "\">\n";
  stream << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
  stream << "<text x=\"" << width / 2.0
         << "\" y=\"28\" text-anchor=\"middle\" "
            "font-family=\"sans-serif\" font-size=\"18\">"
         << title << "</text>\n";

  stream << "<line x1=\"" << left << "\" y1=\"" << top
         << "\" x2=\"" << left << "\" y2=\"" << (height - bottom)
         << "\" stroke=\"black\"/>\n";
  stream << "<line x1=\"" << left << "\" y1=\"" << (height - bottom)
         << "\" x2=\"" << (width - right) << "\" y2=\"" << (height - bottom)
         << "\" stroke=\"black\"/>\n";

  if (y_min <= 0.0 && y_max >= 0.0) {
    const double y0 = py(0.0);
    stream << "<line x1=\"" << left << "\" y1=\"" << y0
           << "\" x2=\"" << (width - right) << "\" y2=\"" << y0
           << "\" stroke=\"gray\" stroke-dasharray=\"4 4\"/>\n";
  }

  stream << "<text x=\"" << width / 2.0 << "\" y=\"" << (height - 18.0)
         << "\" text-anchor=\"middle\" font-family=\"sans-serif\" "
            "font-size=\"13\">x</text>\n";

  for (std::size_t series_index = 0; series_index < series.size(); ++series_index) {
    const auto& definition = series[series_index];
    stream << "<polyline fill=\"none\" stroke=\"black\" stroke-width=\"2\"";
    if (definition.dashed) {
      stream << " stroke-dasharray=\"8 5\"";
    }
    stream << " points=\"";
    for (const auto& sample : samples) {
      stream << px(sample.forces.x) << ',' << py(definition.value(sample.forces))
             << ' ';
    }
    stream << "\"/>\n";

    stream << "<text x=\"" << (left + 10.0 + 120.0 * static_cast<double>(series_index))
           << "\" y=\"" << (top - 15.0)
           << "\" font-family=\"sans-serif\" font-size=\"12\">"
           << definition.label << "</text>\n";
  }

  stream << "<text x=\"" << left << "\" y=\"" << (height - bottom + 20.0)
         << "\" font-family=\"sans-serif\" font-size=\"11\">"
         << x_min << "</text>\n";
  stream << "<text x=\"" << (width - right) << "\" y=\"" << (height - bottom + 20.0)
         << "\" text-anchor=\"end\" font-family=\"sans-serif\" font-size=\"11\">"
         << x_max << "</text>\n";
  stream << "</svg>\n";

  if (!stream) {
    throw std::runtime_error("Failed while writing Beam3D force SVG");
  }
}

template <typename Writer>
void writeSvgFile(
    const std::filesystem::path& path,
    Writer&& writer) {
  std::ofstream stream(path);
  if (!stream) {
    throw std::runtime_error("Could not open Beam3D force SVG output file");
  }
  writer(stream);
}

}  // namespace

void BeamForceDiagramSvgWriter::writeAxial(
    std::ostream& stream,
    const BeamForceSeries& samples) const {
  writeChart(
      stream,
      samples,
      "Axial force N",
      {{"N", [](const BeamSectionForces& f) { return f.N; }, false}});
}

void BeamForceDiagramSvgWriter::writeShear(
    std::ostream& stream,
    const BeamForceSeries& samples) const {
  writeChart(
      stream,
      samples,
      "Shear forces Vy / Vz",
      {
          {"Vy", [](const BeamSectionForces& f) { return f.Vy; }, false},
          {"Vz", [](const BeamSectionForces& f) { return f.Vz; }, true},
      });
}

void BeamForceDiagramSvgWriter::writeBending(
    std::ostream& stream,
    const BeamForceSeries& samples) const {
  writeChart(
      stream,
      samples,
      "Bending moments My / Mz",
      {
          {"My", [](const BeamSectionForces& f) { return f.My; }, false},
          {"Mz", [](const BeamSectionForces& f) { return f.Mz; }, true},
      });
}

void BeamForceDiagramSvgWriter::writeTorsion(
    std::ostream& stream,
    const BeamForceSeries& samples) const {
  writeChart(
      stream,
      samples,
      "Torsion T",
      {{"T", [](const BeamSectionForces& f) { return f.T; }, false}});
}

void BeamForceDiagramSvgWriter::writeSet(
    const std::filesystem::path& directory,
    const std::string& base_name,
    const BeamForceSeries& samples) const {
  if (base_name.empty()) {
    throw std::invalid_argument("Beam3D diagram base name cannot be empty");
  }
  std::filesystem::create_directories(directory);

  writeSvgFile(
      directory / (base_name + "_axial.svg"),
      [&](std::ostream& stream) { writeAxial(stream, samples); });
  writeSvgFile(
      directory / (base_name + "_shear.svg"),
      [&](std::ostream& stream) { writeShear(stream, samples); });
  writeSvgFile(
      directory / (base_name + "_bending.svg"),
      [&](std::ostream& stream) { writeBending(stream, samples); });
  writeSvgFile(
      directory / (base_name + "_torsion.svg"),
      [&](std::ostream& stream) { writeTorsion(stream, samples); });
}

}  // namespace fem
