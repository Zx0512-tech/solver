#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/solver/LinearStaticSolver.hpp"
#include "fem/solver/NewmarkBetaSolver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Options {
  std::vector<std::size_t> sizes{100, 500, 2000};
  std::size_t newmark_steps{50};
  std::size_t repeats{3};
  std::filesystem::path csv_path{"solver-benchmark.csv"};
};

std::vector<std::size_t> parseSizes(const std::string& text) {
  std::vector<std::size_t> sizes;
  std::stringstream stream(text);
  std::string token;
  while (std::getline(stream, token, ',')) {
    if (token.empty()) {
      continue;
    }
    const auto value =
        static_cast<std::size_t>(std::stoull(token));
    if (value == 0U) {
      throw std::invalid_argument(
          "Benchmark element counts must be positive");
    }
    sizes.push_back(value);
  }
  if (sizes.empty()) {
    throw std::invalid_argument(
        "Benchmark requires at least one element count");
  }
  return sizes;
}

Options parseOptions(int argc, char** argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const auto requireValue = [&](const char* name) -> std::string {
      if (i + 1 >= argc) {
        throw std::invalid_argument(
            std::string(name) + " requires a value");
      }
      return argv[++i];
    };

    if (arg == "--sizes") {
      options.sizes = parseSizes(requireValue("--sizes"));
    } else if (arg == "--steps") {
      options.newmark_steps =
          static_cast<std::size_t>(
              std::stoull(requireValue("--steps")));
      if (options.newmark_steps == 0U) {
        throw std::invalid_argument(
            "--steps must be positive");
      }
    } else if (arg == "--repeats") {
      options.repeats =
          static_cast<std::size_t>(
              std::stoull(requireValue("--repeats")));
      if (options.repeats == 0U) {
        throw std::invalid_argument(
            "--repeats must be positive");
      }
    } else if (arg == "--csv") {
      options.csv_path = requireValue("--csv");
    } else {
      throw std::invalid_argument(
          "Unknown benchmark option: " + arg);
    }
  }
  return options;
}

void constrainAxialOnly(fem::Node& node) {
  node.fix(fem::Dof::UY);
  node.fix(fem::Dof::UZ);
  node.fix(fem::Dof::RX);
  node.fix(fem::Dof::RY);
  node.fix(fem::Dof::RZ);
}

fem::Model makeAxialChain(std::size_t element_count) {
  constexpr double element_length = 1.0;
  const fem::LinearElasticMaterial material(
      210.0e9, 0.3, 7850.0);
  const fem::BeamSection section(
      0.01, 2.0e-5, 3.0e-5, 4.0e-5);

  fem::Model model;
  for (std::size_t i = 0; i <= element_count; ++i) {
    auto& node = model.addNode(
        static_cast<fem::NodeId>(i + 1U),
        {static_cast<double>(i) * element_length, 0.0, 0.0});
    constrainAxialOnly(node);
  }
  model.node(1).fix(fem::Dof::UX);

  for (std::size_t i = 0; i < element_count; ++i) {
    model.addElement<fem::Beam3D>(
        static_cast<fem::ElementId>(i + 1U),
        static_cast<fem::NodeId>(i + 1U),
        static_cast<fem::NodeId>(i + 2U),
        material,
        section);
  }

  model.node(
      static_cast<fem::NodeId>(element_count + 1U))
      .addLoad(fem::Dof::UX, 50000.0);
  return model;
}

template <typename Function>
double medianMilliseconds(std::size_t repeats, Function&& function) {
  std::vector<double> samples;
  samples.reserve(repeats);

  for (std::size_t repeat = 0; repeat < repeats; ++repeat) {
    const auto start = Clock::now();
    function();
    const auto stop = Clock::now();
    samples.push_back(
        std::chrono::duration<double, std::milli>(
            stop - start).count());
  }

  std::sort(samples.begin(), samples.end());
  const std::size_t middle = samples.size() / 2U;
  if (samples.size() % 2U != 0U) {
    return samples[middle];
  }
  return 0.5 * (samples[middle - 1U] + samples[middle]);
}

std::size_t approximateSparseBytes(
    const Eigen::SparseMatrix<double>& matrix) {
  using StorageIndex =
      typename Eigen::SparseMatrix<double>::StorageIndex;
  return
      static_cast<std::size_t>(matrix.nonZeros()) *
          (sizeof(double) + sizeof(StorageIndex)) +
      static_cast<std::size_t>(matrix.outerSize() + 1) *
          sizeof(StorageIndex);
}

struct BenchmarkRow {
  std::size_t elements{};
  std::size_t global_dofs{};
  std::size_t free_dofs{};
  std::size_t k_nnz{};
  std::size_t m_nnz{};
  double k_fill_ratio{};
  std::size_t sparse_k_bytes{};
  std::size_t dense_k_bytes{};
  double dense_to_sparse_k_memory_ratio{};
  double assembly_ms{};
  double static_total_ms{};
  double newmark_total_ms{};
  std::size_t newmark_steps{};
  double static_tip_relative_error{};
};

BenchmarkRow runCase(
    std::size_t elements,
    std::size_t steps,
    std::size_t repeats) {
  fem::Model model = makeAxialChain(elements);

  fem::AssembledSystem assembled = model.assemble();
  volatile std::size_t nnz_sink =
      static_cast<std::size_t>(assembled.stiffness.nonZeros());
  (void)nnz_sink;

  const double assembly_ms =
      medianMilliseconds(repeats, [&] {
        auto system = model.assemble();
        volatile std::size_t sink =
            static_cast<std::size_t>(
                system.stiffness.nonZeros() +
                system.mass.nonZeros());
        (void)sink;
      });

  fem::StaticResult static_result =
      fem::LinearStaticSolver{}.solve(model);
  volatile double response_sink =
      static_result.displacementAt(
          static_cast<fem::NodeId>(elements + 1U),
          fem::Dof::UX);
  (void)response_sink;

  const double static_ms =
      medianMilliseconds(repeats, [&] {
        const auto result =
            fem::LinearStaticSolver{}.solve(model);
        volatile double sink =
            result.displacementAt(
                static_cast<fem::NodeId>(elements + 1U),
                fem::Dof::UX);
        (void)sink;
      });

  fem::NewmarkSettings settings;
  settings.time_step = 0.001;
  settings.step_count = steps;

  const std::vector<fem::NodalTimeLoad> time_loads = {
      {static_cast<fem::NodeId>(elements + 1U),
       fem::Dof::UX,
       [](double time) {
         return time > 0.0 ? 1000.0 : 0.0;
       }}};

  const double newmark_ms =
      medianMilliseconds(repeats, [&] {
        const auto result =
            fem::NewmarkBetaSolver{}.solve(
                model, settings, {}, time_loads);
        volatile double sink =
            result.displacementAt(
                settings.step_count,
                static_cast<fem::NodeId>(elements + 1U),
                fem::Dof::UX);
        (void)sink;
      });

  const std::size_t global_dofs =
      static_cast<std::size_t>(assembled.dofs.size());
  const std::size_t free_dofs = elements;
  const std::size_t k_nnz =
      static_cast<std::size_t>(
          assembled.stiffness.nonZeros());
  const std::size_t m_nnz =
      static_cast<std::size_t>(
          assembled.mass.nonZeros());

  const long double dense_entries =
      static_cast<long double>(global_dofs) *
      static_cast<long double>(global_dofs);
  const double fill_ratio =
      dense_entries > 0.0L
          ? static_cast<double>(
                static_cast<long double>(k_nnz) /
                dense_entries)
          : 0.0;

  const std::size_t sparse_k_bytes =
      approximateSparseBytes(assembled.stiffness);
  const std::size_t dense_k_bytes =
      global_dofs * global_dofs * sizeof(double);
  const double memory_ratio =
      sparse_k_bytes > 0U
          ? static_cast<double>(dense_k_bytes) /
                static_cast<double>(sparse_k_bytes)
          : 0.0;

  constexpr double e = 210.0e9;
  constexpr double area = 0.01;
  constexpr double load = 50000.0;
  const double expected_tip =
      load * static_cast<double>(elements) / (e * area);
  const double actual_tip =
      static_result.displacementAt(
          static_cast<fem::NodeId>(elements + 1U),
          fem::Dof::UX);
  const double relative_error =
      std::abs(actual_tip - expected_tip) /
      std::abs(expected_tip);

  if (!std::isfinite(newmark_ms) ||
      !std::isfinite(static_ms) ||
      relative_error > 1.0e-9) {
    throw std::runtime_error(
        "Benchmark correctness check failed");
  }

  return {
      elements,
      global_dofs,
      free_dofs,
      k_nnz,
      m_nnz,
      fill_ratio,
      sparse_k_bytes,
      dense_k_bytes,
      memory_ratio,
      assembly_ms,
      static_ms,
      newmark_ms,
      steps,
      relative_error};
}

void writeHeader(std::ostream& stream) {
  stream
      << "elements,global_dofs,free_dofs,k_nnz,m_nnz,"
      << "k_fill_ratio,sparse_k_bytes,dense_k_bytes,"
      << "dense_to_sparse_k_memory_ratio,"
      << "assembly_ms,static_total_ms,newmark_total_ms,"
      << "newmark_steps,static_tip_relative_error\n";
}

void writeRow(std::ostream& stream, const BenchmarkRow& row) {
  stream << std::setprecision(12)
         << row.elements << ','
         << row.global_dofs << ','
         << row.free_dofs << ','
         << row.k_nnz << ','
         << row.m_nnz << ','
         << row.k_fill_ratio << ','
         << row.sparse_k_bytes << ','
         << row.dense_k_bytes << ','
         << row.dense_to_sparse_k_memory_ratio << ','
         << row.assembly_ms << ','
         << row.static_total_ms << ','
         << row.newmark_total_ms << ','
         << row.newmark_steps << ','
         << row.static_tip_relative_error << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parseOptions(argc, argv);

    if (options.csv_path.has_parent_path()) {
      std::filesystem::create_directories(
          options.csv_path.parent_path());
    }
    std::ofstream csv(options.csv_path);
    if (!csv) {
      throw std::runtime_error(
          "Could not open benchmark CSV output");
    }

    writeHeader(std::cout);
    writeHeader(csv);

    for (const std::size_t elements : options.sizes) {
      const BenchmarkRow row =
          runCase(
              elements,
              options.newmark_steps,
              options.repeats);
      writeRow(std::cout, row);
      writeRow(csv, row);
    }

    std::cout
        << "Benchmark CSV: "
        << options.csv_path.string() << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "Benchmark failed: "
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
