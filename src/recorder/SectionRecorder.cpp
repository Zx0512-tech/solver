#include "fem/recorder/SectionRecorder.hpp"

#include <stdexcept>

namespace fem {
namespace {

Eigen::Matrix<double, 6, 1> forceVector(
    const BeamSectionForces& f) {
  Eigen::Matrix<double, 6, 1> result;
  result << f.N, f.Vy, f.Vz, f.T, f.My, f.Mz;
  return result;
}

void validateNewmarkDisplacement(
    const NewmarkResult& result) {
  if (result.displacement.rows() !=
          static_cast<Eigen::Index>(result.dofs.size()) ||
      result.displacement.cols() !=
          static_cast<Eigen::Index>(result.time.size())) {
    throw std::invalid_argument(
        "Newmark displacement history must match DOF count and time size");
  }
}

}  // namespace

double SectionResponseHistory::forceAt(
    std::size_t step,
    BeamSectionForceComponent component) const {
  if (forces.rows() != 6 ||
      forces.cols() != static_cast<Eigen::Index>(time.size())) {
    throw std::invalid_argument(
        "Section force history must be 6 x time.size()");
  }
  if (step >= time.size()) {
    throw std::out_of_range("Section response step is out of range");
  }

  return forces(
      static_cast<Eigen::Index>(
          static_cast<std::size_t>(component)),
      static_cast<Eigen::Index>(step));
}

ScalarHistory SectionResponseHistory::series(
    BeamSectionForceComponent component) const {
  if (forces.rows() != 6 ||
      forces.cols() != static_cast<Eigen::Index>(time.size())) {
    throw std::invalid_argument(
        "Section force history must be 6 x time.size()");
  }

  ScalarHistory result;
  result.time = time;
  result.value.resize(time.size());

  const Eigen::Index row =
      static_cast<Eigen::Index>(
          static_cast<std::size_t>(component));
  for (std::size_t step = 0; step < time.size(); ++step) {
    result.value[step] =
        forces(row, static_cast<Eigen::Index>(step));
  }
  return result;
}

SectionResponse SectionRecorder::record(
    const Model& model,
    ElementId element_id,
    double x,
    const StaticResult& result,
    BeamSectionSide side) const {
  return {
      element_id,
      x,
      side,
      model.beamSectionForces(
          element_id,
          x,
          result.displacement,
          0.0,
          side)};
}

SectionResponseHistory SectionRecorder::record(
    const Model& model,
    ElementId element_id,
    double x,
    const NewmarkResult& result,
    BeamSectionSide side) const {
  validateNewmarkDisplacement(result);

  SectionResponseHistory history;
  history.element_id = element_id;
  history.x = x;
  history.side = side;
  history.time = result.time;
  history.forces =
      Eigen::MatrixXd::Zero(
          6,
          static_cast<Eigen::Index>(result.time.size()));

  for (std::size_t step = 0; step < result.time.size(); ++step) {
    const BeamSectionForces f =
        model.beamSectionForces(
            element_id,
            x,
            result.displacement.col(
                static_cast<Eigen::Index>(step)),
            result.time[step],
            side);
    history.forces.col(static_cast<Eigen::Index>(step)) =
        forceVector(f);
  }

  return history;
}

ScalarHistory SectionRecorder::recordNormalStress(
    const Model& model,
    ElementId element_id,
    double x,
    double y,
    double z,
    const NewmarkResult& result,
    BeamSectionSide side) const {
  validateNewmarkDisplacement(result);

  ScalarHistory history;
  history.time = result.time;
  history.value.resize(result.time.size());

  for (std::size_t step = 0; step < result.time.size(); ++step) {
    history.value[step] =
        model.beamNormalStressAt(
            element_id,
            x,
            y,
            z,
            result.displacement.col(
                static_cast<Eigen::Index>(step)),
            result.time[step],
            side);
  }

  return history;
}

}  // namespace fem
