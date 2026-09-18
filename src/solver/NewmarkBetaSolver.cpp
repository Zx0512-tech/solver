#include "fem/solver/NewmarkBetaSolver.hpp"
#include "fem/solver/SparseDofReducer.hpp"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem {
namespace {

std::vector<Eigen::Index> freeDofs(
    const Model& model,
    const DofManager& dofs) {
  std::vector<Eigen::Index> result;
  result.reserve(dofs.size());

  for (const NodeId node_id : dofs.nodeOrder()) {
    const Node& node = model.node(node_id);
    for (std::size_t offset = 0;
         offset < kDofsPerFrameNode;
         ++offset) {
      const Dof dof = dofFromOffset(offset);
      if (!node.isFixed(dof)) {
        result.push_back(
            static_cast<Eigen::Index>(
                dofs.equation(node_id, dof)));
      }
    }
  }
  return result;
}

Eigen::VectorXd initialVector(
    const Eigen::VectorXd& supplied,
    Eigen::Index full_size,
    const char* name) {
  if (supplied.size() == 0) {
    return Eigen::VectorXd::Zero(full_size);
  }
  if (supplied.size() != full_size) {
    throw std::invalid_argument(
        std::string(name) +
        " must be empty or match the global DOF count");
  }
  return supplied;
}

Eigen::VectorXd forceAt(
    double time,
    const Model& model,
    const AssembledSystem& system,
    const std::vector<NodalTimeLoad>& time_loads) {
  Eigen::VectorXd force = model.loadVector(time);
  if (force.size() != system.load.size()) {
    throw std::runtime_error(
        "Model load vector size changed during transient analysis");
  }

  for (const auto& load : time_loads) {
    if (!load.value) {
      throw std::invalid_argument(
          "Nodal time load function cannot be empty");
    }
    const std::size_t equation =
        system.dofs.equation(load.node_id, load.dof);
    force[static_cast<Eigen::Index>(equation)] +=
        load.value(time);
  }
  return force;
}

}  // namespace

NewmarkResult NewmarkBetaSolver::solve(
    const Model& model,
    const NewmarkSettings& settings,
    const RayleighDamping& damping,
    const std::vector<NodalTimeLoad>& time_loads,
    const DynamicInitialState& initial) const {
  if (settings.time_step <= 0.0) {
    throw std::invalid_argument(
        "Newmark time step must be positive");
  }
  if (settings.beta <= 0.0 || settings.gamma <= 0.0) {
    throw std::invalid_argument(
        "Newmark beta and gamma must be positive");
  }
  if (settings.step_count == 0U) {
    throw std::invalid_argument(
        "Newmark step count must be at least one");
  }

  AssembledSystem system = model.assemble();
  const auto free = freeDofs(model, system.dofs);
  if (free.empty()) {
    throw std::runtime_error(
        "Model has no free degrees of freedom");
  }

  const Eigen::SparseMatrix<double> stiffness =
      SparseDofReducer::matrix(system.stiffness, free);
  const Eigen::SparseMatrix<double> mass =
      SparseDofReducer::matrix(system.mass, free);
  const Eigen::SparseMatrix<double> damping_matrix =
      damping.matrix(mass, stiffness);

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> mass_factor;
  mass_factor.compute(mass);
  if (mass_factor.info() != Eigen::Success) {
    throw std::runtime_error(
        "Reduced mass matrix factorization failed; "
        "check material density and constraints");
  }

  const Eigen::Index full_size = system.load.size();
  const Eigen::VectorXd u0_full =
      initialVector(
          initial.displacement,
          full_size,
          "Initial displacement");
  const Eigen::VectorXd v0_full =
      initialVector(
          initial.velocity,
          full_size,
          "Initial velocity");

  Eigen::VectorXd displacement =
      SparseDofReducer::vector(u0_full, free);
  Eigen::VectorXd velocity =
      SparseDofReducer::vector(v0_full, free);

  const Eigen::VectorXd f0 =
      SparseDofReducer::vector(
          forceAt(0.0, model, system, time_loads),
          free);

  Eigen::VectorXd acceleration =
      mass_factor.solve(
          f0 -
          damping_matrix * velocity -
          stiffness * displacement);
  if (mass_factor.info() != Eigen::Success ||
      !acceleration.allFinite()) {
    throw std::runtime_error(
        "Initial acceleration solve failed");
  }

  const double dt = settings.time_step;
  const double beta = settings.beta;
  const double gamma = settings.gamma;

  const double a0 = 1.0 / (beta * dt * dt);
  const double a1 = gamma / (beta * dt);
  const double a2 = 1.0 / (beta * dt);
  const double a3 = 1.0 / (2.0 * beta) - 1.0;
  const double a4 = gamma / beta - 1.0;
  const double a5 =
      dt * (gamma / (2.0 * beta) - 1.0);

  Eigen::SparseMatrix<double> effective_stiffness =
      stiffness + a0 * mass + a1 * damping_matrix;
  effective_stiffness.prune(0.0);
  effective_stiffness.makeCompressed();

  // The linear system is time-invariant, so factorize Keff once and reuse
  // the factorization for every Newmark step.
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>
      effective_factor;
  effective_factor.compute(effective_stiffness);
  if (effective_factor.info() != Eigen::Success) {
    throw std::runtime_error(
        "Newmark effective stiffness factorization failed");
  }

  const Eigen::Index columns =
      static_cast<Eigen::Index>(
          settings.step_count + 1U);
  Eigen::MatrixXd displacement_history =
      Eigen::MatrixXd::Zero(full_size, columns);
  Eigen::MatrixXd velocity_history =
      Eigen::MatrixXd::Zero(full_size, columns);
  Eigen::MatrixXd acceleration_history =
      Eigen::MatrixXd::Zero(full_size, columns);
  std::vector<double> times(
      settings.step_count + 1U, 0.0);

  for (std::size_t i = 0; i < free.size(); ++i) {
    const Eigen::Index local =
        static_cast<Eigen::Index>(i);
    displacement_history(free[i], 0) =
        displacement[local];
    velocity_history(free[i], 0) =
        velocity[local];
    acceleration_history(free[i], 0) =
        acceleration[local];
  }

  for (std::size_t step = 0;
       step < settings.step_count;
       ++step) {
    const double time =
        static_cast<double>(step + 1U) * dt;

    const Eigen::VectorXd force =
        SparseDofReducer::vector(
            forceAt(time, model, system, time_loads),
            free);

    const Eigen::VectorXd effective_force =
        force +
        mass *
            (a0 * displacement +
             a2 * velocity +
             a3 * acceleration) +
        damping_matrix *
            (a1 * displacement +
             a4 * velocity +
             a5 * acceleration);

    const Eigen::VectorXd displacement_next =
        effective_factor.solve(effective_force);
    if (effective_factor.info() != Eigen::Success ||
        !displacement_next.allFinite()) {
      throw std::runtime_error(
          "Newmark effective stiffness solve failed");
    }

    const Eigen::VectorXd acceleration_next =
        a0 * (displacement_next - displacement) -
        a2 * velocity -
        a3 * acceleration;
    const Eigen::VectorXd velocity_next =
        velocity +
        dt * (
            (1.0 - gamma) * acceleration +
            gamma * acceleration_next);

    displacement = displacement_next;
    velocity = velocity_next;
    acceleration = acceleration_next;
    times[step + 1U] = time;

    const Eigen::Index column =
        static_cast<Eigen::Index>(step + 1U);
    for (std::size_t i = 0; i < free.size(); ++i) {
      const Eigen::Index local =
          static_cast<Eigen::Index>(i);
      displacement_history(free[i], column) =
          displacement[local];
      velocity_history(free[i], column) =
          velocity[local];
      acceleration_history(free[i], column) =
          acceleration[local];
    }
  }

  return {
      std::move(times),
      std::move(displacement_history),
      std::move(velocity_history),
      std::move(acceleration_history),
      std::move(system.dofs)};
}

}  // namespace fem
