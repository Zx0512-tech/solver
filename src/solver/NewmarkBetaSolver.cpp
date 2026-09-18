#include "fem/solver/NewmarkBetaSolver.hpp"

#include <Eigen/Cholesky>

#include <cmath>
#include <stdexcept>\n#include <string>
#include <vector>

namespace fem {
namespace {

std::vector<Eigen::Index> freeDofs(const Model& model, const DofManager& dofs) {
  std::vector<Eigen::Index> result;
  result.reserve(dofs.size());

  for (const NodeId node_id : dofs.nodeOrder()) {
    const Node& node = model.node(node_id);
    for (std::size_t offset = 0; offset < kDofsPerFrameNode; ++offset) {
      const Dof dof = dofFromOffset(offset);
      if (!node.isFixed(dof)) {
        result.push_back(static_cast<Eigen::Index>(dofs.equation(node_id, dof)));
      }
    }
  }
  return result;
}

Eigen::MatrixXd selectMatrix(const Eigen::MatrixXd& matrix,
                             const std::vector<Eigen::Index>& indices) {
  const Eigen::Index n = static_cast<Eigen::Index>(indices.size());
  Eigen::MatrixXd reduced(n, n);
  for (Eigen::Index i = 0; i < n; ++i) {
    for (Eigen::Index j = 0; j < n; ++j) {
      reduced(i, j) = matrix(indices[static_cast<std::size_t>(i)],
                             indices[static_cast<std::size_t>(j)]);
    }
  }
  return reduced;
}

Eigen::VectorXd selectVector(const Eigen::VectorXd& vector,
                             const std::vector<Eigen::Index>& indices) {
  Eigen::VectorXd reduced(static_cast<Eigen::Index>(indices.size()));
  for (std::size_t i = 0; i < indices.size(); ++i) {
    reduced[static_cast<Eigen::Index>(i)] = vector[indices[i]];
  }
  return reduced;
}

Eigen::VectorXd initialVector(const Eigen::VectorXd& supplied,
                              Eigen::Index full_size,
                              const char* name) {
  if (supplied.size() == 0) {
    return Eigen::VectorXd::Zero(full_size);
  }
  if (supplied.size() != full_size) {
    throw std::invalid_argument(std::string(name) + " must be empty or match the global DOF count");
  }
  return supplied;
}

Eigen::VectorXd forceAt(double time,
                        const AssembledSystem& system,
                        const std::vector<NodalTimeLoad>& time_loads) {
  Eigen::VectorXd force = system.load;
  for (const auto& load : time_loads) {
    if (!load.value) {
      throw std::invalid_argument("Nodal time load function cannot be empty");
    }
    const std::size_t eq = system.dofs.equation(load.node_id, load.dof);
    force[static_cast<Eigen::Index>(eq)] += load.value(time);
  }
  return force;
}

}  // namespace

NewmarkResult NewmarkBetaSolver::solve(const Model& model,
                                       const NewmarkSettings& settings,
                                       const RayleighDamping& damping,
                                       const std::vector<NodalTimeLoad>& time_loads,
                                       const DynamicInitialState& initial) const {
  if (settings.time_step <= 0.0) {
    throw std::invalid_argument("Newmark time step must be positive");
  }
  if (settings.beta <= 0.0 || settings.gamma <= 0.0) {
    throw std::invalid_argument("Newmark beta and gamma must be positive");
  }
  if (settings.step_count == 0U) {
    throw std::invalid_argument("Newmark step count must be at least one");
  }

  AssembledSystem system = model.assemble();
  const auto free = freeDofs(model, system.dofs);
  if (free.empty()) {
    throw std::runtime_error("Model has no free degrees of freedom");
  }

  const Eigen::MatrixXd k = selectMatrix(system.stiffness, free);
  const Eigen::MatrixXd m = selectMatrix(system.mass, free);
  const Eigen::MatrixXd c = damping.matrix(m, k);

  Eigen::LLT<Eigen::MatrixXd> mass_factor(m);
  if (mass_factor.info() != Eigen::Success) {
    throw std::runtime_error(
        "Reduced mass matrix is not positive definite; check material density and constraints");
  }

  const Eigen::Index full_size = system.load.size();
  const Eigen::VectorXd u0_full = initialVector(initial.displacement, full_size, "Initial displacement");
  const Eigen::VectorXd v0_full = initialVector(initial.velocity, full_size, "Initial velocity");
  Eigen::VectorXd u = selectVector(u0_full, free);
  Eigen::VectorXd v = selectVector(v0_full, free);

  const Eigen::VectorXd f0 = selectVector(forceAt(0.0, system, time_loads), free);
  Eigen::VectorXd a = mass_factor.solve(f0 - c * v - k * u);

  const double dt = settings.time_step;
  const double beta = settings.beta;
  const double gamma = settings.gamma;

  const double a0 = 1.0 / (beta * dt * dt);
  const double a1 = gamma / (beta * dt);
  const double a2 = 1.0 / (beta * dt);
  const double a3 = 1.0 / (2.0 * beta) - 1.0;
  const double a4 = gamma / beta - 1.0;
  const double a5 = dt * (gamma / (2.0 * beta) - 1.0);

  const Eigen::MatrixXd effective_k = k + a0 * m + a1 * c;
  Eigen::LDLT<Eigen::MatrixXd> effective_factor(effective_k);
  if (effective_factor.info() != Eigen::Success) {
    throw std::runtime_error("Newmark effective stiffness factorization failed");
  }

  const Eigen::Index columns = static_cast<Eigen::Index>(settings.step_count + 1U);
  Eigen::MatrixXd u_history = Eigen::MatrixXd::Zero(full_size, columns);
  Eigen::MatrixXd v_history = Eigen::MatrixXd::Zero(full_size, columns);
  Eigen::MatrixXd a_history = Eigen::MatrixXd::Zero(full_size, columns);
  std::vector<double> times(settings.step_count + 1U, 0.0);

  for (std::size_t i = 0; i < free.size(); ++i) {
    u_history(free[i], 0) = u[static_cast<Eigen::Index>(i)];
    v_history(free[i], 0) = v[static_cast<Eigen::Index>(i)];
    a_history(free[i], 0) = a[static_cast<Eigen::Index>(i)];
  }

  for (std::size_t step = 0; step < settings.step_count; ++step) {
    const double time = static_cast<double>(step + 1U) * dt;
    const Eigen::VectorXd force = selectVector(forceAt(time, system, time_loads), free);

    const Eigen::VectorXd effective_force =
        force +
        m * (a0 * u + a2 * v + a3 * a) +
        c * (a1 * u + a4 * v + a5 * a);

    const Eigen::VectorXd u_next = effective_factor.solve(effective_force);
    const Eigen::VectorXd a_next = a0 * (u_next - u) - a2 * v - a3 * a;
    const Eigen::VectorXd v_next = v + dt * ((1.0 - gamma) * a + gamma * a_next);

    u = u_next;
    v = v_next;
    a = a_next;
    times[step + 1U] = time;

    const Eigen::Index column = static_cast<Eigen::Index>(step + 1U);
    for (std::size_t i = 0; i < free.size(); ++i) {
      const Eigen::Index local = static_cast<Eigen::Index>(i);
      u_history(free[i], column) = u[local];
      v_history(free[i], column) = v[local];
      a_history(free[i], column) = a[local];
    }
  }

  return {std::move(times),
          std::move(u_history),
          std::move(v_history),
          std::move(a_history),
          std::move(system.dofs)};
}

}  // namespace fem
