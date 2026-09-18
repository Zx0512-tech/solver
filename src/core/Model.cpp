#include "fem/core/Model.hpp"

#include <stdexcept>

namespace fem {

Node& Model::addNode(NodeId id, const Eigen::Vector3d& coordinates) {
  auto [it, inserted] = nodes_.emplace(id, Node{id, coordinates});
  if (!inserted) {
    throw std::invalid_argument("Duplicate node id");
  }
  return it->second;
}

Node& Model::node(NodeId id) {
  const auto it = nodes_.find(id);
  if (it == nodes_.end()) {
    throw std::out_of_range("Unknown node id");
  }
  return it->second;
}

const Node& Model::node(NodeId id) const {
  const auto it = nodes_.find(id);
  if (it == nodes_.end()) {
    throw std::out_of_range("Unknown node id");
  }
  return it->second;
}

AssembledSystem Model::assemble() const {
  DofManager dm(nodes_);
  const Eigen::Index n = static_cast<Eigen::Index>(dm.size());

  Eigen::MatrixXd global_k = Eigen::MatrixXd::Zero(n, n);
  Eigen::MatrixXd global_m = Eigen::MatrixXd::Zero(n, n);
  Eigen::VectorXd global_f = Eigen::VectorXd::Zero(n);

  for (const auto& [node_id, node_ref] : nodes_) {
    for (std::size_t offset = 0; offset < kDofsPerFrameNode; ++offset) {
      const Dof dof = dofFromOffset(offset);
      global_f[static_cast<Eigen::Index>(dm.equation(node_id, dof))] += node_ref.load(dof);
    }
  }

  const NodeResolver resolver = [this](NodeId id) -> const Node& { return node(id); };

  for (const auto& [element_id, element] : elements_) {
    (void)element_id;
    const auto element_dofs = element->dofs();
    const Eigen::MatrixXd ke = element->stiffness(resolver);
    const Eigen::MatrixXd me = element->mass(resolver);
    const Eigen::Index ndof = static_cast<Eigen::Index>(element_dofs.size());

    if (ke.rows() != ndof || ke.cols() != ndof) {
      throw std::runtime_error("Element stiffness size does not match element DOF count");
    }
    if (me.rows() != ndof || me.cols() != ndof) {
      throw std::runtime_error("Element mass size does not match element DOF count");
    }

    for (Eigen::Index a = 0; a < ndof; ++a) {
      const auto [node_a, dof_a] = element_dofs[static_cast<std::size_t>(a)];
      const Eigen::Index ia = static_cast<Eigen::Index>(dm.equation(node_a, dof_a));

      for (Eigen::Index b = 0; b < ndof; ++b) {
        const auto [node_b, dof_b] = element_dofs[static_cast<std::size_t>(b)];
        const Eigen::Index ib = static_cast<Eigen::Index>(dm.equation(node_b, dof_b));
        global_k(ia, ib) += ke(a, b);
        global_m(ia, ib) += me(a, b);
      }
    }
  }

  return {std::move(global_k), std::move(global_m), std::move(global_f), std::move(dm)};
}

}  // namespace fem
