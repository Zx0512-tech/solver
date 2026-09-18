#include "fem/core/Model.hpp"

#include "fem/elements/Beam3D.hpp"
#include "fem/loads/BeamElementLoad3D.hpp"

#include <Eigen/Geometry>

#include <stdexcept>

namespace fem {
namespace {

BeamLoadResultant3D scaledResultant(
    BeamLoadResultant3D result,
    double factor) {
  result.force *= factor;
  result.moment *= factor;
  return result;
}

BeamLoadResultant3D beamLoadResultantTo(
    const ElementLoad& load,
    const Beam3D& beam,
    const NodeResolver& node,
    double x,
    double time,
    BeamSectionSide side) {
  if (const auto* timed =
          dynamic_cast<const TimeDependentElementLoad*>(&load)) {
    return scaledResultant(
        beamLoadResultantTo(
            timed->spatialLoad(), beam, node, x, time, side),
        timed->scaleAt(time));
  }

  const auto* beam_load = dynamic_cast<const BeamElementLoad3D*>(&load);
  if (beam_load == nullptr) {
    throw std::runtime_error(
        "Beam section-force recovery encountered an unsupported element load");
  }
  return beam_load->localResultantTo(beam, node, x, side);
}

}  // namespace

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

Element& Model::element(ElementId id) {
  const auto it = elements_.find(id);
  if (it == elements_.end()) {
    throw std::out_of_range("Unknown element id");
  }
  return *it->second;
}

const Element& Model::element(ElementId id) const {
  const auto it = elements_.find(id);
  if (it == elements_.end()) {
    throw std::out_of_range("Unknown element id");
  }
  return *it->second;
}

Eigen::VectorXd Model::loadVector(double time) const {
  DofManager dm(nodes_);
  Eigen::VectorXd global_f =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(dm.size()));

  for (const auto& [node_id, node_ref] : nodes_) {
    for (std::size_t offset = 0; offset < kDofsPerFrameNode; ++offset) {
      const Dof dof = dofFromOffset(offset);
      global_f[static_cast<Eigen::Index>(dm.equation(node_id, dof))] +=
          node_ref.load(dof);
    }
  }

  const NodeResolver resolver = [this](NodeId id) -> const Node& {
    return node(id);
  };

  for (const auto& load : element_loads_) {
    const Element& target = element(load->elementId());
    const auto element_dofs = target.dofs();
    const Eigen::VectorXd fe =
        load->equivalentNodalLoadAt(target, resolver, time);
    const Eigen::Index ndof = static_cast<Eigen::Index>(element_dofs.size());

    if (fe.size() != ndof) {
      throw std::runtime_error(
          "Element load size does not match element DOF count");
    }

    for (Eigen::Index a = 0; a < ndof; ++a) {
      const auto [node_a, dof_a] =
          element_dofs[static_cast<std::size_t>(a)];
      const Eigen::Index ia =
          static_cast<Eigen::Index>(dm.equation(node_a, dof_a));
      global_f[ia] += fe[a];
    }
  }

  return global_f;
}

AssembledSystem Model::assemble() const {
  DofManager dm(nodes_);
  const Eigen::Index n = static_cast<Eigen::Index>(dm.size());

  Eigen::MatrixXd global_k = Eigen::MatrixXd::Zero(n, n);
  Eigen::MatrixXd global_m = Eigen::MatrixXd::Zero(n, n);

  const NodeResolver resolver = [this](NodeId id) -> const Node& {
    return node(id);
  };

  for (const auto& [element_id, element_ptr] : elements_) {
    (void)element_id;
    const auto element_dofs = element_ptr->dofs();
    const Eigen::MatrixXd ke = element_ptr->stiffness(resolver);
    const Eigen::MatrixXd me = element_ptr->mass(resolver);
    const Eigen::Index ndof = static_cast<Eigen::Index>(element_dofs.size());

    if (ke.rows() != ndof || ke.cols() != ndof) {
      throw std::runtime_error(
          "Element stiffness size does not match element DOF count");
    }
    if (me.rows() != ndof || me.cols() != ndof) {
      throw std::runtime_error(
          "Element mass size does not match element DOF count");
    }

    for (Eigen::Index a = 0; a < ndof; ++a) {
      const auto [node_a, dof_a] =
          element_dofs[static_cast<std::size_t>(a)];
      const Eigen::Index ia =
          static_cast<Eigen::Index>(dm.equation(node_a, dof_a));

      for (Eigen::Index b = 0; b < ndof; ++b) {
        const auto [node_b, dof_b] =
            element_dofs[static_cast<std::size_t>(b)];
        const Eigen::Index ib =
            static_cast<Eigen::Index>(dm.equation(node_b, dof_b));
        global_k(ia, ib) += ke(a, b);
        global_m(ia, ib) += me(a, b);
      }
    }
  }

  return {
      std::move(global_k),
      std::move(global_m),
      loadVector(0.0),
      std::move(dm)};
}

Eigen::VectorXd Model::elementEquivalentLoad(
    ElementId element_id,
    double time) const {
  const Element& target = element(element_id);
  const auto element_dofs = target.dofs();
  Eigen::VectorXd result =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(element_dofs.size()));

  const NodeResolver resolver = [this](NodeId id) -> const Node& {
    return node(id);
  };

  for (const auto& load : element_loads_) {
    if (load->elementId() != element_id) {
      continue;
    }
    const Eigen::VectorXd contribution =
        load->equivalentNodalLoadAt(target, resolver, time);
    if (contribution.size() != result.size()) {
      throw std::runtime_error(
          "Element load size does not match element DOF count");
    }
    result += contribution;
  }
  return result;
}

ElementResponse Model::elementResponse(
    ElementId element_id,
    const Eigen::VectorXd& global_displacement,
    double time) const {
  DofManager dm(nodes_);
  if (global_displacement.size() != static_cast<Eigen::Index>(dm.size())) {
    throw std::invalid_argument(
        "Global displacement size does not match model DOF count");
  }

  const Element& target = element(element_id);
  const auto element_dofs = target.dofs();
  Eigen::VectorXd element_u(static_cast<Eigen::Index>(element_dofs.size()));

  for (std::size_t i = 0; i < element_dofs.size(); ++i) {
    const auto [node_id, dof] = element_dofs[i];
    element_u[static_cast<Eigen::Index>(i)] =
        global_displacement[
            static_cast<Eigen::Index>(dm.equation(node_id, dof))];
  }

  const NodeResolver resolver = [this](NodeId id) -> const Node& {
    return node(id);
  };
  return target.response(
      element_u, elementEquivalentLoad(element_id, time), resolver);
}

std::vector<ElementLoadSamplingLocation> Model::elementLoadSampleLocations(
    ElementId element_id) const {
  const Element& target = element(element_id);
  const NodeResolver resolver = [this](NodeId id) -> const Node& {
    return node(id);
  };

  std::vector<ElementLoadSamplingLocation> locations;
  for (const auto& load : element_loads_) {
    if (load->elementId() != element_id) {
      continue;
    }
    auto load_locations = load->responseSampleLocations(target, resolver);
    locations.insert(
        locations.end(), load_locations.begin(), load_locations.end());
  }
  return locations;
}

BeamSectionForces Model::beamSectionForces(
    ElementId element_id,
    double x,
    const Eigen::VectorXd& global_displacement,
    double time,
    BeamSectionSide side) const {
  const Element& target = element(element_id);
  const auto* beam = dynamic_cast<const Beam3D*>(&target);
  if (beam == nullptr) {
    throw std::invalid_argument(
        "beamSectionForces requires a Beam3D element");
  }

  const NodeResolver resolver = [this](NodeId id) -> const Node& {
    return node(id);
  };
  const double l = beam->length(resolver);
  if (x < 0.0 || x > l) {
    throw std::invalid_argument(
        "Beam section coordinate must lie in [0, L]");
  }

  const ElementResponse end_response =
      elementResponse(element_id, global_displacement, time);

  const Eigen::Vector3d initial_force =
      -end_response.local_end_force.segment<3>(0);
  const Eigen::Vector3d initial_moment =
      -end_response.local_end_force.segment<3>(3);

  BeamLoadResultant3D applied;
  for (const auto& load : element_loads_) {
    if (load->elementId() != element_id) {
      continue;
    }
    const BeamLoadResultant3D contribution =
        beamLoadResultantTo(*load, *beam, resolver, x, time, side);
    applied.force += contribution.force;
    applied.moment += contribution.moment;
  }

  const Eigen::Vector3d ex = Eigen::Vector3d::UnitX();
  const Eigen::Vector3d section_force = initial_force - applied.force;
  const Eigen::Vector3d section_moment =
      initial_moment -
      x * ex.cross(initial_force) -
      applied.moment;

  return {
      x,
      section_force.x(),
      section_force.y(),
      section_force.z(),
      section_moment.x(),
      section_moment.y(),
      section_moment.z()};
}

}  // namespace fem
