#pragma once

#include "fem/core/DofManager.hpp"
#include "fem/core/Element.hpp"
#include "fem/core/ElementLoad.hpp"
#include "fem/core/ElementResponse.hpp"
#include "fem/core/Node.hpp"
#include "fem/core/Types.hpp"
#include "fem/loads/TimeDependentElementLoad.hpp"
#include "fem/response/BeamSectionForces.hpp"

#include <Eigen/Core>

#include <functional>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fem {

struct AssembledSystem {
  Eigen::MatrixXd stiffness;
  Eigen::MatrixXd mass;
  Eigen::VectorXd load;
  DofManager dofs;
};

class Model {
 public:
  Node& addNode(NodeId id, const Eigen::Vector3d& coordinates);

  template <typename ElementType, typename... Args>
  ElementType& addElement(Args&&... args) {
    auto element = std::make_unique<ElementType>(std::forward<Args>(args)...);
    const ElementId id = element->id();
    if (elements_.count(id) != 0U) {
      throw std::invalid_argument("Duplicate element id");
    }
    ElementType& ref = *element;
    elements_.emplace(id, std::move(element));
    return ref;
  }

  template <typename LoadType, typename... Args>
  LoadType& addElementLoad(Args&&... args) {
    auto load = std::make_unique<LoadType>(std::forward<Args>(args)...);
    (void)element(load->elementId());
    LoadType& ref = *load;
    element_loads_.push_back(std::move(load));
    return ref;
  }

  template <typename SpatialLoadType, typename... Args>
  TimeDependentElementLoad& addTimeDependentElementLoad(
      std::function<double(double)> scale,
      Args&&... args) {
    auto spatial =
        std::make_unique<SpatialLoadType>(std::forward<Args>(args)...);
    (void)element(spatial->elementId());
    auto load = std::make_unique<TimeDependentElementLoad>(
        std::move(spatial), std::move(scale));
    TimeDependentElementLoad& ref = *load;
    element_loads_.push_back(std::move(load));
    return ref;
  }

  Node& node(NodeId id);
  const Node& node(NodeId id) const;

  Element& element(ElementId id);
  const Element& element(ElementId id) const;

  AssembledSystem assemble() const;
  Eigen::VectorXd loadVector(double time = 0.0) const;

  Eigen::VectorXd elementEquivalentLoad(
      ElementId element_id,
      double time = 0.0) const;

  ElementResponse elementResponse(
      ElementId element_id,
      const Eigen::VectorXd& global_displacement,
      double time = 0.0) const;

  BeamSectionForces beamSectionForces(
      ElementId element_id,
      double x,
      const Eigen::VectorXd& global_displacement,
      double time = 0.0,
      BeamSectionSide side = BeamSectionSide::Right) const;

  std::vector<ElementLoadSamplingLocation> elementLoadSampleLocations(
      ElementId element_id) const;

  const std::unordered_map<NodeId, Node>& nodes() const noexcept { return nodes_; }

 private:
  std::unordered_map<NodeId, Node> nodes_;
  std::unordered_map<ElementId, std::unique_ptr<Element>> elements_;
  std::vector<std::unique_ptr<ElementLoad>> element_loads_;
};

}  // namespace fem
