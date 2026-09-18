#pragma once

#include "fem/core/ElementResponse.hpp"
#include "fem/core/Node.hpp"
#include "fem/core/Types.hpp"

#include <Eigen/Core>

#include <functional>
#include <utility>
#include <vector>

namespace fem {

using NodeResolver = std::function<const Node&(NodeId)>;
using ElementDof = std::pair<NodeId, Dof>;

class Element {
 public:
  explicit Element(ElementId id) : id_(id) {}
  virtual ~Element() = default;

  ElementId id() const noexcept { return id_; }

  virtual std::vector<ElementDof> dofs() const = 0;
  virtual Eigen::MatrixXd stiffness(const NodeResolver& node) const = 0;
  virtual Eigen::MatrixXd mass(const NodeResolver& node) const = 0;

  // element_global_displacement and element_equivalent_load_global follow dofs().
  // The returned forces are resisting/end forces, i.e. elastic response minus
  // the equivalent nodal load carried by the element itself.
  virtual ElementResponse response(
      const Eigen::VectorXd& element_global_displacement,
      const Eigen::VectorXd& element_equivalent_load_global,
      const NodeResolver& node) const = 0;

 private:
  ElementId id_;
};

}  // namespace fem
