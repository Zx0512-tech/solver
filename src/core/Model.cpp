#include "fem/core/Model.hpp"
#include <stdexcept>
namespace fem {
Node& Model::addNode(NodeId id,const Eigen::Vector3d& coordinates) {
  auto [it,inserted]=nodes_.emplace(id,Node{id,coordinates});
  if (!inserted) throw std::invalid_argument("Duplicate node id");
  return it->second;
}
Node& Model::node(NodeId id) {
  const auto it=nodes_.find(id);
  if (it==nodes_.end()) throw std::out_of_range("Unknown node id");
  return it->second;
}
const Node& Model::node(NodeId id) const {
  const auto it=nodes_.find(id);
  if (it==nodes_.end()) throw std::out_of_range("Unknown node id");
  return it->second;
}
AssembledSystem Model::assemble() const {
  DofManager dm(nodes_);
  const Eigen::Index n=static_cast<Eigen::Index>(dm.size());
  Eigen::MatrixXd k=Eigen::MatrixXd::Zero(n,n);
  Eigen::VectorXd f=Eigen::VectorXd::Zero(n);
  for (const auto& [node_id,nref]:nodes_) {
    for (std::size_t offset=0;offset<kDofsPerFrameNode;++offset) {
      const Dof dof=dofFromOffset(offset);
      f[static_cast<Eigen::Index>(dm.equation(node_id,dof))]+=nref.load(dof);
    }
  }
  const NodeResolver resolver=[this](NodeId id)->const Node& { return node(id); };
  for (const auto& [element_id,element]:elements_) {
    (void)element_id;
    const auto edofs=element->dofs();
    const Eigen::MatrixXd ke=element->stiffness(resolver);
    const Eigen::Index ndof=static_cast<Eigen::Index>(edofs.size());
    if (ke.rows()!=ndof || ke.cols()!=ndof) throw std::runtime_error("Element stiffness size does not match element DOF count");
    for (Eigen::Index a=0;a<ndof;++a) {
      const auto [na,da]=edofs[static_cast<std::size_t>(a)];
      const Eigen::Index ia=static_cast<Eigen::Index>(dm.equation(na,da));
      for (Eigen::Index b=0;b<ndof;++b) {
        const auto [nb,db]=edofs[static_cast<std::size_t>(b)];
        const Eigen::Index ib=static_cast<Eigen::Index>(dm.equation(nb,db));
        k(ia,ib)+=ke(a,b);
      }
    }
  }
  return {std::move(k),std::move(f),std::move(dm)};
}
}  // namespace fem
