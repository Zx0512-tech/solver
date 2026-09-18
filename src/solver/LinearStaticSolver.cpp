#include "fem/solver/LinearStaticSolver.hpp"
#include <Eigen/LU>
#include <stdexcept>
#include <vector>
namespace fem {
StaticResult LinearStaticSolver::solve(const Model& model) const {
  AssembledSystem sys=model.assemble();
  const Eigen::Index n=sys.load.size();
  std::vector<Eigen::Index> free;
  free.reserve(static_cast<std::size_t>(n));
  for(NodeId nid:sys.dofs.nodeOrder()){
    const Node& nd=model.node(nid);
    for(std::size_t o=0;o<kDofsPerFrameNode;++o){
      const Dof d=dofFromOffset(o);
      if(!nd.isFixed(d)) free.push_back(static_cast<Eigen::Index>(sys.dofs.equation(nid,d)));
    }
  }
  if(free.empty()) throw std::runtime_error("Model has no free degrees of freedom");
  const Eigen::Index nf=static_cast<Eigen::Index>(free.size());
  Eigen::MatrixXd kff(nf,nf); Eigen::VectorXd ff(nf);
  for(Eigen::Index i=0;i<nf;++i){
    ff[i]=sys.load[free[static_cast<std::size_t>(i)]];
    for(Eigen::Index j=0;j<nf;++j) kff(i,j)=sys.stiffness.coeff(free[static_cast<std::size_t>(i)],free[static_cast<std::size_t>(j)]);
  }
  Eigen::FullPivLU<Eigen::MatrixXd> lu(kff);
  if(!lu.isInvertible()) throw std::runtime_error("Reduced stiffness matrix is singular; check constraints, connectivity and section properties");
  const Eigen::VectorXd uf=lu.solve(ff);
  Eigen::VectorXd u=Eigen::VectorXd::Zero(n);
  for(Eigen::Index i=0;i<nf;++i) u[free[static_cast<std::size_t>(i)]]=uf[i];
  const Eigen::VectorXd reaction=sys.stiffness*u-sys.load;
  return {std::move(u),reaction,std::move(sys.dofs)};
}
}  // namespace fem
