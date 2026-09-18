#include "fem/elements/Beam3D.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <stdexcept>
namespace fem {
namespace { constexpr double kTol=1.0e-12; }
Beam3D::Beam3D(ElementId id,NodeId ni,NodeId nj,LinearElasticMaterial material,BeamSection section,const Eigen::Vector3d& yhint)
  :Element(id),node_i_(ni),node_j_(nj),material_(material),section_(section),local_y_hint_(yhint) {
  if (node_i_==node_j_) throw std::invalid_argument("Beam3D requires two distinct node ids");
  if (local_y_hint_.norm()<=kTol) throw std::invalid_argument("Beam3D local y-axis hint cannot be zero");
}
std::vector<ElementDof> Beam3D::dofs() const {
  std::vector<ElementDof> out; out.reserve(12);
  for (NodeId n:{node_i_,node_j_}) for (std::size_t i=0;i<kDofsPerFrameNode;++i) out.emplace_back(n,dofFromOffset(i));
  return out;
}
Eigen::MatrixXd Beam3D::stiffness(const NodeResolver& node) const { return globalStiffness(node); }
double Beam3D::length(const NodeResolver& node) const {
  const double l=(node(node_j_).coordinates()-node(node_i_).coordinates()).norm();
  if (l<=kTol) throw std::runtime_error("Beam3D has zero or near-zero length");
  return l;
}
Beam3D::Matrix12d Beam3D::localStiffness(double l) const {
  if (l<=kTol) throw std::invalid_argument("Beam3D length must be positive");
  const double E=material_.youngs_modulus,G=material_.shearModulus();
  const double A=section_.area,Iy=section_.iy,Iz=section_.iz,J=section_.torsion_constant;
  const double L2=l*l,L3=L2*l;
  Matrix12d k=Matrix12d::Zero();
  const double ea=E*A/l, gj=G*J/l;
  k(0,0)=ea;k(0,6)=-ea;k(6,0)=-ea;k(6,6)=ea;
  k(3,3)=gj;k(3,9)=-gj;k(9,3)=-gj;k(9,9)=gj;
  const double a=12*E*Iz/L3,b=6*E*Iz/L2,c=4*E*Iz/l,d=2*E*Iz/l;
  const int bz[4]={1,5,7,11};
  const double vz[4][4]={{a,b,-a,b},{b,c,-b,d},{-a,-b,a,-b},{b,d,-b,c}};
  for(int r=0;r<4;++r) for(int q=0;q<4;++q) k(bz[r],bz[q])=vz[r][q];
  const double e=12*E*Iy/L3,f=6*E*Iy/L2,g=4*E*Iy/l,h=2*E*Iy/l;
  const int by[4]={2,4,8,10};
  const double vy[4][4]={{e,-f,-e,-f},{-f,g,f,h},{-e,f,e,f},{-f,h,f,g}};
  for(int r=0;r<4;++r) for(int q=0;q<4;++q) k(by[r],by[q])=vy[r][q];
  return k;
}
Eigen::Matrix3d Beam3D::rotationToLocal(const NodeResolver& node) const {
  const Eigen::Vector3d x=(node(node_j_).coordinates()-node(node_i_).coordinates()).normalized();
  Eigen::Vector3d yp=local_y_hint_-local_y_hint_.dot(x)*x;
  if (yp.norm()<=kTol) {
    const Eigen::Vector3d axes[3]={Eigen::Vector3d::UnitX(),Eigen::Vector3d::UnitY(),Eigen::Vector3d::UnitZ()};
    int best=0; double min=std::abs(x.dot(axes[0]));
    for(int i=1;i<3;++i){ const double v=std::abs(x.dot(axes[i])); if(v<min){min=v;best=i;} }
    yp=axes[best]-axes[best].dot(x)*x;
  }
  const Eigen::Vector3d y=yp.normalized();
  const Eigen::Vector3d z=x.cross(y).normalized();
  const Eigen::Vector3d yc=z.cross(x).normalized();
  Eigen::Matrix3d r; r.row(0)=x.transpose();r.row(1)=yc.transpose();r.row(2)=z.transpose(); return r;
}
Beam3D::Matrix12d Beam3D::transformation(const NodeResolver& node) const {
  const Eigen::Matrix3d r=rotationToLocal(node);
  Matrix12d t=Matrix12d::Zero();
  for(int b=0;b<4;++b) t.block<3,3>(3*b,3*b)=r;
  return t;
}
Beam3D::Matrix12d Beam3D::globalStiffness(const NodeResolver& node) const {
  const Matrix12d t=transformation(node);
  return t.transpose()*localStiffness(length(node))*t;
}
}  // namespace fem
