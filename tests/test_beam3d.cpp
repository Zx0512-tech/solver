#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/solver/LinearStaticSolver.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
namespace {
void near(double actual,double expected,double tol,const std::string& msg){
  const double scale=std::max(1.0,std::abs(expected));
  if(std::abs(actual-expected)>tol*scale){std::cerr<<"FAIL "<<msg<<" actual="<<actual<<" expected="<<expected<<"\n";std::exit(EXIT_FAILURE);}
}
void fixAll(fem::Node& n){for(std::size_t i=0;i<fem::kDofsPerFrameNode;++i)n.fix(fem::dofFromOffset(i));}
void stiffnessSymmetry(){
  fem::Beam3D beam(1,1,2,{210e9,0.3},{0.01,2e-6,3e-6,4e-6});
  const auto k=beam.localStiffness(2.0);
  near((k-k.transpose()).cwiseAbs().maxCoeff(),0.0,1e-12,"stiffness symmetry");
}
void axial(){
  fem::Model m; auto& r=m.addNode(1,{0,0,0}); auto& t=m.addNode(2,{2,0,0}); fixAll(r);
  const double E=210e9,A=.01,L=2,P=50000;
  m.addElement<fem::Beam3D>(1,1,2,fem::LinearElasticMaterial(E,.3),fem::BeamSection(A,2e-6,3e-6,4e-6));
  t.addLoad(fem::Dof::UX,P); const auto s=fem::LinearStaticSolver{}.solve(m);
  near(s.displacementAt(2,fem::Dof::UX),P*L/(E*A),1e-10,"axial displacement");
  near(s.reactionAt(1,fem::Dof::UX),-P,1e-10,"axial reaction");
}
void bendingZ(){
  fem::Model m; auto& r=m.addNode(1,{0,0,0}); auto& t=m.addNode(2,{3,0,0}); fixAll(r);
  const double E=200e9,Iz=8e-6,L=3,P=-12000;
  m.addElement<fem::Beam3D>(1,1,2,fem::LinearElasticMaterial(E,.3),fem::BeamSection(.01,6e-6,Iz,1e-5));
  t.addLoad(fem::Dof::UY,P); const auto s=fem::LinearStaticSolver{}.solve(m);
  near(s.displacementAt(2,fem::Dof::UY),P*L*L*L/(3*E*Iz),1e-10,"local-z bending tip deflection");
  near(s.displacementAt(2,fem::Dof::RZ),P*L*L/(2*E*Iz),1e-10,"local-z bending tip rotation");
  near(s.reactionAt(1,fem::Dof::UY),-P,1e-10,"local-z bending shear reaction");
  near(s.reactionAt(1,fem::Dof::RZ),-P*L,1e-10,"local-z bending moment reaction");
}
void bendingY(){
  fem::Model m; auto& r=m.addNode(1,{0,0,0}); auto& t=m.addNode(2,{3,0,0}); fixAll(r);
  const double E=200e9,Iy=6e-6,L=3,P=12000;
  m.addElement<fem::Beam3D>(1,1,2,fem::LinearElasticMaterial(E,.3),fem::BeamSection(.01,Iy,8e-6,1e-5));
  t.addLoad(fem::Dof::UZ,P); const auto s=fem::LinearStaticSolver{}.solve(m);
  near(s.displacementAt(2,fem::Dof::UZ),P*L*L*L/(3*E*Iy),1e-10,"local-y bending tip deflection");
  near(s.displacementAt(2,fem::Dof::RY),-P*L*L/(2*E*Iy),1e-10,"local-y bending tip rotation");
  near(s.reactionAt(1,fem::Dof::UZ),-P,1e-10,"local-y bending shear reaction");
  near(s.reactionAt(1,fem::Dof::RY),P*L,1e-10,"local-y bending moment reaction");
}
void torsion(){
  fem::Model m; auto& r=m.addNode(1,{0,0,0}); auto& t=m.addNode(2,{2.5,0,0}); fixAll(r);
  const double E=210e9,nu=.3,J=1.2e-5,L=2.5,T=9000;
  const double G=E/(2*(1+nu));
  m.addElement<fem::Beam3D>(1,1,2,fem::LinearElasticMaterial(E,nu),fem::BeamSection(.01,6e-6,8e-6,J));
  t.addLoad(fem::Dof::RX,T); const auto s=fem::LinearStaticSolver{}.solve(m);
  near(s.displacementAt(2,fem::Dof::RX),T*L/(G*J),1e-10,"torsional rotation");
  near(s.reactionAt(1,fem::Dof::RX),-T,1e-10,"torsional reaction");
}
void rotated(){
  fem::Model m; auto& r=m.addNode(1,{0,0,0}); auto& t=m.addNode(2,{1,1,0}); fixAll(r);
  const double E=70e9,A=.02,P=20000,L=std::sqrt(2.0);
  const Eigen::Vector3d axis=Eigen::Vector3d(1,1,0).normalized();
  m.addElement<fem::Beam3D>(1,1,2,fem::LinearElasticMaterial(E,.33),fem::BeamSection(A,2e-5,3e-5,4e-5),Eigen::Vector3d::UnitZ());
  t.addLoad(fem::Dof::UX,P*axis.x()); t.addLoad(fem::Dof::UY,P*axis.y());
  const auto s=fem::LinearStaticSolver{}.solve(m);
  const Eigen::Vector3d u(s.displacementAt(2,fem::Dof::UX),s.displacementAt(2,fem::Dof::UY),s.displacementAt(2,fem::Dof::UZ));
  near(u.dot(axis),P*L/(E*A),1e-10,"rotated axial deformation");
}
}
int main(){stiffnessSymmetry();axial();bendingZ();bendingY();torsion();rotated();std::cout<<"All Beam3D tests passed.\n";return EXIT_SUCCESS;}
