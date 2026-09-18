#include "fem/core/Material.hpp"
#include "fem/core/Model.hpp"
#include "fem/core/Section.hpp"
#include "fem/elements/Beam3D.hpp"
#include "fem/solver/LinearStaticSolver.hpp"
#include <cstddef>
#include <iomanip>
#include <iostream>
int main(){
  fem::Model model;
  auto& root=model.addNode(1,{0.0,0.0,0.0});
  auto& tip=model.addNode(2,{2.0,0.0,0.0});
  for(std::size_t i=0;i<fem::kDofsPerFrameNode;++i) root.fix(fem::dofFromOffset(i));
  const fem::LinearElasticMaterial steel(210.0e9,0.30,7850.0);
  const fem::BeamSection section(8.0e-3,6.0e-6,8.0e-6,1.0e-5);
  model.addElement<fem::Beam3D>(1,1,2,steel,section);
  tip.addLoad(fem::Dof::UY,-10000.0);
  const auto result=fem::LinearStaticSolver{}.solve(model);
  std::cout<<std::scientific<<std::setprecision(6)
           <<"tip Uy = "<<result.displacementAt(2,fem::Dof::UY)<<" m\n"
           <<"root Fy reaction = "<<result.reactionAt(1,fem::Dof::UY)<<" N\n";
}
