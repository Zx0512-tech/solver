#pragma once

#include <Eigen/Core>

namespace fem {

// Physical resultant of beam element loads accumulated up to a section.
// force is expressed in beam-local axes. moment is the resultant moment of
// those loads about the queried section, also in beam-local axes.
struct BeamLoadResultant3D {
  Eigen::Vector3d force{Eigen::Vector3d::Zero()};
  Eigen::Vector3d moment{Eigen::Vector3d::Zero()};
};

}  // namespace fem
