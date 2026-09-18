#pragma once

#include <vector>

namespace fem {

struct ScalarHistory {
  std::vector<double> time;
  std::vector<double> value;
};

}  // namespace fem
