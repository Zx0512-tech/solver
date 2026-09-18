#pragma once
#include <cstddef>
#include <stdexcept>
namespace fem {
using NodeId = int;
using ElementId = int;
enum class Dof : std::size_t { UX=0, UY=1, UZ=2, RX=3, RY=4, RZ=5 };
constexpr std::size_t kDofsPerFrameNode = 6;
constexpr std::size_t dofOffset(Dof dof) noexcept { return static_cast<std::size_t>(dof); }
inline Dof dofFromOffset(std::size_t offset) {
  if (offset >= kDofsPerFrameNode) throw std::out_of_range("DOF offset must be in [0, 5]");
  return static_cast<Dof>(offset);
}
}  // namespace fem
