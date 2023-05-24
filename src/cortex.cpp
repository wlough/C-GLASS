#include "cglass/cortex.hpp"
#include "cglass/receptor.hpp"

Cortex::Cortex(unsigned long seed) : PointCover(seed) {
  type_ = obj_type::cortex;
}

// Add sites as interactors; not currently used
void Cortex::UpdateInteractors() {
  interactors_.clear();
  for (auto it = sphere_ptrs_.begin(); it != sphere_ptrs_.end(); ++it) {
    interactors_.push_back(*it);
  }
}
