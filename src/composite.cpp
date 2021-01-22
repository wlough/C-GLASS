#include "cglass/composite.hpp"

int Composite::_next_comp_id_ = 0;
std::mutex Composite::_comp_mtx_;

Composite::Composite(unsigned long seed) : Object(seed) {
  InitCompID();
  is_comp_ = true;
  comp_type_ = comp_type::generic;
}
void Composite::InitCompID() {
  std::lock_guard<std::mutex> lk(_comp_mtx_);
  SetMeshID(++_next_comp_id_);
}