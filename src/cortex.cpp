#include "cglass/cortex.hpp"

Cortex::Cortex(unsigned long seed) : Mesh(seed) {
  type_ = obj_type::cortex;
  n_bonds_max_ = 0;
}

void Cortex::InitSiteAt(double *new_pos, double d) {
  n_bonds_max_++;
  InitSiteAt(new_pos, d);
}

void Cortex::UpdateInteractors() {
  interactors_.clear();
  for (auto it = site_ptrs_.begin(); it != site_ptrs_.end(); ++it) {
    interactors_.push_back(*it);
  }
}

void Cortex::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto site_it = sites_.begin(); site_it != sites_.end(); ++site_it) {
    site_it->Draw(graph_array);
  }
}
