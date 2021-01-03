#include "cglass/cortex.hpp"

Cortex::Cortex(unsigned long seed) : Mesh(seed) {
  type_ = obj_type::cortex;
  // No bonds included in Cortex yet
  n_bonds_max_ = 0;
}

// Add a site with diameter d at a given position
void Cortex::InitSiteAt(double *new_pos, double d) {
  n_bonds_max_++;
  Mesh::InitSiteAt(new_pos, d);
}

// Add sites as interactors; not currently used
void Cortex::UpdateInteractors() {
  interactors_.clear();
  for (auto it = site_ptrs_.begin(); it != site_ptrs_.end(); ++it) {
    interactors_.push_back(*it);
  }
}

// Draw all sites on cortex; not currently used
void Cortex::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto site_it = sites_.begin(); site_it != sites_.end(); ++site_it) {
    site_it->Draw(graph_array);
  }
}


