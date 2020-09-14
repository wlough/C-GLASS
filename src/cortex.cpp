#include "cglass/cortex.hpp"

Cortex::Cortex(unsigned long seed) : Mesh(seed) {
  total_sites_ = 0.0;
}

void Cortex::Init(system_parameters *params) {
  params_ = params;
  SetParameters();
  if ((site_concentration_ > 0) && (site_diameter_ > 0)) {
    total_sites_ = (int)round(site_concentration_ * space_->BoundaryArea());
  }
  n_bonds_max_ = total_sites_ - 1;
  AddSites();
}

void Cortex::SetParameters() {
  site_concentration_ = params_->cortex_site_concentration;
  site_diameter_ = params_->cortex_site_diameter;
}

void Cortex::AddSites() {
  double pos[3];
  for (int i = 0; i < total_sites_; ++i) {
    rng_.RandomBoundaryCoordinate(space_, pos);
    InitSiteAt(pos, site_diameter_);
  }
}

void Cortex::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto site_it = sites_.begin(); site_it != sites_.end(); ++site_it) {
    site_it->Draw(graph_array);
  }
}
