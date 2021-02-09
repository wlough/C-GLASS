#include "cglass/receptor_species.hpp"

ReceptorSpecies::ReceptorSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::receptor);
}

// Use default initialization, plus include a concentration of receptors
void ReceptorSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  concentration_ = sparams_.concentration;
}

/* To generalize, could pass vector of species and search for species
 * name. Would need to develop local random position on species' surface
 * and calculate area for each componenti, and static cast to mesh. */
void ReceptorSpecies::SetMesh(Cortex* cx) {
  if (sparams_.component.compare("cortex") == 0) {
    mesh_ = cx;
  } else {
    Logger::Error("Receptors on species not yet implemented in Receptor::Set"
                  "Component");
  }
}

void ReceptorSpecies::Reserve() {
  // Concentration overrides number of species
  if (concentration_ > 0) {
    sparams_.num = (int)round(concentration_ * mesh_->GetArea());
  }
  members_.reserve(sparams_.num);
  Logger::Debug("Reserving memory for %d members in Receptor Species",
                sparams_.num);
  mesh_->SetNBondsMax(sparams_.num - 1);
}

/* Add receptor as a member of Receptor species and as a site on
 * the cell component. */
void ReceptorSpecies::AddMember() {
  Species::AddMember();
  mesh_->AddSitePtr(&(members_.back()));
  double pos[3];
  double u[3] = {1,0,0};
  if (sparams_.component.compare("cortex") == 0) {
   rng_.RandomBoundaryCoordinate(space_, pos);
  } else {
    Logger::Error("Receptors on species not yet implemented in Receptor::Set"
                  "Component");
  }
  members_.back().InsertAt(pos, u);
}

