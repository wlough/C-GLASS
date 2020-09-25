#include "cglass/receptor_species.hpp"

ReceptorSpecies::ReceptorSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::receptor);
}
void ReceptorSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  concentration_ = sparams_.concentration;
}

/* To generalize, could pass vector of species and search for species
 * name. Would need to develop local random position on species' surface
 * and calculate area for each componenti, and static cast to mesh. */
void ReceptorSpecies::SetComponent(Cortex* cx) {
  if (sparams_.component.compare("cortex") == 0) {
    component_ = cx;
  } else {
    Logger::Error("Receptors on species not yet implemented in Receptor::Set"
                  "Component");
  }
}

void ReceptorSpecies::Reserve() {
  if (concentration_ > 0) {
    sparams_.num = (int)round(concentration_ * component_->GetArea());
  }
  members_.reserve(sparams_.num);
  Logger::Debug("Reserving memory for %d members in Receptor Species",
                sparams_.num);
  component_->SetNBondsMax(sparams_.num - 1);
}

// Add receptor as a member of Receptor species and as a site on
// the cell component.
void ReceptorSpecies::AddMember() {
  Species::AddMember();
  component_->AddSitePtr(&(members_.back()));
  double pos[3];
  if (sparams_.component.compare("cortex") == 0) {
    rng_.RandomBoundaryCoordinate(space_, pos);
  } else {
    Logger::Error("Receptors on species not yet implemented in Receptor::Set"
                  "Component");
  }
  members_.back().SetPosition(pos);
}

