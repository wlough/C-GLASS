#include "cglass/chromosome_species.hpp"

ChromosomeSpecies::ChromosomeSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::chromosome);
}

void ChromosomeSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  if (GetNInsert() <= 0) {
    return;
  }
}

void ChromosomeSpecies::Reserve() {
  SpeciesBase::Reserve();
  members_.reserve(GetNInsert());
}

void ChromosomeSpecies::PopMember() { Species::PopMember(); }

void ChromosomeSpecies::AddMember() { Species::AddMember(); }

void ChromosomeSpecies::UpdatePositions() {
  for (auto &&chromo : members_) {
    chromo.UpdatePosition();
  }
}

void ChromosomeSpecies::GetInteractors(std::vector<Object *> &ixors) {
  for (auto &&chromo : members_) {
    chromo.GetInteractors(ixors);
  }
}