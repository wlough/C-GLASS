#include "cglass/centrosome_species.hpp"

CentrosomeSpecies::CentrosomeSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::centrosome);
  printf("NEW centrosome!\n");
}

void CentrosomeSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  if (GetNInsert() <= 0) {
    return;
  }
  printf("Initializing centrosome\n");
}

void CentrosomeSpecies::Reserve() {
  SpeciesBase::Reserve();
  members_.reserve(GetNInsert());
}

void CentrosomeSpecies::PopMember() { Species::PopMember(); }
void CentrosomeSpecies::AddMember() { Species::AddMember(); }

void CentrosomeSpecies::UpdatePositions() {
  for (auto &&centro : members_) {
    centro.UpdatePosition();
  }
}

void CentrosomeSpecies::GetInteractors(std::vector<Object *> &ixors) {
  for (auto &&centro : members_) {
    centro.GetInteractors(ixors);
  }
}