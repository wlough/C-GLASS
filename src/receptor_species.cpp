#include "cglass/receptor_species.hpp"

ReceptorSpecies::ReceptorSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::receptor);
}

// Use default initialization, plus include a concentration of receptors
void ReceptorSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  concentration_ = sparams_.concentration;
  component_ = sparams_.component;
}

/* To generalize, could pass vector of species and search for species
 * name. Would need to develop local random position on species' surface
 * and calculate area for each component, and static cast to mesh. */
void ReceptorSpecies::SetPC(Cortex* cx, std::vector<SpeciesBase *> &species) {
  if (component_.compare("cortex") == 0) {
    pc_ = cx;
  } else {
    for (auto sp = species.begin(); sp != species.end(); ++sp) {
      if ((*sp)->GetSpeciesName().compare(component_) == 0) {
        pc_ = (*sp)->GetPC();
        pc_species_ = *sp;
      }
    }
    if (!pc_) {
      Logger::Error("Receptor component name is not a species name.");
    }
  }
}

void ReceptorSpecies::Reserve() {
  // Concentration overrides number of species
  if (concentration_ > 0) {
    sparams_.num = (int)round(concentration_ * pc_->GetArea());
  }
  members_.reserve(sparams_.num);
  Logger::Debug("Reserving memory for %d members in Receptor Species",
                sparams_.num);
}

/* Add receptor as a member of Receptor species and as a site on
 * the cell component. */
void ReceptorSpecies::AddMember() {
  Species::AddMember();
  pc_->AddSpherePtr(&(members_.back()));
  double pos[3] =  {0,0,0};
  double u[3] = {1,0,0};
  if (component_.compare("cortex") == 0) {
    rng_.RandomBoundaryCoordinate(space_, pos);
  } else {
    s_ = pc_species_->GetSpecLength() * (rng_.RandomUniform()-0.5);
    i_ = rng_.RandomInt(pc_species_->GetNMembers());
    pc_species_->CalcPCPosition(i_, s_, pos);
  }
  members_.back().InsertAt(pos, u);
  members_.back().SetCompPtr(pc_);
}

