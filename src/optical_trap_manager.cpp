/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : optical_trap_manager
 * @created     : Tuesday Jul 07, 2020 11:19:30 MDT
 */

#include "cglass/optical_trap_manager.hpp"

void OpticalTrapManager::Init(system_parameters *params, space_struct *space) {
  update_ = false;
  params_ = params;
  space_ = space;
}

void OpticalTrapManager::InitSpecies(sid_label &slab, ParamsParser &parser,
                                     unsigned long seed) {
  if (otrap_species_.size() == 0) {
    otrap_species_.reserve(parser.GetNOpticalTrapSpecies());
  }
  otrap_species_.push_back(new OpticalTrapSpecies(seed));
  otrap_species_.back()->Init(slab.second, parser);
  // Delete all optical trap species if there are no optical traps.
  if (otrap_species_.back()->GetNInsert() <= 0) {
    delete otrap_species_.back();
    otrap_species_.pop_back();
  }
}

void OpticalTrapManager::InsertOpticalTraps(
    std::vector<SpeciesBase *> *species) {
  for (auto it = otrap_species_.begin(); it != otrap_species_.end(); ++it) {
    (*it)->InsertOpticalTraps(species);
  }
}

void OpticalTrapManager::UpdateOpticalTraps() {
  for (auto ot_spec : otrap_species_) {
    ot_spec->ApplyOpticalTrapForces();
  }
}

void OpticalTrapManager::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto &ot_spec : otrap_species_) {
    ot_spec->Draw(graph_array);
  }
}
