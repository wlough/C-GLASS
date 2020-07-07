/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : optical_trap_species
 * @created     : Tuesday Jul 07, 2020 09:41:45 MDT
 */

#include "cglass/optical_trap_species.hpp"

OpticalTrapSpecies::OpticalTrapSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::optical_trap);
}

void OpticalTrapSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
}
