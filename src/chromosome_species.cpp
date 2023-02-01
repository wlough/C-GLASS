#include "cglass/chromosome_species.hpp"

ChromosomeSpecies::ChromosomeSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::chromosome);
  printf("Initializing chromosome!!\n");
}