#include "cglass/chromosome.hpp"
#include <iostream>

Chromosome::Chromosome(unsigned long seed) : Object(seed) {
  SetSID(species_id::chromosome);
}