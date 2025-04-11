#pragma once

#include "membrane.hpp"
#include "species.hpp"

class MembraneSpecies : public Species<Membrane, species_id::membrane> {
public:
  MembraneSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
};
