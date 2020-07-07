/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : optical_trap_species
 * @created     : Tuesday Jul 07, 2020 09:43:46 MDT
 */

#ifndef _CGLASS_OPTICAL_TRAP_SPECIES_HPP_

#define _CGLASS_OPTICAL_TRAP_SPECIES_HPP_

#include "optical_trap.hpp"
#include "species.hpp"

class OpticalTrapSpecies
    : public Species<OpticalTrap, species_id::optical_trap> {
public:
  OpticalTrapSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);

private:
  /* private data */
};

#endif

