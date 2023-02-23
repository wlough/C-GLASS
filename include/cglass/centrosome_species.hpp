#ifndef _CGLASS_CENTROSOME_SPECIES_H_
#define _CGLASS_CENTROSOME_SPECIES_H_

#include "centrosome.hpp"
#include "species.hpp"

class CentrosomeSpecies : public Species<Centrosome, species_id::centrosome> {
protected:
  bool midstep_;

public:
  CentrosomeSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void PopMember();
  void AddMember();

  void Reserve();

  void UpdatePositions();

  void GetInteractors(std::vector<Object *> &ixors);

  void GetLastInteractors(std::vector<Object *> &ixors) {}
};

#endif
