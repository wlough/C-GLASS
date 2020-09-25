#ifndef _CGLASS_RECEPTOR_SPECIES_H_
#define _CGLASS_RECEPTOR_SPECIES_H_

#include "receptor.hpp"
#include "species.hpp"

class ReceptorSpecies: public Species<Receptor, species_id::receptor> {
private:
  double concentration_;
//  int n_insert_;
  Mesh* component_;

public:
  ReceptorSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void SetComponent(Cortex *cortex);
  void Reserve();
  void AddMember();
//  const int const GetNInsert() {return n_insert_;}
};

#endif
