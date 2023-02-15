#ifndef _CGLASS_CHROMOSOME_SPECIES_H_
#define _CGLASS_CHROMOSOME_SPECIES_H_

#include "KMC/kmc.hpp"
#include "chromosome.hpp"
#include "species.hpp"

class ChromosomeSpecies : public Species<Chromosome, species_id::chromosome> {
private:
  bool *update_;
  std::string checkpoint_file_;
  std::vector<Object *> *objs_;

  Tracker *tracker_ = nullptr;
  LookupTable lut_;

public:
private:
public:
  // ChromosomeSpecies(){};
  ChromosomeSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void PopMember();
  void AddMember();

  void Reserve();

  void UpdatePositions();

  void GetInteractors(std::vector<Object *> &ixors);

  void GetLastInteractors(std::vector<Object *> &ixors) {
    ixors.push_back(&members_.back());
  }
};

#endif