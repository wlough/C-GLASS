#ifndef _CGLASS_CHROMOSOME_MANAGER_H_
#define _CGLASS_CHROMOSOME_MANAGER_H_

#include "output_manager.hpp"
#include "rng.hpp"
#include "species_factory.hpp"
#include <numeric>

class ChromosomeOutputManager : public OutputManagerBase<ChromosomeSpecies> {
  void WriteThermo() {}
  void ReadThermo() {}
  void InitThermo(std::string fname) {}
};

class ChromosomeManager {
private:
  system_parameters *params_;
  ChromosomeOutputManager output_mgr_;
  double obj_size_;
  double rcutoff_ = 0; // Cutoff for binding any crosslink and bond
  bool update_;

  std::vector<ChromosomeSpecies *> chromo_species_;
  std::vector<Object *> *objs_;
  std::map<std::string, ChromosomeSpecies *> species_map_;
  SpaceBase *space_;
  Tracker *tracker_ = nullptr;
  RNG *rng_;

public:
private:
public:
};

#endif