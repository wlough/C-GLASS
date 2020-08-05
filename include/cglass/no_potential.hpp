#ifndef _CGLASS_NO_POTENTIAL_H_
#define _CGLASS_NO_POTENTIAL_H_

#include "auxiliary.hpp"
#include "interaction.hpp"
#include "potential_base.hpp"

class NoPotential : public PotentialBase {
public:
  NoPotential() {}
  void CalcPotential(Interaction &ix) {
    for (int i = 0; i < n_dim_; ++i) {
      ix.force[i] = 0;
    }
    std::fill(ix.force, ix.force + n_dim_, 0);
    std::fill(ix.stress, ix.stress + n_dim_*n_dim_, 0);
    ix.pote = 0;
  }

  void InitPotentialParams(system_parameters *params) {
    rcut_ = 0;
    rcut2_ = 0; // should never interact
  }
};

#endif
