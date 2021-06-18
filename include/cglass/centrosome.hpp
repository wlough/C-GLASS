#ifndef _CGLASS_CENTROSOME_H_
#define _CGLASS_CENTROSOME_H_

#include "species.hpp"
#include "filament.hpp"

class Centrosome : public Object {
 protected:
  bool alignment_potential_, fixed_spacing_;
  int n_filaments_, n_filaments_min_, n_filaments_max_;
  double k_spring_, k_align_, spring_length_, anchor_distance_, gamma_trans_,
      gamma_rot_, diffusion_;
  std::vector<Filament> filaments_;
  std::vector<Anchor> anchors_;
  void ApplyForcesTorques();
  void ApplyBoundaryForces();
  void InsertCentrosome();
  void GenerateAnchorSites();
  void SetDiffusion();
  void Translate();
  void Rotate();
  void Integrate();
  void RandomizeAnchorPosition(int i_fil);

 public:
  Centrosome();
  void Init();
  void UpdatePosition() {}
  void UpdatePosition(bool midstep);
  virtual std::vector<Object *> GetInteractors();
  virtual int GetCount();
  virtual void Draw(std::vector<graph_struct *> *graph_array);
  virtual void ZeroForce();
};

class CentrosomeSpecies : public Species<Centrosome> {
 protected:
  bool midstep_;

 public:
  CentrosomeSpecies() : Species() {
    SetSID(species_id::centrosome);
    midstep_ = true;
  }
  void Init(system_parameters *params, SpaceBase *space, long seed) {
    Species::Init(params, space, seed);
    sparams_ = &(params_->centrosome);
    name_ = sparams_->name;
  }
  void UpdatePositions() {
    for (auto it = members_.begin(); it != members_.end(); ++it) {
      it->UpdatePosition(midstep_);
    }
    midstep_ = !midstep_;
  }
};

#endif  // _CGLASS_CENTROSOME_H_
