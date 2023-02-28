#ifndef _CGLASS_CENTROSOME_H_
#define _CGLASS_CENTROSOME_H_

#include "anchor.hpp"
#include "filament.hpp"
#include <armadillo>

class Centrosome : public Object {
  // Direction vectors that define plane of SPB
  // (r = position_ and u = orentation_)

protected:
  bool zero_temperature_ = false;

  bool alignment_potential_, fixed_spacing_;
  int n_filaments_, n_filaments_min_, n_filaments_max_;
  double k_spring_, k_align_, spring_length_, anchor_distance_, gamma_trans_,
      gamma_rot_, diffusion_;

  centrosome_parameters *sparams_;

  double gamma_perp_ = 0;
  double noise_tr_ = 0.0;
  double noise_rot_ = 0.0;
  double diffusion_rot_ = 0.0;
  double body_frame_[6];

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
  double r_[3];
  double u_[3];
  double v_[3];
  double w_[3];

public:
  Centrosome(unsigned long seed);
  void Init(centrosome_parameters *sparams);
  void SetPosition(double const *const new_pos) {}
  void Draw(std::vector<graph_struct *> &graph_array) {}

  void GetInteractors(std::vector<Object *> &ixors) {}

  double GetR(int i_dim) { return r_[i_dim]; }
  double GetU(int i_dim) { return u_[i_dim]; }
  double GetV(int i_dim) { return v_[i_dim]; }
  double GetW(int i_dim) { return w_[i_dim]; }
  void UpdatePosition() {}
  void UpdatePosition(double *r_new, double *u_new, double *v_new,
                      double *w_new) {

    for (int i{0}; i < 3; i++) {
      position_[i] = r_[i] = r_new[i];
      orientation_[i] = u_[i] = u_new[i];
      v_[i] = v_new[i];
      w_[i] = w_new[i];
    }
  }
};

#endif // _CGLASS_CENTROSOME_H_
