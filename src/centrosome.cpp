#include "cglass/centrosome.hpp"
#include <iostream>

Centrosome::Centrosome(unsigned long seed) : Object(seed) {
  SetSID(species_id::centrosome);
  // FIXME anchors and filament insertion when??
  printf("NEW centrosome\n");
}

void Centrosome::Init(centrosome_parameters *sparams) {

  sparams_ = sparams;
  zero_temperature_ = sparams_->zero_temperature;
  color_ = sparams->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = sparams->length;
  diameter_ = sparams->diameter;
  gamma_trans_ = 1.0 / (diameter_);
  gamma_rot_ = 3.0 / CUBE(diameter_);
  noise_tr_ = sparams->translational_noise;
  noise_rot_ = sparams->rotational_noise;
  diffusion_ = noise_tr_ * sqrt(24.0 * diameter_ / delta_);
  //   printf("diameter: %g\n", diameter_);
  //   printf("noise: %g\n", noise_tr_);
  //   printf("diffusion: %g\n\n", diffusion_);
  diffusion_rot_ = noise_rot_ * sqrt(8.0 * CUBE(diameter_) / delta_);
  Logger::Trace("Inserting object %d randomly", GetOID());

  n_filaments_ = sparams->num_filaments_ea;

  // We want the SPB to point directly at the center of the spherical cell wall
  double r{params_->system_radius};
  double pos[3] = {r, 0, 0};
  double u[3] = {1, 0, 0}; // Initially align with z
  InsertAt(pos, u);
  // flip orientation fector (neg. sign.) cuz we want it to point inwards
  for (int i{0}; i < 3; i++) {
    r_[i] = position_[i] = pos[i];
    u_[i] = orientation_[i] = -u[i];
  }
  // Calculate vectors that define plane of SPB
  double vect1[3] = {1.0, 0, 0.0};
  double vect2[3] = {0.0, 1.0, 0.0};
  if (1.0 - ABS(u[0]) > 1e-2)
    cross_product(u, vect1, v_, 3);
  else
    cross_product(u, vect2, v_, 3);
  double norm_factor = sqrt(1.0 / dot_product(3, v_, v_));
  for (int i = 0; i < 3; ++i)
    v_[i] *= norm_factor;
  cross_product(u, v_, w_, 3);

  // Initialize AnchorSite properties
  anchors_.resize(n_filaments_);
  for (int i_anchor{0}; i_anchor < n_filaments_; i_anchor++) {
    anchors_[i_anchor].r0 = 2.0;
    anchors_[i_anchor].k_ = 100;
    anchors_[i_anchor].kr_ = 1000;
  }

  //   for (int i_sis{0}; i_sis < sisters_.size(); i_sis++) {
  //     double pos_sis[3]{pos[0], pos[1], pos[2]};
  //     for (int i{0}; i < 3; i++) {
  //       double sign{i_sis == 0 ? -1.0 : 1.0};
  //       pos_sis[i] += 0.5 * sign * diameter_;
  //     }
  //     sisters_[i_sis].Init(sparams);
  //     sisters_[i_sis].InsertAt(pos_sis, u);
  //     printf("chromatid %i inserted @ (%g, %g, %g)\n", i_sis,
  //            sisters_[i_sis].GetPosition()[0], sisters_[i_sis].GetPosition()[1],
  //            sisters_[i_sis].GetPosition()[2]);
  //   }
}

/*
void Centrosome::InitFilaments(filament_parameters *fparams) {
  // Set pointer to filament parameters
  fparams_ = fparams;
  // Create filament objects
  n_filaments_ = 0;
  filaments_.reserve(sparams_->num_filaments_ea);
  anchors_.resize(sparams_->num_filaments_ea);
  for (int i_fil{0}; i_fil < sparams_->num_filaments_ea; i_fil++) {
    Filament newfila{Filament(rng_.GetSeed())};
    filaments_.push_back(newfila);
    filaments_.back().SetSID(+species_id::filament);
    // INSERTION HAPPENS HERE
    filaments_.back().Init(fparams_);
    anchors_[i_fil].filament_ = &filaments_.back();
    n_filaments_++;
  }
  if (n_filaments_ != sparams_->num_filaments_ea) {
    printf("Error in Centrosome::InitFilaments()!\n");
    exit(1);
  }
  // Couple filaments to anchors and update their positions
  for (int i_fil{0}; i_fil < n_filaments_; i_fil++) {
    const double *const u = GetOrientation();
    const double *const r = GetPosition();
  }
}
*/