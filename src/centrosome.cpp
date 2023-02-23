#include "cglass/centrosome.hpp"
#include <iostream>

Centrosome::Centrosome(unsigned long seed) : Object(seed) {
  SetSID(species_id::centrosome);
  // FIXME anchors and filament insertion when??
}

void Centrosome::Init(centrosome_parameters *sparams) {

  printf("CALLEd IT\n");

  sparams_ = sparams;
  zero_temperature_ = sparams_->zero_temperature;
  color_ = sparams->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = sparams->length;
  diameter_ = sparams->diameter;
  gamma_trans_ = 1.0 / (diameter_);
  gamma_rot_ = 3.0 / CUBE(diameter_);
  printf("diameter: %g\n", diameter_);
  noise_tr_ = sparams->translational_noise;
  printf("noise: %g\n", noise_tr_);
  noise_rot_ = sparams->rotational_noise;
  diffusion_ = noise_tr_ * sqrt(24.0 * diameter_ / delta_);
  printf("diffusion: %g\n\n", diffusion_);
  diffusion_rot_ = noise_rot_ * sqrt(8.0 * CUBE(diameter_) / delta_);
  Logger::Trace("Inserting object %d randomly", GetOID());
  //   double u[3] = {0, 0, 0};
  //   rng_.RandomUnitVector(n_dim_, u);
  double r{params_->system_radius};
  //   double pos[3] = {};
  double u[3] = {1, 0, 0}; // Initially align with x
  double pos[3] = {r, 0, 0};
  InsertAt(pos, u);

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