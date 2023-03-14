#include "cglass/centrosome.hpp"
#include <iostream>

size_t Centrosome::i_spb_ = 0;

Centrosome::Centrosome(unsigned long seed) : Object(seed) {
  SetSID(species_id::centrosome);
  // FIXME anchors and filament insertion when??
  printf("NEW centrosome\n");
}

void Centrosome::Init(centrosome_parameters *sparams) {

  int n_dim{params_->n_dim};
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
  diffusion_rot_ = noise_rot_ * sqrt(8.0 * CUBE(diameter_) / delta_);
  Logger::Trace("Inserting object %d randomly", GetOID());
  printf("I_SPB: %zu\n", i_spb_);
  double sign{i_spb_++ == 0 ? 1.0 : -1.0};
  // Insert centrosome
  double pos[3] = {0, sign * params_->system_radius, 0}; // position of C.O.M.
  double u[3] = {0, sign * -1.0,
                 0};       // orientation; normal points to boundary center
  double v[3] = {0, 0, 0}; // aux vector; defines plane of SPB together with w
  double w[3] = {0, 0, 0}; // aux vector; defines plane of SPB together with v
  // SF: not sure if this is general or b/c we initialize aligned to zhat
  double xhat[3] = {1.0, 0.0, 0.0};
  double yhat[3] = {0.0, 1.0, 0.0};
  if (1.0 - ABS(u[0]) > 1e-2) {
    cross_product(u, xhat, v, 3);
  } else {
    cross_product(u, yhat, v, 3);
  }
  double norm_factor{sqrt(1.0 / dot_product(3, v, v))};
  for (int i_dim{0}; i_dim < n_dim; i_dim++) {
    v[i_dim] *= norm_factor;
  }
  cross_product(u, v, w, 3);
  InsertAt(pos, u);
  UpdatePosition(pos, u, v, w);

  // Insert associated anchors
  n_filaments_ = sparams->num_filaments_ea; // FIXME: change to n_anchors
  // Initialize AnchorSite properties
  anchors_.resize(n_filaments_);
  for (int i_anchor{0}; i_anchor < n_filaments_; i_anchor++) {
    // Set basic physical properties -- SF TODO put in YAML files
    anchors_[i_anchor].k_ = 100;
    anchors_[i_anchor].kr_ = 1000;
    anchors_[i_anchor].r0_ = 0.5;
    // Centrosome surface is 2-D,so we can define position w/ R and phi
    // For now, distribute them uniformly along circumference of a circle
    double r_on_spb = 1.25;
    double delta_phi = 2.0 * M_PI / double(n_filaments_);
    double base_pos[3];
    for (int i_dim{0}; i_dim < n_dim; i_dim++) {
      base_pos[i_dim] = GetPosition()[i_dim] +
                        r_on_spb * sin(i_anchor * delta_phi) * v_[i_dim] +
                        r_on_spb * cos(i_anchor * delta_phi) * w_[i_dim];
    }
    // Set anchor absolute position and orientation (same as SPB for now)
    for (int i_dim{0}; i_dim < n_dim; i_dim++) {
      anchors_[i_anchor].pos_[i_dim] =
          base_pos[i_dim] + anchors_[i_anchor].r0_ * u[i_dim];
      anchors_[i_anchor].u_[i_dim] = GetOrientation()[i_dim];
    }
    // Set anchor relative position (should be r0 for now ...)
    double r_rel[3];
    for (int i_dim{0}; i_dim < n_dim; i_dim++) {
      r_rel[i_dim] = anchors_[i_anchor].pos_[i_dim] - GetPosition()[i_dim];
    }
    double u_proj = dot_product(n_dim, r_rel, u_);
    double v_proj = dot_product(n_dim, r_rel, v_);
    double w_proj = dot_product(n_dim, r_rel, w_);
    for (int i_dim{0}; i_dim < n_dim; i_dim++) {
      anchors_[i_anchor].pos_rel_[i_dim] =
          u_[i_dim] * u_proj + v_[i_dim] * v_proj + w_[i_dim] * w_proj;
    }
    // no angular springs for now (mode 0 from NAB)
    double u_anchor[3];
    for (int i_dim{0}; i_dim < n_dim; i_dim++) {
      u_anchor[i_dim] = anchors_[i_anchor].u_[i_dim];
    }
    double u_proj_rel = dot_product(n_dim, u_anchor, u_);
    double v_proj_rel = dot_product(n_dim, u_anchor, v_);
    double w_proj_rel = dot_product(n_dim, u_anchor, w_);
    anchors_[i_anchor].u_rel_[0] = u_proj_rel;
    anchors_[i_anchor].u_rel_[1] = v_proj_rel;
    anchors_[i_anchor].u_rel_[2] = w_proj_rel;
    // Position and orientation set; still need to anchor w/ filament
  }
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