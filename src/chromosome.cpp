#include "cglass/chromosome.hpp"
#include <iostream>

Chromosome::Chromosome(unsigned long seed) : Object(seed) {
  SetSID(species_id::chromosome);
  sisters_.emplace_back(Chromatid(seed));
  sisters_.emplace_back(Chromatid(seed));
}

void Chromosome::Init(chromosome_parameters *sparams) {

  sparams_ = sparams;
  zero_temperature_ = sparams_->zero_temperature;
  color_ = sparams->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = 0.0; //sparams->length;
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
  double pos[3] = {0, 0, 0};
  double u[3] = {0, 0, 0};
  Logger::Trace("Inserting object %d randomly", GetOID());
  rng_.RandomCoordinate(space_, pos, diameter_);
  rng_.RandomUnitVector(n_dim_, u);
  for (int i_sis{0}; i_sis < sisters_.size(); i_sis++) {
    double pos_sis[3]{pos[0], pos[1], pos[2]};
    for (int i{0}; i < 3; i++) {
      double sign{i_sis == 0 ? -1.0 : 1.0};
      pos_sis[i] += 0.8 * sign * diameter_;
    }
    sisters_[i_sis].Init(sparams);
    sisters_[i_sis].InsertAt(pos_sis, u);
    printf("chromatid %i inserted @ (%g, %g, %g)\n", i_sis,
           sisters_[i_sis].GetPosition()[0], sisters_[i_sis].GetPosition()[1],
           sisters_[i_sis].GetPosition()[2]);
  }
  InsertAt(pos, u);
}

void Chromosome::GetBodyFrame() {
  if (n_dim_ == 2) {
    body_frame_[0] = orientation_[1];
    body_frame_[1] = -orientation_[0];
  } else {
    double vect1[3] = {1.0, 0.0, 0.0};
    double vect2[3] = {0.0, 1.0, 0.0};
    if (1.0 - ABS(orientation_[0]) > 1e-2)
      cross_product(orientation_, vect1, &(body_frame_[0]), n_dim_);
    else
      cross_product(orientation_, vect2, &(body_frame_[0]), n_dim_);
    normalize_vector(&(body_frame_[0]), n_dim_);
    cross_product(orientation_, &(body_frame_[0]), &(body_frame_[3]), n_dim_);
  }
}

void Chromosome::UpdatePosition() {
  if (!zero_temperature_) {
    GetBodyFrame();
    // printf("boink\n");
    if (diffusion_ > 0) {
      // printf("bonk\n");
      double kick_force[n_dim_];
      for (int i = 0; i < n_dim_; ++i) {
        double kick = rng_.RandomUniform() - 0.5;
        kick_force[i] = kick * diffusion_;
        force_[i] += kick_force[i];
        for (auto &&sis : sisters_) {
          sis.force_[i] += kick_force[i];
        }
        // printf("force_[%i] = %g\n", i, force_[i]);
      }
      double dr[3];
      double fmag = 0;
      for (int i = 0; i < n_dim_; ++i) {
        fmag += force_[i] * force_[i];
        dr[i] = force_[i] * delta_ * gamma_trans_;
        position_[i] += dr[i];
        for (auto &&sis : sisters_) {
          sis.position_[i] += dr[i];
        }
      }
    }
    if (diffusion_rot_ > 0) {
      // Now handle the random orientation update
      for (int j = 0; j < n_dim_ - 1; ++j) {
        double mag = rng_.RandomNormal(diffusion_rot_);
        for (int i = 0; i < n_dim_; ++i) {
          orientation_[i] += mag * body_frame_[n_dim_ * j + i];
          for (auto &&sis : sisters_) {
            sis.orientation_[i] += mag * body_frame_[n_dim_ * j + i];
          }
        }
      }
      normalize_vector(orientation_, n_dim_);
      for (auto &&sis : sisters_) {
        normalize_vector(sis.orientation_, n_dim_);
      }
    }
  }
}