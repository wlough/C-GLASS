#include "cglass/chromosome.hpp"
#include <iostream>

Chromosome::Chromosome(unsigned long seed) : Object(seed) {
  SetSID(species_id::chromosome);
}

void Chromosome::Init(chromosome_parameters *sparams) {
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
  Object::InsertRandom();
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
      for (int i = 0; i < n_dim_; ++i) {
        double kick = rng_.RandomUniform() - 0.5;
        force_[i] += kick * diffusion_;
        // printf("force_[%i] = %g\n", i, force_[i]);
      }
      double dr[3];
      double fmag = 0;
      for (int i = 0; i < n_dim_; ++i) {
        fmag += force_[i] * force_[i];
        dr[i] = force_[i] * delta_ * gamma_trans_;
        position_[i] += dr[i];
      }
    }
    if (diffusion_rot_ > 0) {
      // Now handle the random orientation update
      for (int j = 0; j < n_dim_ - 1; ++j) {
        double mag = rng_.RandomNormal(diffusion_rot_);
        for (int i = 0; i < n_dim_; ++i) {
          orientation_[i] += mag * body_frame_[n_dim_ * j + i];
        }
      }
      normalize_vector(orientation_, n_dim_);
    }
  }
}