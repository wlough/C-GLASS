#include "cglass/chromosome.hpp"
#include <iostream>

Chromosome::Chromosome(unsigned long seed) : Object(seed) {
  SetSID(species_id::chromosome);
}

void Chromosome::Init(chromosome_parameters *sparams) {
  sparams_ = sparams;
  zero_temperature_ = sparams_->zero_temperature;
  color_ = sparams->color;
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

void Chromosome::UpdatePosition() {
  if (!zero_temperature_) {
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
  }
}