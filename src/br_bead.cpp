#include "cglass/br_bead.hpp"

BrBead::BrBead(unsigned long seed) : Object(seed) {
  SetSID(species_id::br_bead);
}

void BrBead::Init(br_bead_parameters *sparams) {
  sparams_ = sparams;
  name_ = sparams_->name;
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  diameter_ = sparams_->diameter;
  driving_factor_ = sparams_->driving_factor;
  zero_temperature_ = params_->zero_temperature;
  chiral_handedness_ = sparams_->chiral_handedness;
  alignment_interaction_ = sparams_->alignment_interaction;
  alignment_torque_ = sparams_->alignment_torque;
  noise_tr_ = sparams_->translational_noise;
  noise_rot_ = sparams_->rotational_noise;
  if (sparams_->draw_shape.compare("sphere") != 0) {
    draw_arrow_ = true;
  }
  if (sparams_->randomize_handedness) {
    chiral_handedness_ = (rng_.RandomUniform() > 0.5 ? 1 : -1);
  }
  if (sparams_->highlight_handedness) {
    draw_ = draw_type::fixed;
    color_ = (chiral_handedness_ > 0 ? sparams_->color : sparams_->color + M_PI);
  }
  driving_torque_ = sparams_->driving_torque;
  SetDiffusion();
  InsertBrBead();
}

int BrBead::GetCount() { return 1; }

void BrBead::ZeroForce() { Object::ZeroForce(); }

void BrBead::InsertBrBead() {
  if (sparams_->insertion_type.compare("random") == 0) {
    InsertRandom();
  } else if (sparams_->insertion_type.compare("random_oriented") == 0) {
    InsertRandom();
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = 1.0;
  } else if (sparams_->insertion_type.compare("centered_random") == 0) {
    std::fill(position_, position_ + 3, 0.0);
    rng_.RandomUnitVector(n_dim_, orientation_);
  } else if (sparams_->insertion_type.compare("centered_oriented") == 0) {
    std::fill(position_, position_ + 3, 0.0);
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = 1.0;
  } else if (sparams_->insertion_type.compare("custom") == 0) {
    // Nothing to do
  } else {
    Logger::Error("BrBead insertion type not recognized!");
  }
}

void BrBead::UpdatePosition() {
  SetPrevPosition(GetPosition());
  SetPrevOrientation(GetOrientation());
  ApplyForcesTorques();
  if (!params_->on_midstep && !sparams_->stationary_flag)
    Integrate();
  UpdatePeriodic();
}

void BrBead::ApplyForcesTorques() {
  // Add random thermal kick to the bead
  if (!zero_temperature_) {
    if (diffusion_ > 0) {
      for (int i = 0; i < n_dim_; ++i) {
        double kick = rng_.RandomUniform() - 0.5;
        force_[i] += kick * diffusion_;
      }
    }
    if (diffusion_rot_ > 0) {
      if (n_dim_ == 2) {
        double kick = rng_.RandomUniform() - 0.5;
        torque_[2] += diffusion_rot_ * kick;
      } else {
        double kicks[3] = {0};
        rng_.RandomUnitVector(n_dim_, kicks);
        for (int i=0; i<n_dim_; ++i) {
          torque_[i] += diffusion_rot_ * kicks[i];
        }
      }
    }
  }
  if (driving_factor_ > 0) {
    for (int i = 0; i < n_dim_; ++i) {
      force_[i] += driving_factor_ * orientation_[i];
    }
  }
  if (driving_torque_ > 0 && n_dim_ == 2) {
    torque_[2] += chiral_handedness_ * driving_torque_;
  }
  if (alignment_interaction_ && n_dim_ == 2) {
    for (auto ix = ixs_.begin(); ix != ixs_.end(); ++ix) {
      if (ix->first->pause_interaction)
        continue;
      double dr2 = ix->first->dr_mag2;
      double cp[3] = {0, 0, 0};
      double align_t = alignment_torque_;
      if (ix->second) {
        const double * const u = ix->first->obj2->GetOrientation();
        cross_product(orientation_, u, cp, n_dim_);
        int sign = SIGNOF(cp[2]);
        align_t *= sign*0.5*(1 - dot_product(n_dim_, u, orientation_));
      } else {
        const double * const u = ix->first->obj1->GetOrientation();
        cross_product(orientation_, u, cp, n_dim_);
        int sign = SIGNOF(cp[2]);
        align_t *= sign*0.5*(1 - dot_product(n_dim_, u, orientation_));
      }
      cp[2] = align_t/dr2;
      AddTorque(cp);
    }
  }
}

void BrBead::SetDiffusion() {
  gamma_trans_ = 1.0 / (diameter_);
  gamma_rot_ = 3.0 / CUBE(diameter_);
  diffusion_ = noise_tr_ * sqrt(24.0 * diameter_ / delta_);
  diffusion_rot_ = noise_rot_ * sqrt(8.0 * CUBE(diameter_) / delta_);
}

void BrBead::Translate() {
  double dr[3];
  double fmag = 0;
  for (int i = 0; i < n_dim_; ++i) {
    fmag += force_[i] * force_[i];
    dr[i] = force_[i] * delta_ * gamma_trans_;
    position_[i] += dr[i];
  }
}

void BrBead::Rotate() {
  double unit_torque[3], temp[3], r_rel[3];
  double domega, cos_domega, sin_domega, torque_mag;
  // First rotate orientation vector of sphere
  if (n_dim_ == 2) {
    domega = torque_[2] * delta_ * gamma_rot_;
    cos_domega = cos(domega);
    sin_domega = sin(domega);
    std::copy(orientation_, orientation_ + 3, temp);
    orientation_[0] = cos_domega * temp[0] - sin_domega * temp[1];
    orientation_[1] = sin_domega * temp[0] + cos_domega * temp[1];
  } else if (n_dim_ == 3) {
    torque_mag = 0.0;
    for (int i = 0; i < 3; ++i)
      torque_mag += torque_[i];
    for (int i = 0; i < 3; ++i)
      unit_torque[i] = torque_[i] / torque_mag;
    domega = torque_mag * delta_ * gamma_rot_;
    rotate_vector(orientation_, unit_torque, domega, n_dim_);
  }
  normalize_vector(orientation_, n_dim_);
}

void BrBead::Integrate() {
  Translate();
  Rotate();
}

void BrBead::GetInteractors(std::vector<Object *> &ix) { ix.push_back(this); }

void BrBead::Draw(std::vector<graph_struct *> &graph_array) {
  if (draw_arrow_) {
    if (draw_ == +draw_type::orientation) {
      g_.color = atan2(orientation_[1], orientation_[0]); 
    } else {
      g_.color = color_;
    }
    g2_.color = g_.color;
    std::copy(scaled_position_, scaled_position_ + 3, g_.r);
    std::copy(scaled_position_, scaled_position_ + 3, g2_.r);
    for (int i = space_->n_periodic; i < n_dim_; ++i) {
      g_.r[i] = position_[i];
      g2_.r[i] = position_[i];
    }
    std::copy(orientation_, orientation_ + 3, g_.u);
    std::copy(orientation_, orientation_ + 3, g2_.u);
    double theta = 0.1*M_PI;
    rotate_vector(g_.u, orientation_, theta, 2);
    rotate_vector(g2_.u, orientation_, -theta, 2);
    g_.diameter = 0.5*diameter_;
    g2_.diameter = 0.5*diameter_;
    g_.length = diameter_;
    g2_.length = diameter_;
    g_.draw = draw_type::fixed;
    g2_.draw = draw_type::fixed;
    for (int i=0; i<space_->n_periodic; ++i) {
      g_.r[i] -= 0.3*g_.length*g_.u[i]/space_->radius;
      g2_.r[i] -= 0.3*g2_.length*g2_.u[i]/space_->radius;
      g_.r[i] -= NINT(g_.r[i]);
      g2_.r[i] -= NINT(g2_.r[i]);
    }
    graph_array.push_back(&g_);
    graph_array.push_back(&g2_);
  } else {
    Object::Draw(graph_array);
  }
}
void BrBead::WriteSpec(std::fstream &ospec) {
  Logger::Trace("Writing br_bead specs, object id: %d", GetOID());
  Object::WriteSpec(ospec);
  ospec.write(reinterpret_cast<char *>(&chiral_handedness_), sizeof(int));
}

void BrBead::WriteSpecTextHeader(std::fstream &otext) {
  otext << "position[0] position[1] position[2] scaled_position[0] "
        << "scaled_position[1] scaled_position[2] orientation[0] "
        << "orientation[1] orientation[2] diameter length "
        << "chiral_handedness" << std::endl;
}

void BrBead::ConvertSpec(std::fstream &ispec, std::fstream &otext) {
  double position[3], scaled_position[3], orientation[3];
  double diameter, length;
  int chiral_handedness;
  if (ispec.eof())
    return;
  for (auto &posit : position)
    ispec.read(reinterpret_cast<char *>(&posit), sizeof(posit));
  for (auto &spos : scaled_position)
    ispec.read(reinterpret_cast<char *>(&spos), sizeof(spos));
  for (auto &u : orientation)
    ispec.read(reinterpret_cast<char *>(&u), sizeof(u));
  ispec.read(reinterpret_cast<char *>(&diameter), sizeof(diameter));
  ispec.read(reinterpret_cast<char *>(&length), sizeof(length));
  ispec.read(reinterpret_cast<char *>(&chiral_handedness), sizeof(int));
  otext << position[0] << " " << position[1] << " " << position[2] << " " 
        << scaled_position[0] << " " << scaled_position[1] << " "
        << scaled_position[2] << " " << orientation[0] << " " << orientation[1] 
        << " " << orientation[2] << " " << diameter << " " << length << " "
        << chiral_handedness << std::endl;
}

void BrBead::ReadSpec(std::fstream &ispec) {
  Object::ReadSpec(ispec);
  ispec.read(reinterpret_cast<char *>(&chiral_handedness_), sizeof(int));
  if (sparams_->highlight_handedness) {
    draw_ = draw_type::fixed;
    color_ = (chiral_handedness_ > 0 ? sparams_->color : sparams_->color + M_PI);
  }
}


