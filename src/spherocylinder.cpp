#include "cglass/spherocylinder.hpp"

Spherocylinder::Spherocylinder(unsigned long seed) : BrRod(seed) {
  SetSID(species_id::spherocylinder);
}

void Spherocylinder::SetParameters() {
  color_ = sparams_->color;
  name_ = sparams_->name;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  diameter_ = sparams_->diameter;
  length_ = sparams_->length;
  std::fill(body_frame_, body_frame_ + 6, 0.0);
  SetDiffusion();
}

void Spherocylinder::Init(spherocylinder_parameters *sparams) {
  sparams_ = sparams;
  SetParameters();
  InsertSpherocylinder();
  interactors_.push_back(this);
}

void Spherocylinder::InsertSpherocylinder() {
  InsertRod(sparams_->insertion_type);
}

void Spherocylinder::UpdatePosition() {
  SetPrevPosition(position_);
  ApplyForcesTorques();
  if (!params_->on_midstep && !sparams_->stationary_flag)
    Integrate();
    //position_[0]-=0.001;
  UpdatePeriodic();
}

// TODO: Interacting spheros
void Spherocylinder::ApplyForcesTorques() {}

