#include "cglass/receptor.hpp"

Receptor::Receptor(unsigned long seed) : Site(seed) {
  SetSID(species_id::receptor);
}

void Receptor::SetParameters() {
  /* Set parameters from struct*/
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = sparams_->length;
  diameter_ = sparams_->diameter;
}

void Receptor::Init(receptor_parameters *sparams) {
  sparams_ = sparams;
  SetParameters();
  interactors_.push_back(this);
}
