#include "cglass/receptor.hpp"

Receptor::Receptor(unsigned long seed) : Sphere(seed) {
  SetSID(species_id::receptor);
}

// Set parameters from sparams struct
void Receptor::SetParameters() {
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = sparams_->length;
  diameter_ = sparams_->diameter;
}

// Copy parameters struct as a member, set parameters, and add as an interactor
void Receptor::Init(receptor_parameters *sparams) {
  sparams_ = sparams;
  SetParameters();
  UpdatePeriodic(); // Set KMC parameters and BC's
  interactors_.push_back(this); // Receptor can interact with other objects
}

// i/o functions- Write/read species file and convert to text
void Receptor::WriteSpec(std::fstream &ospec) {
  ospec.write(reinterpret_cast<char*>(&n_anchored_), sizeof(int));
}

void Receptor::WriteSpecTextHeader(std::fstream &otext) {
  otext << "n_anchored" << std::endl;
}

void Receptor::ConvertSpec(std::fstream &ispec, std::fstream &otext) {
  if (ispec.eof()) return;
  int n_anchored;
  ispec.read(reinterpret_cast<char*>(&n_anchored), sizeof(int));
  otext << n_anchored << std::endl;
}

void Receptor::ReadSpec(std::fstream &ispec) {
  ispec.read(reinterpret_cast<char*>(&n_anchored_), sizeof(int));
}
