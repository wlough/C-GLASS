#include "cglass/receptor.hpp"

Receptor::Receptor(unsigned long seed) : Object(seed) {
  shape_ = shape::sphere;
  SetSID(species_id::receptor);
}

// Copy parameters struct as a member, set parameters, and add as an interactor
void Receptor::Init(receptor_parameters *sparams) {
  sparams_ = sparams;
  SetParameters();
  UpdatePeriodic();             // Set KMC parameters and BC's
  interactors_.push_back(this); // Receptor can interact with other objects
}

void Receptor::SetParameters() {
  color_ = sparams_->color;
  name_ = sparams_->name;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = sparams_->length;
  diameter_ = sparams_->diameter;
  induces_catastrophe_ = sparams_->induce_catastrophe;
}

void Receptor::SetSphereS(double s) { s_ = s; }

void Receptor::SetLocations(int i, double s) {
  i_ = i;
  s_ = s;
}

void Receptor::SetStepSize(double step_size) { step_size_ = step_size; }

void Receptor::SetPCSpecies(SpeciesBase *pc_species) {}

void Receptor::SetPCObject(Object *pc_object) { pc_object_ = pc_object; }

void Receptor::SetPCObjectForSphere(Object *pc_object) {
  pc_object_ = pc_object;
}

void Receptor::SetNeighbors(Receptor *prevN, Receptor *nextN) {
  prev_r_ = prevN;
  next_r_ = nextN;
}

double Receptor::GetSphereS() { return s_; }

double Receptor::GetStepSize() { return step_size_; }

Object *Receptor::GetPCObject() { return pc_object_; }

Object *Receptor::GetPCObjectForSphere() { return pc_object_; }

Receptor *Receptor::GetPlusNeighbor() { return next_r_; }

Receptor *Receptor::GetMinusNeighbor() { return prev_r_; }

// Use PointCover object positions to update
void Receptor::UpdatePosition() {
  // Update position only on the midstep to save time, since Receptors
  // currently can't lie on Filaments. If Receptors on Filaments is implemented,
  // change to update on main step & midstep if attached to Filament (which updates
  // on both), and only on main step if attached to anything else.
  if (params_->on_midstep || sparams_->stationary_flag)
    return;
  // Check that the PointCover is associated with a species
  if (pc_species_) {
    pc_species_->CalcPCPosition(i_, s_, position_);
  }
  // Rescale position for periodic BC's
  UpdatePeriodic();
}

void Receptor::CalcTorque() {
  // Calculate torque by using the length along object.
  double r_par[3] = {0, 0, 0};
  const double *o = pc_object_->GetOrientation();
  for (int i = 0; i < n_dim_; ++i) {
    r_par[i] = s_ * o[i];
  }
  cross_product(r_par, force_, torque_, 3);
}

void Receptor::AddForce(const double *const force) {
  if (pc_object_) {
    Object::AddForce(force);
    pc_object_->AddForce(force);
  }
}

void Receptor::AddTorque(const double *const torque) {
  if (pc_object_) {
    CalcTorque();
    pc_object_->AddTorque(torque_);
  }
}

void Receptor::SubForce(const double *const force) {
  if (pc_object_) {
    Object::SubForce(force);
    pc_object_->SubForce(force);
  }
}

void Receptor::SubTorque(const double *const torque) {
  if (pc_object_) {
    CalcTorque();
    pc_object_->SubTorque(torque_);
  }
}

void Receptor::ReadSpec(std::fstream &ispec) {
  ispec.read(reinterpret_cast<char *>(&n_anchored_), sizeof(int));
}

void Receptor::WriteSpec(std::fstream &ospec) {
  ospec.write(reinterpret_cast<char *>(&n_anchored_), sizeof(int));
}

void Receptor::ConvertSpec(std::fstream &ispec, std::fstream &otext) {
  if (ispec.eof())
    return;
  int n_anchored;
  ispec.read(reinterpret_cast<char *>(&n_anchored), sizeof(int));
  otext << n_anchored << std::endl;
}

void Receptor::WriteSpecTextHeader(std::fstream &otext) {
  otext << "n_anchored" << std::endl;
}
