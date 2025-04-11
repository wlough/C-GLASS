#include <cglass/membrane.hpp>

Membrane::Membrane(unsigned long seed) : Mesh(seed) {
  printf("NEW membrane\n");
  throw std::runtime_error(
      "************************************************************************"
      "********************Membrane::Membrane");
  SetSID(species_id::membrane);
}

void Membrane::SetParameters() {
  /* Read parameters from membrane parameters */
  bending_modulus_ = sparams_->bending_modulus;
  spontaneous_curvature_ = sparams_->spontaneous_curvature;
  splay_modulus_ = sparams_->splay_modulus;
  dimensionless_tether_repulsive_singularity_ =
      sparams_->dimensionless_tether_repulsive_singularity;
  dimensionless_tether_repulsive_onset_ =
      sparams_->dimensionless_tether_repulsive_onset;
  dimensionless_tether_attractive_onset_ =
      sparams_->dimensionless_tether_attractive_onset;
  dimensionless_tether_attractive_singularity_ =
      sparams_->dimensionless_tether_attractive_singularity;
  tether_stiffness_ = sparams_->tether_stiffness;
  fix_target_face_area_ = sparams_->fix_target_face_area;
  area_stiffness_ = sparams_->area_stiffness;
  fix_target_volume_ = sparams_->fix_target_volume;
  volume_stiffness_ = sparams_->volume_stiffness;
  node_drag_coefficient_ = sparams_->node_drag_coefficient;
  kBT_ = sparams_->kBT;
  enable_flipping_ = sparams_->enable_flipping;
  enable_fluctuations_ = sparams_->enable_fluctuations;
  dt_flip_ = sparams_->dt_flip;
  flipping_probability_ = sparams_->flipping_probability;
}

void Membrane::Init(membrane_parameters *sparams) {
  sparams_ = sparams;
  SetParameters();
}