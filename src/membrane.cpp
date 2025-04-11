#include <cglass/membrane.hpp>
#include <filesystem>

Membrane::Membrane(unsigned long seed) : Mesh(seed) {
  printf("NEW membrane\n");
  //   throw std::runtime_error("Membrane::Membrane");
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
  ply_path_ = sparams_->ply_path;
  use_surface_tension_constant_ = sparams_->use_surface_tension_constant;
  use_surface_tension_penalty_local_ =
      sparams_->use_surface_tension_penalty_local;
  surface_tension_constant_ = sparams_->surface_tension_constant;
  dt0_ = delta_;
  radius_vertex_ = sparams_->radius_vertex;
  //   diameter_ = 2 * radius_vertex_;
}

void Membrane::Init(membrane_parameters *sparams) {
  sparams_ = sparams;
  SetParameters();
  LoadPly();
  init();
  SyncWithMats();
  for (auto &&v : vrts_) {
    v.SetDiameter(2 * radius_vertex_);
    v.SetColor(1.5, draw_type::fixed);
  }
}

void Membrane::Draw(std::vector<graph_struct *> &graph_array) {
  //   throw std::runtime_error("Membrane::Draw");
  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    vrts_[i_vrt].Draw(graph_array);
  }
}

void Membrane::LoadPly() {
  // printf("Loading ply file\n");
  printf("  Loading ply file %s\n", ply_path_.c_str());
  std::filesystem::path path(ply_path_);
  std::string directory = path.parent_path().string();
  std::string filename = path.filename().string();
  meshbrane::Membrane m = meshbrane::Membrane(ply_path_);
  xyz_coord_V_ = m.xyz_coord_V_;
  h_out_V_ = m.h_out_V_;
  v_origin_H_ = m.v_origin_H_;
  h_next_H_ = m.h_next_H_;
  h_twin_H_ = m.h_twin_H_;
  f_left_H_ = m.f_left_H_;
  h_right_F_ = m.h_right_F_;
  h_negative_B_ = m.h_negative_B_;
  //   update_vef_from_he();
  //   printf("  Saving mesh data to output/%s\n", filename.c_str());
  //   write_he_ply("output/" + filename);
}

void Membrane::SyncWithMats() {
  //   printf("SyncWithMats()\n");
  size_t num_vertices = get_num_vertices();
  //   size_t num_faces = get_num_faces();
  //   size_t num_edges = get_num_edges();
  //   size_t num_half_edges = get_num_half_edges();
  //   size_t num_boundaries = get_num_boundaries();

  // delete existing data
  vrts_.clear();

  // reserve space for new data
  vrts_.reserve(num_vertices);
  // Define vertices
  // * Assign: index_, pos_
  // printf("  Initializing vertices\n");
  for (size_t _v = 0; _v < num_vertices; _v++) {
    vrts_.emplace_back(_v, xyz_coord_v(_v));
    vrts_.back().SetDiameter(2 * radius_vertex_);
  }

  //   printf("Done SyncWithMats()\n");
}

void Membrane::UpdatePosition() {
  // printf("UpdatePositions()\n");
  //   xyz_coord_V_ *= 1.001;
  //   printf("  UpdatePositions()\n");
  //   printf("  t_= %f\n", t_);
  //   printf("  dt0_ %f\n", dt0_);
  //   printf("  dt_ %f\n", dt_);
  evolve_until(t_ + dt0_);
  SyncWithMats();
}