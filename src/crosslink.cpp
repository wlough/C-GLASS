#include "cglass/crosslink.hpp"
#include <iostream>

Crosslink::Crosslink(unsigned long seed) : Object(seed) {
  SetSID(species_id::crosslink);
}

void Crosslink::Init(crosslink_parameters *sparams) {
  sparams_ = sparams;
  name_ = sparams_->name;
  length_ = -1;
  diameter_ = sparams_->tether_diameter;
  color_ = sparams_->tether_color;
  draw_ = draw_type::_from_string(sparams_->tether_draw_type.c_str());
  rest_length_ = sparams_->rest_length;
  static_flag_ = sparams_->static_flag;
  k_spring_ = sparams_->k_spring;
  k_spring_compress_ =
      sparams_->k_spring_compress; // spring const for compression if
                                   // asymmetric_spring_flag is true
  // asymmetric_spring_flag_ = sparams_->asymmetric_spring_flag;
  static_flag_ = sparams_->static_flag;
  k_align_ = sparams_->k_align;
  // rcapture_ = sparams_->r_capture;
  bind_site_density_ = sparams_->bind_site_density;
  e_dep_factor_ = sparams_->energy_dep_factor;
  fdep_length_ = sparams_->force_dep_length;
  polar_affinity_ = sparams_->polar_affinity;
  use_bind_file_ = sparams_->anchors[0].bind_file.compare("none");
  
  Anchor anchor1(rng_.GetSeed());
  Anchor anchor2(rng_.GetSeed());
  anchors_.push_back(anchor1);
  anchors_.push_back(anchor2);
  anchors_[0].Init(sparams_, 0);
  anchors_[1].Init(sparams_, 1);
  SetSingly(bound_anchor_);
  Logger::Trace("Initializing crosslink %d with anchors %d and %d", GetOID(),
                anchors_[0].GetOID(), anchors_[1].GetOID());
}

void Crosslink::InitInteractionEnvironment(LookupTable *lut, Tracker *tracker, 
                                           std::map<Sphere *, std::pair<std::vector<double>, 
                                           std::vector<Anchor*> > > *bound_curr) { 
  lut_ = lut;
  tracker_ = tracker;
  bound_curr_ = bound_curr;
}

/* Function used to set anchor[0] position etc to xlink position etc */
void Crosslink::UpdatePosition() {}

void Crosslink::GetAnchors(std::vector<Object *> &ixors) {
  if (IsUnbound())
    return;
  if (!static_flag_) {
    ixors.push_back(&anchors_[bound_anchor_]);
  }
  if (IsDoubly()) {
    ixors.push_back(&anchors_[(int)!bound_anchor_]);
  }
}

/* Perform kinetic monte carlo step of protein with 1 head attached. */
void Crosslink::SinglyKMC() {

  double roll = rng_.RandomUniform();
  int head_bound = 0;
  // Set up KMC objects and calculate probabilities
  double unbind_prob = anchors_[bound_anchor_].GetOffRate() * delta_;
  if (static_flag_) {
    unbind_prob = 0;
  }
  tracker_->TrackSU(unbind_prob);
  int n_neighbors_rod = anchors_[bound_anchor_].GetNNeighborsRod();
  int n_neighbors_sphere = anchors_[bound_anchor_].GetNNeighborsSphere();
  int n_neighbors = n_neighbors_rod + n_neighbors_sphere;
 
  /* Initialize KMC calculation */
  KMC<Rod, Sphere> kmc_bind(anchors_[bound_anchor_].pos, n_neighbors_rod, n_neighbors_sphere, delta_, lut_);

  /* Initialize periodic boundary conditions */
  kmc_bind.SetPBCs(n_dim_, space_->n_periodic, space_->unit_cell);

  /* Calculate probability to bind */
  double kmc_bind_prob = 0;
  std::vector<double> bind_factors(n_neighbors);
  const std::vector<const Rod*>& rod_nbr_list = anchors_[bound_anchor_].GetNeighborListMemRods();
  const std::vector<const Sphere*>& sphere_nbr_list = anchors_[bound_anchor_].GetNeighborListMemSpheres();
  if (use_bind_file_) {
    for (int i = 0; i < rod_nbr_list.size(); ++i) {
      std::string name = rod_nbr_list[i]->GetName();
      bind_factors[i] = bind_param_map_->at((int)!bound_anchor_)[name].k_on_d
                        * bind_param_map_->at((int)!bound_anchor_)[name].bind_site_density;
    }
    for (int i = 0; i < sphere_nbr_list.size(); ++i) {
      std::string name = sphere_nbr_list[i]->GetName();
      bind_factors[rod_nbr_list.size() + i] =
               bind_param_map_->at((int)!bound_anchor_)[name].k_on_d
               * bind_param_map_->at((int)!bound_anchor_)[name].bind_site_density;
    }
  } else {
  double bind_factor_rod = anchors_[(int)!bound_anchor_].GetOnRate() * bind_site_density_;
  double bind_factor_sphere = anchors_[(int)!bound_anchor_].GetOnRate() * bind_site_density_;

  /* Fill vector of bind factors with rod factors, then sphere factors */
    std::fill(bind_factors.begin(), bind_factors.begin() + n_neighbors_rod, bind_factor_rod);
    std::fill(bind_factors.begin() + n_neighbors_rod, bind_factors.end(), bind_factor_sphere);
  }
 
  if (n_neighbors > 0) {
    if (!static_flag_ && polar_affinity_ != 1.0) {
      anchors_[bound_anchor_].CalculatePolarAffinity(bind_factors);
    }
    /* Use auto-filter populated with 1's for every neighbor.
    We already guarantee uniqueness, so we won't overcount. */
    kmc_bind.LUCalcTotProbsSD(anchors_[bound_anchor_].GetNeighborListMemRods(), 
                              anchors_[bound_anchor_].GetNeighborListMemSpheres(), 
                              anchors_[bound_anchor_].GetBoundOID(), bind_factors); 
    kmc_bind_prob = kmc_bind.getTotProb();
    tracker_->TrackSD(kmc_bind_prob);
  } // Find out whether we bind, unbind, or neither.
  int head_activate = choose_kmc_double(unbind_prob, kmc_bind_prob, roll);
  // Change status of activated head
  if (head_activate == 0) {
    // Unbind bound head
    // Track unbinding
    tracker_->UnbindSU();
    anchors_[bound_anchor_].Unbind();
    SetUnbound();
    Logger::Trace("Crosslink %d came unbound", GetOID());
  } else if (head_activate == 1) {
    // Bind unbound head
    // Track binding
    tracker_->BindSD();
    /* Position on rod where protein will bind with respect to center of rod,
     * passed by reference */
    double bind_lambda;
    /* Find out which rod we are binding to */
    int i_bind = kmc_bind.whichObjBindSD(bind_lambda, roll);
    if (i_bind < 0) {
      printf("i_bind = %d\nbind_lambda = %2.2f\n", i_bind, bind_lambda);
      Logger::Error("kmc_bind.whichRodBindSD in Crosslink::SinglyKMC"
                    " returned an invalid result!");
    }
    if (i_bind < n_neighbors_rod) {
      Rod *bind_obj = anchors_[bound_anchor_].GetRodNeighbor(i_bind);
      double obj_length = bind_obj->GetLength();
      /* KMC returns bind_lambda to be with respect to center of rod. We want 
      it to be specified from the tail of the rod to be consistent */
      bind_lambda += 0.5 * obj_length;
      /* KMC can return values that deviate a very small amount from the true 
      rod length. Bind to ends if lambda < 0 or lambda > bond_length. */
      if (bind_lambda > obj_length) {
        bind_lambda = obj_length;
      } else if (bind_lambda < 0) {
        bind_lambda = 0;
      }
      anchors_[(int)!bound_anchor_].AttachObjLambda(bind_obj, bind_lambda);
      SetDoubly();
      Logger::Trace("Crosslink %d became doubly bound to obj %d", GetOID(),
                  bind_obj->GetOID());
    } else {
      Sphere *bind_obj = anchors_[bound_anchor_].GetSphereNeighbor(i_bind - n_neighbors_rod);
      (*bound_curr_)[bind_obj].first.push_back(kmc_bind.getProb(i_bind));
      (*bound_curr_)[bind_obj].second.push_back(&anchors_[(int)!bound_anchor_]);
      anchors_[(int)!bound_anchor_].AttachObjCenter(bind_obj);
      bind_obj->DecrementNAnchored(); // For knockout loop- allow collisions
      SetDoubly();
      Logger::Trace("Crosslink %d became doubly bound to obj %d", GetOID(),
                  bind_obj->GetOID());
    }
  }
}

/* Perform kinetic monte carlo step of protein with 2 heads of protein
 * object attached. */
void Crosslink::DoublyKMC() {
  /* Calculate force-dependent unbinding for each head */
  double tether_stretch = length_ - rest_length_;
  // For asymmetric springs apply second spring constant for compression
  double e_dep = e_dep_factor_;
  double f_dep = fdep_length_;
  if (k_spring_compress_ >= 0 && tether_stretch < 0) {
    e_dep *= 0.5 * k_spring_compress_ * SQR(tether_stretch);
    f_dep *= k_spring_compress_ * tether_stretch;
  } else {
    e_dep *= 0.5 * k_spring_ * SQR(tether_stretch);
    f_dep *= k_spring_ * tether_stretch;
  }
  std::vector<double> unbind_prob;
  for (int i = 0; i < 2; i++) {
    unbind_prob.push_back(anchors_[i].GetOffRate() * delta_ * exp(e_dep + f_dep));
  }
  tracker_->TrackDS(unbind_prob[0]); // Richelle modify to track full prob
  double roll = rng_.RandomUniform();
  int head_activate = -1;
  if (static_flag_) {
    head_activate = choose_kmc_double(0, unbind_prob[1], roll);
  } else {
    // Probability of unbinding follows a poisson process but assume that only
    // one head can unbind during a time step.
    head_activate = choose_kmc_double(unbind_prob[0], unbind_prob[1], roll);
  }
  if (head_activate > -1) {
    tracker_->UnbindDS();
    Logger::Trace("Doubly-bound crosslink %d came unbound from %d", GetOID(),
                  anchors_[head_activate].GetBoundOID());
    anchors_[head_activate].Unbind();
    SetSingly((int)!head_activate);
  }
}

void Crosslink::CalculateBinding() {
  if (IsSingly()) {
    SinglyKMC();
  } else if (IsDoubly()) {
    DoublyKMC();
  }
  ClearNeighbors();
}

/* Only singly-bound crosslinks interact */
void Crosslink::GetInteractors(std::vector<Object *> &ixors) {
  ClearNeighbors();
  if (IsSingly()) {
    ixors.push_back(&anchors_[bound_anchor_]);
  }
}

void Crosslink::ClearNeighbors() { anchors_[bound_anchor_].ClearNeighbors(); }

void Crosslink::UpdateAnchorsToMesh() {
  anchors_[0].UpdateAnchorPositionToMesh();
  anchors_[1].UpdateAnchorPositionToMesh();
}

void Crosslink::UpdateAnchorPositions() {
  if (!sparams_->stationary_flag) {
    anchors_[0].UpdatePosition();   
    anchors_[1].UpdatePosition();
  }
}

void Crosslink::ApplyTetherForces() {
  if (!IsDoubly())
    return;
  anchors_[0].ApplyAnchorForces();
  anchors_[1].ApplyAnchorForces();
}

void Crosslink::UpdateCrosslinkForces() {
  /* Update anchor positions in space to calculate tether forces */
  UpdateAnchorsToMesh();
  /* Check if an anchor became unbound due to diffusion, etc */
  UpdateXlinkState();
  /* If we are doubly-bound, calculate and apply tether forces */
  CalculateTetherForces();
}

void Crosslink::UpdateCrosslinkPositions() {
  /* Have anchors diffuse/walk along mesh */
  UpdateAnchorPositions();
  /* Check if an anchor became unbound do to diffusion, etc */
  UpdateXlinkState();
  /* Check for binding/unbinding events using KMC */
  CalculateBinding();
}

/* This function ensures that singly-bound crosslinks have anchor[0] bound and
   anchor[1] unbound. */
void Crosslink::UpdateXlinkState() {
  if (!anchors_[0].IsBound() && !anchors_[1].IsBound()) {
    SetUnbound();
    return;
  }
  if (IsDoubly() && !anchors_[1].IsBound()) {
    SetSingly(0);
  } else if (IsDoubly() && !anchors_[0].IsBound()) {
    SetSingly(1);
  }
  if (IsSingly() && (anchors_[1].IsBound() && anchors_[0].IsBound())) {
    SetDoubly();
  }
}

void Crosslink::ZeroForce() {
  std::fill(force_, force_ + 3, 0.0);
  anchors_[0].ZeroForce();
  anchors_[1].ZeroForce();
}

void Crosslink::CalculateTetherForces() {
  ZeroForce();
  if (!IsDoubly())
    return;
  Interaction ix(&anchors_[0], &anchors_[1]);
  MinimumDistance mindist;
  mindist.ObjectObject(ix);
  length_ = sqrt(ix.dr_mag2);
  double stretch = length_ - rest_length_;
  for (int i = 0; i < params_->n_dim; ++i) {
    orientation_[i] = ix.dr[i] / length_;
    position_[i] = ix.midpoint[i];
  }
  // If compression spring constant is set, use that when spring is compressed.
  if (k_spring_compress_ >= 0 && stretch < 0)
    tether_force_ = k_spring_compress_ * stretch;
  else
    tether_force_ = k_spring_ * stretch;

  for (int i = 0; i < params_->n_dim; ++i) {
    force_[i] = tether_force_ * orientation_[i];
  }
  anchors_[0].AddForce(force_);
  anchors_[1].SubForce(force_);

  // If one anchor induces catastrophe and the other is attached to a filament, depolymerize
  // attached filament.
  if (anchors_[0].InducesCatastrophe() && anchors_[1].AttachedToFilament()) {
    anchors_[1].InduceCatastrophe();
  } else if (anchors_[1].InducesCatastrophe() && anchors_[0].AttachedToFilament()) {
    anchors_[0].InduceCatastrophe();
  }
  
  // Update xlink's position (for drawing)
  UpdatePeriodic();
}

/* Attach a crosslink anchor to object in a random fashion */
void Crosslink::AttachObjRandom(std::pair<Object*, int> obj_index) {
  /* Attaching to random obj implies first anchor binding from solution, so
   * this crosslink should be new and should not be singly or doubly bound */
  if ((obj_index.first->GetShape() == +shape::rod) || (obj_index.first->GetShape() == +shape::sphere)) {
    int roll = rng_.RandomUniform();
    bound_anchor_ = obj_index.second;
    anchors_[obj_index.second].AttachObjRandom(obj_index.first);
    SetCompID(obj_index.first->GetCompID());
  } else {
    Logger::Error("Crosslink binding to %s shaped objects not yet implemented.", obj_index.first->GetShape()._to_string());
  }
}

void Crosslink::Draw(std::vector<graph_struct *> &graph_array) {
  /* Draw anchors */
  anchors_[0].Draw(graph_array);
  anchors_[1].Draw(graph_array);
  /* Draw tether */
  if (IsDoubly() && length_ > 0) {
    std::copy(scaled_position_, scaled_position_ + 3, g_.r);
    for (int i = space_->n_periodic; i < n_dim_; ++i) {
      g_.r[i] = position_[i];
    }
    // std::copy(position_, position_+3, g_.r);
    std::copy(orientation_, orientation_ + 3, g_.u);
    g_.color = color_;
    if (params_->graph_diameter > 0) {
      g_.diameter = params_->graph_diameter;
    } else {
      g_.diameter = diameter_;
    }
    g_.length = length_;
    g_.draw = draw_;
    graph_array.push_back(&g_);
  }
}

void Crosslink::SetDoubly() {
  state_ = bind_state::doubly;
  SetAnchorStates();
}

void Crosslink::SetSingly(int bound_anchor) {
  state_ = bind_state::singly;
  bound_anchor_ = bound_anchor;
  SetAnchorStates();
}

void Crosslink::SetUnbound() {
  state_ = bind_state::unbound;
  SetAnchorStates();
}

void Crosslink::SetAnchorStates() {
  anchors_[0].SetState(state_);
  anchors_[1].SetState(state_);
};

const bool Crosslink::IsDoubly() const { return state_ == +bind_state::doubly; }
const bool Crosslink::IsSingly() const { return state_ == +bind_state::singly; }
const bool Crosslink::IsUnbound() const {
  return state_ == +bind_state::unbound;
}

void Crosslink::WriteSpec(std::fstream &ospec) {
  if (IsUnbound()) {
    Logger::Error("Unbound crosslink tried to WriteSpec!");
  }
  bool is_doubly = IsDoubly();
  ospec.write(reinterpret_cast<char *>(&is_doubly), sizeof(bool));
  ospec.write(reinterpret_cast<char *>(&diameter_), sizeof(double));
  ospec.write(reinterpret_cast<char *>(&length_), sizeof(double));
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&position_[i]), sizeof(double));
  }
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&orientation_[i]), sizeof(double));
  }
  ospec.write(reinterpret_cast<char *>(&oid_), sizeof(int));
  anchors_[0].WriteSpec(ospec);
  anchors_[1].WriteSpec(ospec);
}

void Crosslink::WriteSpecTextHeader(std::fstream &otext) {
  otext << "is_doubly diameter length position[0] position[1] position[2] "
        << "orientation[0] orientation[1] orientation[2] oid_" << std::endl;
}

void Crosslink::ConvertSpec(std::fstream &ispec, std::fstream &otext) {
  if (ispec.eof())
    return;
  bool is_doubly;
  double diameter, length;
  double position[3], orientation[3];
  int oid;
  // Read in all data from spec file ispec
  ispec.read(reinterpret_cast<char *>(&is_doubly), sizeof(bool));
  ispec.read(reinterpret_cast<char *>(&diameter), sizeof(double));
  ispec.read(reinterpret_cast<char *>(&length), sizeof(double));
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&position[i]), sizeof(double));
  }
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&orientation[i]), sizeof(double));
  }
  ispec.read(reinterpret_cast<char *>(&oid), sizeof(int));
  // Write out data to SpecText file otext
  otext << is_doubly << " " << diameter << " " << length << " " << position[0] << " " 
        << position[1] << " " << position[2] << " " << orientation[0] 
        << " " << orientation[1] << " " << orientation[2] << " " << oid << std::endl;
  // Convert anchor data
  Anchor::WriteSpecTextHeader(otext);
  for (int i = 0; i < 2; ++i) {
    Anchor::ConvertSpec(ispec, otext);
  }
}

void Crosslink::ReadSpec(std::fstream &ispec) {
  if (ispec.eof())
    return;
  SetSingly(bound_anchor_);
  bool is_doubly;
  ispec.read(reinterpret_cast<char *>(&is_doubly), sizeof(bool));
  ispec.read(reinterpret_cast<char *>(&diameter_), sizeof(double));
  ispec.read(reinterpret_cast<char *>(&length_), sizeof(double));
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&position_[i]), sizeof(double));
  }
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&orientation_[i]), sizeof(double));
  }
  ispec.read(reinterpret_cast<char *>(&oid_), sizeof(int));
  UpdatePeriodic();
  anchors_[0].ReadSpec(ispec);
  anchors_[1].ReadSpec(ispec);
  if (is_doubly) {
    SetDoubly();
  }
}

void Crosslink::WriteCheckpoint(std::fstream &ocheck) {
  Object::WriteCheckpoint(ocheck);
  anchors_[0].WriteCheckpointHeader(ocheck);
  anchors_[1].WriteCheckpointHeader(ocheck);
}

void Crosslink::ReadCheckpoint(std::fstream &icheck) {
  Object::ReadCheckpoint(icheck);
  anchors_[0].ReadCheckpointHeader(icheck);
  anchors_[1].ReadCheckpointHeader(icheck);
  Logger::Trace("Reloading anchor from checkpoint with cid %d",
                anchors_[bound_anchor_].GetCompID());
  if (IsDoubly()) {
    Logger::Trace("Reloading anchor from checkpoint with cid %d",
                  anchors_[(int)!bound_anchor_].GetCompID());
  }
}

const double Crosslink::GetDrTot() {
  if (IsSingly()) {
    return anchors_[bound_anchor_].GetDrTot();
  } else if (IsDoubly()) {
    double dr1 = anchors_[0].GetDrTot();
    double dr2 = anchors_[1].GetDrTot();
    if (dr1 > dr2) {
      return dr1;
    } else {
      return dr2;
    }
  } else {
    return 0;
  }
}

void Crosslink::ZeroDrTot() {
  anchors_[bound_anchor_].ZeroDrTot();
  if (IsDoubly()) {
    anchors_[(int)!bound_anchor_].ZeroDrTot();
  }
}

void Crosslink::InsertAt(double const *const new_pos, double const *const u) {
  static_flag_ = true;
  anchors_[bound_anchor_].InsertAt(new_pos, u);
  anchors_[bound_anchor_].SetBound();
  anchors_[bound_anchor_].SetStatic(true);
  SetSingly(bound_anchor_);
}

void Crosslink::SetBindParamMap(std::vector<std::map<std::string, bind_params> > *bind_param_map) {
  bind_param_map_ = bind_param_map;
  anchors_[0].SetBindParamMap(bind_param_map_);
  anchors_[1].SetBindParamMap(bind_param_map_);
}

void Crosslink::SetObjSize(double *obj_size) {
  if (!obj_size) Logger::Error("Crosslink received nullptr obj_size");
  obj_size_ = obj_size;
  anchors_[0].SetObjSize(obj_size);
  anchors_[1].SetObjSize(obj_size);
}

void Crosslink::SetBindRate(double *bind_rate) {
  if (!bind_rate) Logger::Warning("Crosslink received nullptr bind_rate");
  bind_rate_ = bind_rate;
  anchors_[0].SetBindRate(bind_rate);
  anchors_[1].SetBindRate(bind_rate);
}

const double* const Crosslink::GetObjSize() {
  if (!obj_size_) Logger::Warning("Crosslink sent nullptr obj_size");
  return obj_size_;
}

const int Crosslink::GetNNeighbors() const {
  return anchors_[bound_anchor_].GetNNeighbors();
}

const double *const Crosslink::GetPosition() {
  return anchors_[bound_anchor_].GetPosition();
}

const double *const Crosslink::GetOrientation() {
  return anchors_[bound_anchor_].GetOrientation();
}
