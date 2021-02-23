#include "cglass/crosslink_species.hpp"

CrosslinkSpecies::CrosslinkSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::crosslink);
}
void CrosslinkSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  k_on_ = sparams_.k_on_s;
  linear_bind_site_density_ = sparams_.linear_bind_site_density;
  surface_bind_site_density_ = sparams_.surface_bind_site_density;
  begin_with_bound_crosslinks_ = sparams_.begin_with_bound_crosslinks;
  xlink_concentration_ = sparams_.concentration;
  infinite_reservoir_flag_ = sparams_.infinite_reservoir_flag;
  sparams_.num = (int)round(sparams_.concentration * space_->volume);
}

void CrosslinkSpecies::LoadBindingSpecies() {
  YAML::Node node = YAML::LoadFile(sparams_.bind_file);
}
void CrosslinkSpecies::AddMember() {
  Species::AddMember();
  members_.back().InitInteractionEnvironment(&lut_, tracker_, bound_curr_);
  members_.back().SetObjArea(obj_area_);
  *update_ = true;
}

void CrosslinkSpecies::InitInteractionEnvironment(
    std::vector<Object *> *objs, double *obj_len, double *obj_area,
    Tracker *tracker, bool *update,
    std::map<Sphere *, std::pair<std::vector<double>, std::vector<Anchor *>>>
        *bound_curr) {
  objs_ = objs;
  obj_length_ = obj_len;
  obj_area_ = obj_area;
  bound_curr_ = bound_curr;
  update_ = update;
  tracker_ = tracker;
  LUTFiller *lut_filler_ptr = MakeLUTFiller();
  lut_ = LookupTable(lut_filler_ptr, sparams_.use_binding_volume);
  TestKMCStepSize();
  delete lut_filler_ptr;
}

/*! \brief Test simulation step size to make sure KMC will work properly.
 *
 * \return void, will end program if dt is too high.
 */
void CrosslinkSpecies::TestKMCStepSize() {
  /* TODO: Make this acceptable for C-GLASS <24-06-20, ARL> */
  KMC<int, char> kmc_diag(params_->delta, &lut_);

  const double prob_thresh = 1e-4;

  double probs[4];
  double k_on_s = sparams_.k_on_s;
  double k_on_d = sparams_.k_on_d;
  double k_off_s = sparams_.k_off_s;
  double k_off_d = sparams_.k_off_d;

  // Constant rate factors for (un)binding. Change if bind model changes.
  double u_s_fact = xlink_concentration_ * k_on_s *
                    linear_bind_site_density_; //per length basis
  double s_u_fact = k_off_s;
  double s_d_fact = k_on_d * linear_bind_site_density_;
  double d_s_fact = k_off_d;

  if (sparams_.energy_dep_factor == 0 && sparams_.force_dep_length == 0) {
    kmc_diag.Diagnostic(u_s_fact, s_u_fact, s_d_fact, d_s_fact, probs);
  } else {
    kmc_diag.DiagnosticUnBindDep(u_s_fact, s_u_fact, s_d_fact, probs);
  }

  // If any probabilities of a double binding event occuring is too high,
  // throw a warning.
  bool throw_error = false;
  if (probs[0] > prob_thresh || std::isnan(probs[0])) {
    Logger::Warning(
        "Probability of double event (U->S->U = %f) is too high. Try "
        "decreasing delta, concentration, or singly (un)binding parameters.",
        probs[0]);
    throw_error = true;
  }
  if (probs[1] > prob_thresh || std::isnan(probs[1])) {
    Logger::Warning(
        "Probability of double event (U->S->D = %f) is too high. Try "
        "decreasing delta, concentration, or binding parameters.",
        probs[1]);
    throw_error = true;
  }
  if (probs[2] > prob_thresh || std::isnan(probs[2])) {
    Logger::Warning("Probability of double event (S->D->S = %f) is too high. "
                    "Try decreasing delta or doubly (un)binding parameters.",
                    probs[2]);
    throw_error = true;
  }
  if (probs[3] > prob_thresh || std::isnan(probs[3])) {
    Logger::Warning(" !!!Probability of double event (D->S->U = %f) is too "
                    "high. Try decreasing delta or unbinding parameters.",
                    probs[3]);
    throw_error = true;
  }

  if (throw_error) {
    Logger::Error("The likelyhood of a double KMC event is too high. Fix "
                  "before continuing.");
  }
}

LUTFiller *CrosslinkSpecies::MakeLUTFiller() {
  int grid_num = sparams_.lut_grid_num;
  if (sparams_.k_spring_compress >= 0) {
    Logger::Warning("!!!Asymmetric springs are being used. This has not been "
                    "fully tested. Use at your own risk!");
    LUTFillerAsym *lut_filler_ptr = new LUTFillerAsym(grid_num, grid_num);
    lut_filler_ptr->Init(sparams_.k_spring_compress, sparams_.k_spring,
                         sparams_.energy_dep_factor, sparams_.force_dep_length,
                         sparams_.rest_length, 1);
    return lut_filler_ptr;

  } else if (sparams_.force_dep_length == 0) {
    LUTFillerEdep *lut_filler_ptr = new LUTFillerEdep(grid_num, grid_num);
    lut_filler_ptr->Init(sparams_.k_spring * .5 *
                             (1. - sparams_.energy_dep_factor),
                         sparams_.rest_length, 1);
    return lut_filler_ptr;
  } else {
    LUTFillerFdep *lut_filler_ptr = new LUTFillerFdep(grid_num, grid_num);
    lut_filler_ptr->Init(sparams_.k_spring, sparams_.energy_dep_factor,
                         sparams_.force_dep_length, sparams_.rest_length, 1);
    return lut_filler_ptr;
  }
}

void CrosslinkSpecies::InsertCrosslinks() {
  // Random insertion (default) implies crosslinks in solution
  if (sparams_.insertion_type.compare("random") == 0) {
    sparams_.static_flag = false;
    return;
  } else if (sparams_.insertion_type.compare("centered") == 0) {
    sparams_.num = 1;
    sparams_.static_flag = true;
    sparams_.infinite_reservoir_flag = false;
    AddMember();
    double pos[3] = {0, 0, 0};
    double u[3] = {0, 0, 0};
    u[params_->n_dim - 1] = 1.0;
    members_.back().InsertAt(pos, u);
  } else if (sparams_.insertion_type.compare("random_grid") == 0) {
    sparams_.static_flag = true;
    sparams_.infinite_reservoir_flag = false;
    if (params_->n_dim == 3) {
      if (space_->type == +boundary_type::none ||
          space_->type == +boundary_type::box) {
        sparams_.num = (int)round(4 * space_->radius * space_->radius *
                                  xlink_concentration_);
      } else if (space_->type == +boundary_type::sphere) {
        sparams_.num = (int)round(M_PI * space_->radius * space_->radius *
                                  xlink_concentration_);
      } else if (space_->type == +boundary_type::budding) {
        double R = space_->radius;
        double r = space_->bud_radius;
        double d = space_->bud_height;
        sparams_.num = (int)round(
            M_PI * SQR(R) + M_PI * SQR(r) -
            SQR(r) * acos((SQR(d) + SQR(r) - SQR(R)) / (2 * d * r)) -
            SQR(R) * acos((SQR(d) + SQR(R) - SQR(r)) / (2 * d * R)) +
            0.5 * sqrt((R + r - d) * (r + d - R) * (R + d - r) * (R + r + d)));
      } else {
        Logger::Error("Boundary type not recognized in CrosslinkSpecies");
      }
    }
    for (int i = 0; i < sparams_.num; ++i) {
      AddMember();
      members_.back().InsertRandom();
      /* If in 3D, zero out the third dimension so the anchor is on a plane */
      if (params_->n_dim == 3) {
        double projected_pos[3] = {0};
        double same_u[3] = {0};
        const double *const pos = members_.back().GetPosition();
        const double *const u = members_.back().GetOrientation();
        for (int j = 0; j < 3; ++j) {
          projected_pos[j] = pos[j];
          same_u[j] = u[j];
        }
        projected_pos[2] = 0;
        members_.back().InsertAt(projected_pos, same_u);
      }
    }
  } else if (sparams_.insertion_type.compare("random_boundary") == 0) {
    sparams_.infinite_reservoir_flag = false;
    sparams_.num = (int)round(space_->BoundaryArea() * xlink_concentration_);
    double pos[3] = {0};
    double u[3] = {1, 0, 0};
    for (int i = 0; i < sparams_.num; ++i) {
      rng_.RandomBoundaryCoordinate(space_, pos);
      AddMember();
      members_.back().InsertAt(pos, u);
    }
  } else {
    Logger::Error("Insertion type %s not implemented yet for crosslinks",
                  sparams_.insertion_type.c_str());
  }
}

//Adds in Crosslinkers for begin_with_bound_crosslinks flag
void CrosslinkSpecies::InsertAttachedCrosslinksSpecies() {
  if (begin_with_bound_crosslinks_ <= 0) {
    return;
  }
  Crosslink xlink(rng_.GetSeed());
  xlink.Init(&sparams_);
  xlink.InitInteractionEnvironment(&lut_, tracker_, bound_curr_);
  xlink.SetSID(GetSID());
  members_.resize(begin_with_bound_crosslinks_, xlink);
  UpdateBoundCrosslinks();
  // Begin with bound crosslinks currently just implemented to start on rods
  for (int i = 0; i < begin_with_bound_crosslinks_; ++i) {
    BindCrosslink(shape::rod);
  }
}

// Calculate and bind crosslinkers from solution implicitly
void CrosslinkSpecies::CalculateBindingFree() {
  /* Static crosslinks are never free */
  if (sparams_.static_flag) {
    return;
  }
  /* Check crosslink binding */
  double free_concentration;
  if (infinite_reservoir_flag_) { // Have a constant concentration of
                                  // crosslinkers binding from solution
    free_concentration = xlink_concentration_;
  } else { // Have a constant number of crosslinkers in a space
    free_concentration = (sparams_.num - n_members_) / space_->volume;
  }
  double expected_lin_bind_n = linear_bind_site_density_ * free_concentration *
                               (*obj_length_) * k_on_ * params_->delta;
  double expected_surf_bind_n = surface_bind_site_density_ *
                                free_concentration * (*obj_area_) * k_on_ *
                                params_->delta;
  int linear_bind_num = rng_.RandomPoisson(expected_lin_bind_n);
  int surface_bind_num = rng_.RandomPoisson(expected_surf_bind_n);
  // Track US probabilities
  tracker_->TrackUS(expected_lin_bind_n + expected_surf_bind_n);
  tracker_->BindUS(linear_bind_num + surface_bind_num);
  // Use a Poisson distribution to calculate the number of particles
  // binding from distribution
  for (int i = 0; i < linear_bind_num; ++i) {
    /* Create a new crosslink and bind an anchor to a random rod
     * in the system */
    BindCrosslink(shape::rod);
  }
  for (int i = 0; i < surface_bind_num; ++i) {
    /* Create a new crosslink and bind an anchor to a random sphere
     * in the system */
    BindCrosslink(shape::sphere);
  }
}

/* Returns a random object with selection probability proportional to object
   length */
Object *CrosslinkSpecies::GetRandomObject(shape sh) {
  double roll = rng_.RandomUniform();
  switch (sh) {
  case shape::rod:
    roll *= (*obj_length_);
    break;
  case shape::sphere:
    roll *= (*obj_area_);
    break;
  default:
    Logger::Error("Binding to object type %s not yet implemented in "
                  "CrosslinkSpecies::GetRandomObject",
                  sh._to_string());
  }
  double vol = 0;

  // Search through interactors to find an object of type type
  for (auto obj = objs_->begin(); obj != objs_->end(); ++obj) {
    if ((*obj)->GetShape() == sh) {
      switch (sh) {
      case shape::rod:
        vol += (*obj)->GetLength();
        break;
      case shape::sphere:
        if ((*obj)->GetNAnchored() == 0) {
          vol += (*obj)->GetArea();
        }
        break;
      default:
        Logger::Error("Binding to object type %s not yet implemented in "
                      "CrosslinkSpecies::GetRandomObject",
                      sh._to_string());
      }
      if (vol > roll) {
        Logger::Trace(
            "Binding free crosslink to random object: xl %d -> obj %d",
            members_.back().GetOID(), (*obj)->GetOID());
        return *obj;
      }
    }
  }
  Logger::Error("CrosslinkSpecies::GetRandomObject should never get here!");
  return nullptr;
}

/* A crosslink binds to an object from solution */
void CrosslinkSpecies::BindCrosslink(shape sh) {
  /* Create crosslink object and initialize. Crosslink will
   * initially be singly-bound. */
  AddMember();
  members_.back().AttachObjRandom(GetRandomObject(sh));
}

/* Return singly-bound anchors, for finding neighbors to bind to */
void CrosslinkSpecies::GetInteractors(std::vector<Object *> &ixors) {
  for (auto xlink = members_.begin(); xlink != members_.end(); ++xlink) {
    xlink->GetInteractors(ixors);
  }
}

/* Returns all anchors, not just singly-bound anchors. Used for reassigning
   bound anchors to objects upon a checkpoint reload */
void CrosslinkSpecies::GetAnchorInteractors(std::vector<Object *> &ixors) {
  for (auto xlink = members_.begin(); xlink != members_.end(); ++xlink) {
    xlink->GetAnchors(ixors);
  }
}

void CrosslinkSpecies::UpdatePositions() {
  /* Only do this every other step (assuming flexible filaments with midstep)
   */
  if (params_->no_midstep) {
    UpdateBoundCrosslinks();
    CalculateBindingFree();
  } else if (params_->i_step % 2 == 0) {
    /* First update bound crosslinks state and positions */
    UpdateBoundCrosslinks();
    /* Calculate implicit binding of crosslinks from solution */
    CalculateBindingFree();
  } else {
    /* Apply tether forces from doubly-bound crosslinks onto anchored objects.
       We do this every half step only, because the fullstep tether forces are
       handled by the crosslink update on the full step */
    ApplyCrosslinkTetherForces();
  }
  midstep_ = !midstep_;
  // for (auto it=members_.begin(); it!=members_.end(); ++it) {
  // it->SanityCheck();
  //}
}

void CrosslinkSpecies::UpdateObjectArea() {
  *obj_area_ = 0.0;
  for (auto obj = objs_->begin(); obj != objs_->end(); ++obj) {
    if ((*obj)->GetShape() == +shape::sphere) {
      if ((*obj)->GetNAnchored() == 0)
        *obj_area_ += (*obj)->GetArea();
    }
  }
}

void CrosslinkSpecies::UpdateBoundCrosslinks() {
  n_members_ = 0;
  /* Update anchor positions to their attached meshes and calculate anchor
     forces */
  UpdateBoundCrosslinkForces();
  /* Apply anchor forces on bound objects sequentially */
  ApplyCrosslinkTetherForces();
  /* Update anchor positions from diffusion, walking */
  UpdateBoundCrosslinkPositions();
  /* Remove crosslinks that came unbound */
  if (!sparams_.static_flag) {
    members_.erase(std::remove_if(members_.begin(), members_.end(),
                                  [](Crosslink x) { return x.IsUnbound(); }),
                   members_.end());
  }
  /* Get the number of bound crosslinks so we know what the current
     concentration of free crosslinks is */
  n_members_ = members_.size();
}

/* This must be done sequentially to avoid racy conditions when accessing
   bound object's forces */
void CrosslinkSpecies::ApplyCrosslinkTetherForces() {
  for (auto xlink = members_.begin(); xlink != members_.end(); ++xlink) {
    xlink->ApplyTetherForces();
  }
}

void CrosslinkSpecies::UpdateBoundCrosslinkForces() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  xlink_chunk_vector chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = members_.size() / max_threads;
  xlink_iterator cur_iter = members_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    xlink_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks, update_)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto xlink = chunks[i].first; xlink != chunks[i].second; ++xlink) {
        bool init_state = xlink->IsSingly();
        if (sparams_.static_flag && init_state && xlink->GetNNeighbors() == 0) {
          continue;
        }
        xlink->UpdateCrosslinkForces();
        if (xlink->IsSingly() != init_state) {
          *update_ = true;
        }
      }
    }
  }
#else
  for (xlink_iterator xlink = members_.begin(); xlink != members_.end();
       ++xlink) {
    bool init_state = xlink->IsSingly();
    if (sparams_.static_flag && init_state && xlink->GetNNeighbors() == 0) {
      continue;
    }
    xlink->UpdateCrosslinkForces();
    if (xlink->IsSingly() != init_state) {
      *update_ = true;
    }
  }
#endif
}

void CrosslinkSpecies::UpdateBoundCrosslinkPositions() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  xlink_chunk_vector chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = members_.size() / max_threads;
  xlink_iterator cur_iter = members_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    xlink_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks, update_)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto xlink = chunks[i].first; xlink != chunks[i].second; ++xlink) {
        bool init_state = xlink->IsSingly();
        if (sparams_.static_flag && init_state && xlink->GetNNeighbors() == 0) {
          continue;
        }
        xlink->UpdateCrosslinkPositions();
        /* Xlink is no longer bound, return to solution */
        if (xlink->IsUnbound()) {
          if (sparams_.static_flag) {
            Logger::Error("Static crosslinks became unbound");
          }
          *update_ = true;
          /* If a crosslink enters or leaves the singly state, we need to
           * update xlink interactors */
        } else if (xlink->IsSingly() != init_state) {
          *update_ = true;
        }
      }
    }
  }
#else
  for (xlink_iterator xlink = members_.begin(); xlink != members_.end();
       ++xlink) {
    bool init_state = xlink->IsSingly();
    if (sparams_.static_flag && init_state && xlink->GetNNeighbors() == 0) {
      continue;
    }
    xlink->UpdateCrosslinkPositions();
    /* Xlink is no longer bound, return to solution */
    if (xlink->IsUnbound()) {
      if (sparams_.static_flag) {
        Logger::Error("Static crosslinks became unbound");
      }
      *update_ = true;
      /* If a crosslink enters or leaves the singly state, we need to update
       * xlink interactors */
    } else if (xlink->IsSingly() != init_state) {
      *update_ = true;
    }
  }
#endif
}

void CrosslinkSpecies::CleanUp() { members_.clear(); }

void CrosslinkSpecies::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->Draw(graph_array);
  }
}

void CrosslinkSpecies::ReadSpecs() {
  if (ispec_file_.eof()) {
    if (HandleEOF()) {
      return;
    } else {
      Logger::Info("EOF reached in spec file for %s %s", GetSID()._to_string(),
                   GetSpeciesName().c_str());
      early_exit = true;
      return;
    }
  }
  if (!ispec_file_.is_open()) {
    Logger::Warning("ERROR. Spec file unexpectedly not open! Exiting early.");
    early_exit = true;
    return;
  }
  n_members_ = -1;
  ispec_file_.read(reinterpret_cast<char *>(&n_members_), sizeof(int));
  /* For some reason, we can't catch the EOF above. If size == -1 still, then
     we caught a EOF here */
  if (n_members_ == -1) {
    if (HandleEOF()) {
      return;
    } else {
      Logger::Info("EOF reached in spec file for %s %s", GetSID()._to_string(),
                   GetSpeciesName().c_str());
      early_exit = true;
      return;
    }
  }
  if (n_members_ == 0) {
    members_.clear();
  } else if (n_members_ != members_.size()) {
    Crosslink xlink(rng_.GetSeed());
    xlink.Init(&sparams_);
    xlink.InitInteractionEnvironment(&lut_, tracker_, bound_curr_);
    xlink.SetSID(GetSID());
    members_.resize(n_members_, xlink);
  }
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->ReadSpec(ispec_file_);
  }
}

const int CrosslinkSpecies::GetDoublyBoundCrosslinkNumber() const {
  int num = 0;
  for (const auto &xl : members_) {
    if (xl.IsDoubly())
      ++num;
  }
  return num;
}

const double CrosslinkSpecies::GetConcentration() const {
  return sparams_.concentration;
}
const double CrosslinkSpecies::GetRCutoff() const { return lut_.getLUCutoff(); }
