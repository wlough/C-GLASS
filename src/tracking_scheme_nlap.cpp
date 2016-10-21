// Implementation for tracking scheme with a neighbor list underneath (all pairs)

#include "tracking_scheme_nlap.h"

// Init functionality
void TrackingSchemeNeighborListAllPairs::Init(int pModuleID,
                                              space_struct *pSpace,
                                              PotentialBase *pPotentialBase,
                                              std::vector<interaction_t> *pInteractions,
                                              std::vector<SpeciesBase*> *pSpecies,
                                              std::vector<Simple*> *pSimples,
                                              std::unordered_map<int, int> *pOIDMap,
                                              YAML::Node *pNode) {
  TrackingScheme::Init(pModuleID, pSpace, pPotentialBase, pInteractions, pSpecies, pSimples, pOIDMap, pNode);

  // We also need our information about skin depth and rcut_
  YAML::Node node = *pNode; 

  // Grab the cutoff from the potential we're associated with
  rcut_ = pbase_->GetRCut();
  skin_ = node["skin"].as<double>();

  CreateTrackingScheme();
  last_time_ = std::chrono::high_resolution_clock::now();
  avg_update_time_ = 0.0;
  avg_occupancy_ = 0.0;
}

// Print functionality
void TrackingSchemeNeighborListAllPairs::Print() {
  TrackingScheme::Print();
  std::cout << "   rcut:     " << std::setprecision(16) << rcut_ << std::endl;
  std::cout << "   skin:     " << std::setprecision(16) << skin_ << std::endl;
  std::cout << "   eff cut:  " << std::setprecision(16) << sqrt(rcs2_) << std::endl;
  std::cout << "   dmax:     " << std::setprecision(16) << sqrt(half_skin2_) << std::endl;
}

// Print statistics
void TrackingSchemeNeighborListAllPairs::PrintStatistics() {
  GenerateStatistics();
  std::cout << "********\n";
  std::cout << name_ << std::endl;
  std::cout << "   {" << SIDToString(sid0_) << ", " << SIDToString(sid1_) << "}\n";
  std::cout << "   type: " << PtypeToString(type_) << std::endl;
  std::cout << "   nupdates: " << nupdates_ << std::endl;
  std::cout << "   rcut:     " << std::setprecision(16) << rcut_ << std::endl;
  std::cout << "   skin:     " << std::setprecision(16) << skin_ << std::endl;
  std::cout << "   eff cut:  " << std::setprecision(16) << sqrt(rcs2_) << std::endl;
  std::cout << "   dmax:     " << std::setprecision(16) << sqrt(half_skin2_) << std::endl;
  std::cout << "   avg time between updates: " << std::setprecision(8) << avg_update_time_/nupdates_ << " microseconds\n";
  std::cout << "   avg occpancy:             " << std::setprecision(8) << avg_occupancy_/nupdates_ << " particles\n";
}

// Generate the neighbor list tracking scheme (if needed)
void TrackingSchemeNeighborListAllPairs::CreateTrackingScheme() {
  rcut2_ = rcut_*rcut_;
  skin2_ = skin_*skin_;
  rcs2_ = skin2_ + rcut2_ + 2*rcut_*skin_;
  half_skin2_ = 0.25*skin2_;
}

// Overload the LoadSimples functionality to generate neighbor list
void TrackingSchemeNeighborListAllPairs::LoadSimples() {
  if (debug_trace) {
    std::cout << "TrackingNeighborListAllPairs::LoadSimples\n";
  }
  TrackingScheme::LoadSimples();
}

// Generate the interactions
void TrackingSchemeNeighborListAllPairs::GenerateInteractions(bool pForceUpdate) {
  nl_update_ = pForceUpdate;

  // Check accumulators
  double dr2 = 0.0;
  dr2 = spec0_->GetDrMax();
  nl_update_ = nl_update_ || (dr2 > half_skin2_);
  dr2 = spec1_->GetDrMax();
  nl_update_ = nl_update_ || (dr2 > half_skin2_);
  
  if (nl_update_) {
    LoadSimples();
    UpdateNeighborList();
    GenerateStatistics();

    // Reset accumulators
    spec0_->ZeroDr();
    spec1_->ZeroDr();
  }
  interactions_->insert(interactions_->end(), m_interactions_.begin(), m_interactions_.end());
}

// Actually update the neighbor list
void TrackingSchemeNeighborListAllPairs::UpdateNeighborList() {
  if (debug_trace) {
    std::cout << "[" << moduleid_ << "]TrackingSchemeNeighborListAllPairs Updating NL\n";
  }
  nupdates_++;
  m_interactions_.clear();

  // Find the unique rigids
  unique_rids_->clear();
  for (int i = 0; i < nmsimples_; ++i) {
    unique_rids_->insert(m_simples_[i]->GetRID());
    mneighbors_[i].clear();
  }
  maxrigid_ = *(unique_rids_->rbegin());

  // Loop over particles, generate the interactions
  rid_self_check_->clear();
  rid_self_check_->resize(maxrigid_+1);
  std::fill(rid_self_check_->begin(), rid_self_check_->end(), false);

  #ifdef ENABLE_OPENMP
  #pragma omp parallel
  #endif
  {
    int tid = 0;
    #ifdef ENABLE_OPENMP
    tid = omp_get_thread_num();
    #else
    tid = 0;
    #endif

    #ifdef ENABLE_OPENMP
    #pragma omp for schedule(runtime) nowait
    #endif
    for (int idx = 0; idx < nmsimples_ - 1; ++idx) {
      auto p1 = m_simples_[idx];
      int rid1 = p1->GetRID();

      bool should_exit = false;
      #ifdef ENABLE_OPENMP
      #pragma omp critical
      #endif
      {
        should_exit = (*rid_self_check_)[rid1];
        (*rid_self_check_)[rid1] = true;
      }
      if (should_exit) {
        continue;
      }

      std::unordered_set<int>* rid_check_local = rid_check_local_[tid];
      rid_check_local->clear();

      for (int jdx = idx; jdx < nmsimples_; ++jdx) {
        if (idx == jdx) continue;
        auto p2 = m_simples_[jdx];
        int rid2 = p2->GetRID();
        if (rid1 == rid2) continue;

        // XXX FIXME this might not be optimal, but check anyway....
        auto sid0 = p1->GetSID();
        auto sid1 = p2->GetSID();
        if (!(sid0 == sid0_ && sid1 == sid1_) &&
            !(sid1 == sid0_ && sid0 == sid1_)) continue;

        // We are guranteed for an interaction
        if (rid_check_local->count(rid2)) {
          continue;
        }
        // Only insert if we've found an interaction
        // otherwise, we miss stuff

        // Minimum distance calc
        interactionmindist idm;
        MinimumDistance(p1, p2, idm, ndim_, nperiodic_, space_);

        // Check this out
        if (idm.dr_mag2 < rcs2_) {
          // Create the interaction
          //interaction_t new_interaction;
          //new_interaction.idx_ = (*oid_position_map_)[p1->GetOID()];
          //new_interaction.jdx_ = (*oid_position_map_)[p2->GetOID()];
          //new_interaction.type_ = type_;
          //new_interaction.pot_ = pbase_;
          //new_interaction.kmc_track_module_ = moduleid_;
          //
          //// KMC specifics
          //if (type_ == ptype::kmc) {
          //  new_interaction.kmc_target_ = kmc_target_;
          //}
          //
          //#ifdef ENABLE_OPENMP
          //#pragma omp critical
          //#endif
          //{
          //  m_interactions_.push_back(new_interaction);
          //}
          neighbor_kmc_t new_neighbor;
          new_neighbor.idx_ = jdx;
          mneighbors_[idx].push_back(new_neighbor);
          rid_check_local->insert(rid2);
        }
      } // inner loop for second particle
    } // omp for schedule(runtime) nowait

    #ifdef ENABLE_OPENMP
    #pragma omp barrier
    #endif
  } // omp parallel

  // Serialize the new interactions to the neighbor list
  for (int idx = 0; idx < nmsimples_; ++idx) {
    auto p1 = m_simples_[idx];
    nl_kmc_list *mlist = &mneighbors_[idx];
    for (auto nldx = mlist->begin(); nldx != mlist->end(); ++nldx) {
      auto p2 = m_simples_[nldx->idx_];
      interaction_t new_interaction;
      new_interaction.idx_ = (*oid_position_map_)[p1->GetOID()];
      new_interaction.jdx_ = (*oid_position_map_)[p2->GetOID()];
      new_interaction.type_ = type_;
      new_interaction.pot_ = pbase_;
      new_interaction.kmc_track_module_ = moduleid_;

      // KMC specifics
      if (type_ == ptype::kmc) {
        new_interaction.kmc_target_ = kmc_target_;
      }

      m_interactions_.push_back(new_interaction);
    }
  } //serialized add
}
