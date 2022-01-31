#include "cglass/crosslink_manager.hpp"

void CrosslinkManager::Init(system_parameters *params, SpaceBase *space, 
                            Tracker *tracker, std::vector<Object *> *objs) {
  objs_ = objs;
  update_ = false;
  obj_size_ = 0.0;
  params_ = params;
  space_ = space;
  tracker_ = tracker;
}

void CrosslinkManager::InitSpecies(sid_label &slab, ParamsParser &parser,
                                   unsigned long seed) {
  if (xlink_species_.size() == 0) {
    xlink_species_.reserve(parser.GetNCrosslinkSpecies());
  }
  xlink_species_.push_back(new CrosslinkSpecies(seed));
  xlink_species_.back()->Init(slab.second, parser);
  if (xlink_species_.back()->GetNInsert() <= 0) {
    delete xlink_species_.back();
    xlink_species_.pop_back();
  } else {
    xlink_species_.back()->InitInteractionEnvironment(objs_, &obj_size_, tracker_,
                                                      &update_, &bound_curr_);
    rcutoff_ = xlink_species_.back()->GetRCutoff();
  }
}

/* Keep track of volume of objects in the system. Affects the
 * probability of a free crosslink binding to an object. */
void CrosslinkManager::UpdateObjsSize() {
  obj_size_ = 0;
  for (auto it = objs_->begin(); it != objs_->end(); ++it) {
    switch ((*it)->GetShape()) {
      case shape::rod:
        obj_size_ += (*it)->GetLength();
        break;
      case shape::sphere:
        if ((*it)->GetNAnchored() == 0) obj_size_ += (*it)->GetArea();
        break;
      default:
        break;
    }
  }
}

/* Whether to reinsert anchors into the interactors list */
bool CrosslinkManager::CheckUpdate() {
  if (update_) {
    update_ = false;
    return true;
  }
  return false;
}

/* Return singly-bound anchors, for finding neighbors to bind to */
void CrosslinkManager::GetInteractors(std::vector<Object *> &ixors) {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->GetInteractors(ixors);
  }
}

/* Returns all anchors, not just singly-bound anchors. Used for reassigning
   bound anchors to bonds upon a checkpoint reload */
void CrosslinkManager::GetAnchorInteractors(std::vector<Object *> &ixors) {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->GetAnchorInteractors(ixors);
  }
}

void CrosslinkManager::UpdateCrosslinks() {
  update_ = false;
  UpdateObjsSize();
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->UpdatePositions();
    (*it)->UpdateBindRate();
  }
  Knockout();
  //Use global check for cross so we don't need to check individual crosslinkers every time step
  if (global_check_for_cross == true) {
    CheckForCross();
    global_check_for_cross = false;
  }
}

//Check if any crosslinkers are crossing
void CrosslinkManager::CheckForCross() {
  Logger::Trace("Checking if crosslinkers are crossing");
  //Sorry this is a mess to look at, this is going over over crosslink in every species 
  //and comparing it to every other crosslink in every species to see if the crosslinks
  //are crossing
  bool same_microtubles;
  for (auto spec_one = xlink_species_.begin(); spec_one != xlink_species_.end(); ++spec_one) {
  int member_num_ = (*spec_one) -> GetNMembers(); 
    for (int i = 0; i <= (member_num_-1); i++) {
      Crosslink* link_one = (*spec_one) -> GetCrosslink(i);
      // Check is we need to check for crossing for this crosslinker 
      if (link_one -> IsDoubly() && link_one -> ReturnCheckForCross() == true) {
        for (auto spec_two = xlink_species_.begin(); spec_two != xlink_species_.end(); ++spec_two) {
          int member_num_two_ = (*spec_two) -> GetNMembers();  
          for (int j = 0; j <= (member_num_two_-1); j++) {
            Crosslink* link_two = (*spec_two) -> GetCrosslink(j);
              //If link two is double bound and the crosslinkers aren't the same crosslinker
            if(link_two -> IsDoubly() && link_one -> GetOID() != link_two -> GetOID()) {
              //Get how far the crosslinker anchors are along the filament
              std::vector<double> link_one_s = link_one -> GetAnchorS();
              std::vector<double> link_two_s = link_two -> GetAnchorS();
              std::vector<int> rec_ids_one = link_one -> GetReceptorIDs();
              std::vector<int> rec_ids_two = link_two -> GetReceptorIDs();
              //If same index deosn't refer to sites on the same filament switch indexs so we are comparing same microtubule
              if (rec_ids_one[0] != rec_ids_two[0]) {
                std::swap(link_two_s[0], link_two_s[1]);
                std::swap(rec_ids_two[0], rec_ids_two[1]);
              } 
              //If the sites are still on different filaments then crosslinkers aren't binding the same microtubule 
              if (rec_ids_one[0] == rec_ids_two[0]) {
                same_microtubules = true;
              }
              else {
              same_microtubules = false;
              }

              if (same_microtubules && ((link_one_s[0] > link_two_s[0] && link_one_s[1] < link_two_s[1]) 
                     || (link_one_s[0] < link_two_s[0] && link_one_s[1] > link_two_s[1]))) { 
                Logger::Trace("Crosslinkers are crossing, %f, %f, %f, %f", link_one_s[0], link_one_s[1], link_two_s[0], link_two_s[1]);
                link_one -> UnbindCrossing();
                break;
              } else {
                link_one -> SetCheckForCross();
              } 
            }
          }
        }    
      }
    }
  }
}



// Loop over all spheres bound in the last dt. If multiple xlinks want to bind
// to a site, roll based on their relative probabilities, and unbind all of the
// "losers".
void CrosslinkManager::Knockout() {
  for (auto it = bound_curr_.begin(); it != bound_curr_.end(); ++it) {
    // Add an anchor now that loop over xlinks is complete (collisions have
    // already happened in this dt).
    it->first->IncrementNAnchored();
    if ((it->second.first).size()>1) {
      Logger::Error("Anchor bound twice in Knockout");
    }
  }
  bound_curr_.clear();
}

void CrosslinkManager::InsertCrosslinks() {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->InsertCrosslinks();
    (*it)->SetCheckForCrossPointer(&global_check_for_cross);
  } 
}

void CrosslinkManager::InsertAttachedCrosslinks() {
  // Need to do this for GetRandomObject to work with spheres
  UpdateObjsSize();
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->UpdateBindRate();
    (*it)->InsertAttachedCrosslinksSpecies();
  }
}

void CrosslinkManager::Clear() {
  output_mgr_.Close();
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->CleanUp();
    delete (*it);
  }
  xlink_species_.clear();
}

void CrosslinkManager::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->Draw(graph_array);
  }
}

void CrosslinkManager::AddNeighborToAnchor(Object *anchor, Object *neighbor) {
  if (anchor->GetSID() != +species_id::crosslink) {
    Logger::Error(
        "AddNeighborToAnchor expected crosslink object, got generic object.");
  }
  Anchor *a = dynamic_cast<Anchor *>(anchor);
  a->AddNeighbor(neighbor);
}

void CrosslinkManager::InitOutputs(bool reading_inputs, run_options *run_opts) {
  output_mgr_.Init(params_, &xlink_species_, space_, reading_inputs, run_opts);
}

void CrosslinkManager::WriteOutputs() { output_mgr_.WriteOutputs(); }

void CrosslinkManager::ZeroDrTot() {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->ZeroDrTot();
  }
}

const int CrosslinkManager::GetDoublyBoundCrosslinkNumber() const {
  int num = 0;
  for (const auto &xl_spec : xlink_species_) {
    num += xl_spec->GetDoublyBoundCrosslinkNumber();
  }
  return num;
}

const double CrosslinkManager::GetDrMax() {
  double dr_max = 0;
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    double dr = (*it)->GetDrMax();
    if (dr > dr_max) {
      dr_max = dr;
    }
  }
  return dr_max;
}

void CrosslinkManager::LoadCrosslinksFromCheckpoints(
    std::string run_name, std::string checkpoint_run_name) {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->LoadFromCheckpoints(run_name, checkpoint_run_name);
  }
}

void CrosslinkManager::ReadInputs() { output_mgr_.ReadInputs(); }
void CrosslinkManager::Convert() { output_mgr_.Convert(); }
