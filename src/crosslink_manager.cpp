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
  //for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
  
  //  (*it)->CheckForCross();
  //}
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
}

//void CrosslinkManager::CheckForCross() {
//	int size=0;
//for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
//	size=((*it)->GetSize());
//	Logger::Warning("size= %i", size);
//	for (auto i = 0; i != size; ++i) {i
//		if (*it)->GetSize()
//		Logger::Warning("i=%i", i);
//	}
//}
  
//}

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
  }
}

void CrosslinkManager::InsertAttachedCrosslinks(std::vector<Object *> vO_, std::vector<Object *> vT_) {
  // Need to do this for GetRandomObject to work with spheres
  UpdateObjsSize();
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->UpdateBindRate();
    (*it)->InsertAttachedCrosslinksSpecies(vO_, vT_);
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
