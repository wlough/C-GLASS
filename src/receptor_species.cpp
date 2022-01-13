#include "cglass/receptor_species.hpp"

ReceptorSpecies::ReceptorSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::receptor);
}

// Use default initialization, plus include a concentration of receptors
void ReceptorSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  concentration_ = sparams_.concentration;
  component_ = sparams_.component;
}

/* To generalize, could pass vector of species and search for species
 * name. Would need to develop local random position on species' surface
 * and calculate area for each component, and static cast to mesh. */
void ReceptorSpecies::SetPC(Cortex* cx, std::vector<SpeciesBase *> &species) {
  if (component_.compare("cortex") == 0) {
    pc_ = cx;
  } else {
    for (auto sp = species.begin(); sp != species.end(); ++sp) {
      if ((*sp)->GetSpeciesName().compare(component_) == 0) {
        std::cout << "Follow this command: " << (*sp)->GetSpeciesName();
        std::cout <<"\n";
        pc_ = (*sp)->GetPC();
        pc_species_ = *sp;
      }
    }
    if (!pc_) {
      Logger::Error("Receptor component name is not a species name.");
    }
  }
}

void ReceptorSpecies::Reserve() {
  // Concentration overrides number of species
  if (concentration_ > 0) {
    sparams_.num = (int)round(concentration_ * pc_->GetArea());
  }
  members_.reserve(sparams_.num);
  Logger::Debug("Reserving memory for %d members in Receptor Species",
                sparams_.num);
}

/* Add receptor as a member of Receptor species and as a site on
 * the cell component. */
void ReceptorSpecies::AddMember() {
  // Initialize pointcover locations/variables
  if ((pc_species_) && (members_.size() == 0)) {
    smax_ = pc_species_->GetSpecLength() / 2.0;
    n_members_pc_ = pc_species_->GetNMembers();
    spacing_ = 2 * n_members_pc_ * smax_ / sparams_.num;
    s0_ = spacing_ / 2.0 - smax_;
    i_ = 0;
    s_ = s0_;
  }
  Species::AddMember();
  pc_->AddSpherePtr(&(members_.back()));
  double pos[3] =  {0,0,0};
  double u[3] = {1,0,0};
  if (sparams_.insertion_type.compare("random") == 0) {
    if (!pc_species_) {
      rng_.RandomBoundaryCoordinate(space_, pos);
    } else {
      // TO-DO: make more efficient for non-rod objects
      s_ = 2 * smax_ * (rng_.RandomUniform()-0.5);
      i_ = rng_.RandomInt(n_members_pc_);
      pc_species_->CalcPCPosition(i_, s_, pos);
      members_.back().SetLocations(i_, s_);
      members_.back().SetPCObject(pc_species_->GetMember(i_));
      //PC object seperatly for sphere class
      members_.back().SetPCObjectForSphere(pc_species_->GetMember(i_));
      members_.back().SetPCSpecies(pc_species_);
    }
  } else if (sparams_.insertion_type.compare("grid") == 0) {
    if (!pc_species_) {
      Logger::Error("Grid insertion not set up with cortex receptor component.");
    } else {
      if (s_ > smax_) {
        i_++;
        if (i_ >= n_members_pc_) {
          Logger::Error("Receptor species tried to insert at invalid member index.");
        }
        s_ = s0_;
      }
      pc_species_->CalcPCPosition(i_, s_, pos);
      members_.back().SetLocations(i_, s_);
      members_.back().SetPCObject(pc_species_->GetMember(i_));
      members_.back().SetPCObjectForSphere(pc_species_->GetMember(i_));
      members_.back().SetStepSize(spacing_);
      Logger::Warning("Spacing is %f, step is %f", spacing_, members_.back().GetStepSize());
      members_.back().SetPCSpecies(pc_species_);
      s_ += spacing_;
    }
  } else if (sparams_.insertion_type.compare("custom") != 0) {
    Logger::Error("Insertion type not recognized in receptor_species.cpp.");
  }
  members_.back().InsertAt(pos, u);
  members_.back().SetCompPtr(pc_);
}

// Overload to recognize grid arrangement
void ReceptorSpecies::ArrangeMembers() {
  if (GetInsertionType().compare("custom") == 0) {
    CustomInsert();
  } else if (GetInsertionType().compare("grid") != 0) {
    Logger::Error("Arrangement not recognized in receptor_species.cpp.");
  }
}

//Each receptor has it's neighbors set, this way diffusing/walking
//motors know which receptor to move to
void ReceptorSpecies::SetAllNeighbors() {
  int i=0;
  //const double *const vec = ((*pc_species_).members_[0])->GetOrientation();
  //std::cout << "Still laive " << vec[1];
  //std::cout <<"\n";
  //Cycle through each receptor
  while(i != members_.size()) {
		//If receptor is first receptor on filament it only has one neighbor 
    //set other neighbor to null pointer
    if (i==0) {
      members_[i].SetNeighbors(nullptr, &members_[i+1]);
    }
		//If receptor is last receptor on filament it only has one neighbor 
    //set other neighbor to null pointer
    else if (i==members_.size()-1) {
      members_[i].SetNeighbors(&members_[i-1],nullptr);
    }
    //Set both neighors if not on end
    else {
      members_[i].SetNeighbors(&members_[i-1], &members_[i+1]);
    }
    i+=1;
  }
}
