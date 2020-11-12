#ifndef _CGLASS_FILAMENT_FLOCKING_ANALYSIS_H_
#define _CGLASS_FILAMENT_FLOCKING_ANALYSIS_H_

#include "analysis.hpp"
#include "filament.hpp"
#include "minimum_distance.hpp"

class BindingAnalysis : public Analysis<Receptor, species_id::receptor> {
protected:
  void InitOutput() {
    SetAnalysisName("binding");
    Analysis::InitOutput();
  }
  void InitAnalysis() {
    output_ << "diam nrecept\n";
    auto it = members_->begin();
    double d = it->GetDiameter();
    double nr = members_->size();
    output_ << d << " " << nr << std::endl; 
    output_ << "time n_bound";
  }
  void RunAnalysis() {
    int n_bound = 0;
    for (auto it = members_->begin(); it != members_->end(); ++it) {
      n_bound+=it->GetNAnchored();
    }
    if (output_.is_open()) {
      output_ << time_ << " " << n_bound << std::endl;
    } else {
      Logger::Error("Problem opening file in FlockingAnalysis!");
    }
  }

  void EndAnalysis() { }
};

#endif // _CGLASS_FILAMENT_FLOCKING_ANALYSIS_H_
