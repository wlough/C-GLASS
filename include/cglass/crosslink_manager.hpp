#ifndef _CGLASS_CROSSLINK_MANAGER_H_
#define _CGLASS_CROSSLINK_MANAGER_H_

#include "output_manager.hpp"
#include "rng.hpp"
#include "species_factory.hpp"
#include <numeric>

class CrosslinkOutputManager : public OutputManagerBase<CrosslinkSpecies> {
  /* Do not write thermo - base output manager will handle that */
  void WriteThermo() {}
  void ReadThermo() {}
  void InitThermo(std::string fname) {}
};

class CrosslinkManager {
private:
  system_parameters *params_;
  CrosslinkOutputManager output_mgr_;
  double obj_size_;
  double rcutoff_ = 0; // Cutoff for binding any crosslink and bond
  bool update_;
  std::vector<CrosslinkSpecies *> xlink_species_;
  std::vector<Object *> *objs_;
  std::map<Sphere *, std::pair<std::vector<double>,
                               std::vector<std::pair<Anchor *, std::string>>>>
      bound_curr_;
  std::map<std::string, CrosslinkSpecies *> species_map_;
  SpaceBase *space_;
  Tracker *tracker_ = nullptr;
  bool global_check_for_cross = false;
  bool same_microtubules = true;
  void KnockoutBind(Sphere *receptor, int winner);
  RNG *rng_;
  std::vector<std::pair<CrosslinkSpecies *, Sphere *>> cl_to_add_this_step_;
  void AddKnockoutCrosslinks();

public:
  void Init(system_parameters *params, SpaceBase *space, Tracker *tracker,
            std::vector<Object *> *objs, unsigned long seed);
  void GetInteractors(std::vector<Object *> &ixors);
  void UpdateCrosslinks();
  void UpdateObjsSize();
  bool CheckUpdate();
  void CheckForCross();
  void Clear();
  void Draw(std::vector<graph_struct *> &graph_array);
  void AddNeighborToAnchor(Object *anchor, Object *neighbor);
  void WriteOutputs();
  void Knockout();
  void InitOutputs(bool reading_inputs = false,
                   run_options *run_opts = nullptr);
  void GetAnchorInteractors(std::vector<Object *> &ixors);
  void InitSpecies(sid_label &slab, ParamsParser &parser, unsigned long seed);
  void LoadCrosslinksFromCheckpoints(std::string run_name,
                                     std::string checkpoint_run_name);
  const int GetDoublyBoundCrosslinkNumber() const;
  void ZeroDrTot();
  const double GetDrMax();
  void ReadInputs();
  void Convert();
  void InsertCrosslinks();
  void
  InsertAttachedCrosslinks(std::vector<std::vector<Object *>> receptor_list);
  const double GetRCutoff() const {
    Logger::Trace("Crosslink rcutoff is %2.2f", rcutoff_);
    return rcutoff_;
  }

  CrosslinkSpecies *GetSPBCrosslinkSpecies(std::string name) {
    auto map_entry{species_map_.find(name)};
    if (map_entry == species_map_.end()) {
      return nullptr;
    }
    return map_entry->second;
  }
};

#endif
