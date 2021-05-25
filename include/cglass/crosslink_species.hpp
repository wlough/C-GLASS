#ifndef _CGLASS_CROSSLINK_SPECIES_H_
#define _CGLASS_CROSSLINK_SPECIES_H_

#include "crosslink.hpp"
#include "species.hpp"
#include <KMC/kmc.hpp>

typedef std::vector<std::pair<std::vector<Crosslink>::iterator,
                              std::vector<Crosslink>::iterator>>
    xlink_chunk_vector;
typedef std::vector<Crosslink>::iterator xlink_iterator;

class CrosslinkSpecies : public Species<Crosslink, species_id::crosslink> {
private:
  bool *update_;
  bool midstep_ = true;
  std::string checkpoint_file_;
  double *obj_length_; // Total length of all the objects in the system
  double *obj_area_; // Total area of all the objects in the system
  double xlink_concentration_;
  int begin_with_bound_crosslinks_;
  double linear_bind_site_density_;
  double surface_bind_site_density_;
  bool infinite_reservoir_flag_;
  double k_on_;
  bool static_flag_;
  std::vector<std::string> bind_species_;
  Tracker *tracker_ = nullptr;
  LookupTable lut_;
  std::vector<Object *> *objs_;
  std::map<Sphere *, std::pair<std::vector<double>, std::vector<Anchor*> > > *bound_curr_;

  LUTFiller *MakeLUTFiller();
  void CalculateBindingFree();
  void BindCrosslink(shape sh);
  void UpdateBoundCrosslinks();
  void UpdateBoundCrosslinkForces();
  void UpdateBoundCrosslinkPositions();
  void ApplyCrosslinkTetherForces();
  Object *GetRandomObject(shape sh);

public:
  CrosslinkSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void LoadBindingSpecies();
  void InitInteractionEnvironment(std::vector<Object *> *objs, double *obj_len,
                                  double *obj_area, Tracker *tracker, bool *update,
                                  std::map<Sphere *, std::pair<std::vector<double>, std::vector<Anchor*> > > *bound_curr);
  void TestKMCStepSize();
  void GetInteractors(std::vector<Object *> &ixors);
  void UpdatePositions();
  void UpdateObjectArea();
  void CleanUp();
  void Draw(std::vector<graph_struct *> &graph_array);
  void BindCrosslinkObj(Object *obj);
  void AddNeighborToAnchor(Object *anchor, Object *neighbor);
  void AddMember();
  void InsertAttachedCrosslinksSpecies();
  void GetAnchorInteractors(std::vector<Object *> &ixors);
  void ReadSpecs();
  void InsertCrosslinks();
  const int GetDoublyBoundCrosslinkNumber() const;
  const double GetConcentration() const;
  const double GetRCutoff() const;
};

#endif