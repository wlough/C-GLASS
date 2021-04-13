#ifndef _CGLASS_CROSSLINK_H_
#define _CGLASS_CROSSLINK_H_

//#include "species.hpp"
#include "anchor.hpp"
#include "minimum_distance.hpp"
#include "tracker.hpp"
#include <KMC/kmc.hpp>
#include <KMC/kmc_choose.hpp>

// enum xstate { unbound, singly, doubly };

/* Class that represents a two-headed crosslink that can create tethers between
 * two objects in the simulation. Governs binding and unbinding of crosslink
 * heads, and tether forces between bound heads. */
class Crosslink : public Object {
private:
  crosslink_parameters *sparams_;
  draw_type draw_;
  bind_state state_;
  LookupTable *lut_;
  bool static_flag_ = false;
  bool asymmetric_spring_flag_ = false;
  double k_spring_;
  double k_spring_compress_;
  double k_align_;
  double rest_length_;
  double rcapture_;
  double bind_site_density_;
  double tether_force_;
  double e_dep_factor_;
  double fdep_length_;
  double polar_affinity_;
  bool use_bind_file_;
  int bound_anchor_ = 0; // Index of anchor that is bound if Singly
  std::map<Sphere *, std::pair<std::vector<double>, std::vector<Anchor*> > > *bound_curr_;
  std::vector<std::map<std::string, bind_params> > *bind_param_map_;
  double *bind_rate_;
  std::vector<Anchor> anchors_;
  void CalculateTetherForces();
  void CalculateBinding();
  void SinglyKMC();
  void DoublyKMC();
  void UpdateAnchorsToMesh();
  void UpdateAnchorPositions();
  void UpdateXlinkState();
  double *obj_area_ = nullptr;
  Tracker *tracker_ = nullptr;
public:
  Crosslink(unsigned long seed);
  void Init(crosslink_parameters *sparams);
  void InitInteractionEnvironment(LookupTable *lut, Tracker *tracker, 
                                  std::map<Sphere *, std::pair<std::vector<double>, 
                                  std::vector<Anchor*> > > *bound_curr);
  void SetBindParamMap(std::vector<std::map<std::string, bind_params> > 
                       *bind_param_map);
  void AttachObjRandom(Object *obj);
  void UpdateCrosslinkForces();
  void UpdateCrosslinkPositions();
  void GetAnchors(std::vector<Object *> &ixors);
  void GetInteractors(std::vector<Object *> &ixors);
  void Draw(std::vector<graph_struct *> &graph_array);
  void SetDoubly();
  void SetSingly(int bound_anchor);
  void SetUnbound();
  void SetAnchorStates();
  const bool IsDoubly() const;
  const bool IsUnbound() const;
  const bool IsSingly() const;
  void UpdatePosition();
  void WriteSpec(std::fstream &ospec);
  void WriteCheckpoint(std::fstream &ocheck);
  void ReadSpec(std::fstream &ispec);
  void ReadCheckpoint(std::fstream &icheck);
  void ClearNeighbors();
  void ZeroForce();
  void ApplyTetherForces();
  void ZeroDrTot();
  const double GetDrTot();
  void InsertAt(double const *const new_pos, double const *const u);
  const int GetNNeighbors() const;
  void SetObjArea(double *obj_area);
  void SetBindRate(double *bind_rate);
  void SetSpheresBoundCurr(double *obj_area);
  const double* const GetObjArea();
  const double *const GetPosition();
  const double *const GetOrientation();

  // Convert binary data to text. Static to avoid needing to istantiate
  // species members.
  static void ConvertSpec(std::fstream &ispec, std::fstream &otext);
  static void WriteSpecTextHeader(std::fstream &otext);
};

#endif
