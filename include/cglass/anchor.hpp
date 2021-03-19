#ifndef _CGLASS_ANCHOR_H_
#define _CGLASS_ANCHOR_H_

#include "filament.hpp"
#include "neighbor_list.hpp"

struct bind_params {
  double k_on_s;
  double k_on_d;
  density_type dens_type;
  double bind_site_density;
  bool single_occupancy;
};

/* Class for bound crosslink heads (called anchors). Tracks and updates its
   absolute position in space and relative position to its bound object. */
class Anchor : public Object {
 private:
  bool bound_;
  bool static_flag_;
  bool active_;
  bool plus_end_pausing_;
  bool minus_end_pausing_;
  crosslink_parameters *sparams_;
  int step_direction_;

  double rod_length_;
  double bond_lambda_;
  double mesh_length_;
  double mesh_lambda_;
  double max_velocity_s_;
  double max_velocity_d_;
  double diffusion_s_;
  double diffusion_d_;
  double kick_amp_s_;
  double kick_amp_d_;
  double k_on_s_;
  double k_on_d_;
  double k_off_s_;
  double k_off_d_;
  double polar_affinity_;
  double f_stall_;
  double force_dep_vel_flag_;
  bool use_bind_file_;
  std::map<std::string, bind_params> *bind_param_map_;

  double input_tol = 1e-8; // Tolerance for comparing inputs to 0

  bind_state state_;

  NeighborList neighbors_;

  Rod *rod_ = nullptr;
  Sphere *sphere_ = nullptr;
  Composite *comp_ = nullptr;

  // Retain these for walking behaviors
  Bond *bond_ = nullptr;
  Mesh *mesh_ = nullptr;

  // Richelle- completely replace obj_area with bind_rate.
  double *obj_area_ = nullptr;
  double *bind_rate_ = nullptr;

  int mesh_n_bonds_;

  // Helper functions
  void UpdateAnchorPositionToRod();
  void Diffuse();
  void Walk();
  bool CheckMesh();
  bool CalcRodLambda();

 public:
  Anchor(unsigned long seed);
  void Init(crosslink_parameters *sparams, int index, std::map<std::string, bind_params> *bind_param_map);
  bool IsBound();
  void UpdatePosition();
  void Activate();
  void Deactivate();
  void ApplyAnchorForces();
  void UpdateAnchorPositionToMesh();
  void SetDiffusion();
  void AttachObjRandom(Object *o);
  void AttachObjLambda(Object *o, double lambda);
  void AttachObjCenter(Object *o);
  void AttachObjMeshLambda(Object *o, double mesh_lambda);
  void AttachObjMeshCenter(Object *o);
  void CalculatePolarAffinity(std::vector<double> &doubly_binding_rates);
  void SetRodLambda(double l);
  void SetMeshLambda(double ml);
  void SetBound();
  void Unbind();
  int const GetBoundOID();
  void Draw(std::vector<graph_struct *> &graph_array);
  void AddNeighbor(Object *neighbor);
  void ClearNeighbors();
  const Object *const *GetNeighborListMem();
  const std::vector<Rod*>& GetNeighborListMemRods();
  const std::vector<Sphere*>& GetNeighborListMemSpheres();
  void WriteSpec(std::fstream &ospec);
  void ReadSpec(std::fstream &ispec);
  void BindToPosition(double *bind_pos);
  void SetStatic(bool static_flag);
  void SetState(bind_state state);

  double const GetMeshLambda();
  double const GetBondLambda();
  Object *GetNeighbor(int i_neighbor);
  Sphere *GetSphereNeighbor(int i_neighbor);
  Rod *GetRodNeighbor(int i_neighbor);
  const int GetNNeighbors() const;
  const int GetNNeighborsSphere() const;
  const int GetNNeighborsRod() const;
  const double GetOnRate() const;
  const double GetOnRateBindFile(std::string &name) const;
  const double GetOffRate() const;
  const double GetMaxVelocity() const;
  const double GetDiffusionConst() const;
  const double GetKickAmplitude() const;
  const double GetKonS() const;
  const double GetLinearBindSiteDensity() const;
  const double* const GetObjArea();
  void SetObjArea(double* obj_area);
  const double* const GetBindRate();
  void SetBindRate(double* bind_rate);
  bool InducesCatastrophe();
  bool AttachedToFilament();
  void InduceCatastrophe();

  // Convert binary data to text. Static to avoid needing to instantiate
  // species members in conversion mode.
  static void ConvertSpec(std::fstream &ispec, std::fstream &otext);
  static void WriteSpecTextHeader(std::fstream &otext);
};

#endif
