#ifndef _CGLASS_OBJECT_H_
#define _CGLASS_OBJECT_H_

#include "auxiliary.hpp"
#include "interaction.hpp"
#include "rng.hpp"
#include <mutex>

class Object {
private:
  int comp_id_;
  static int _next_oid_;
  static std::mutex _obj_mtx_;
  void InitOID();
  Object *comp_ptr_; // If part of a composite

protected:
  int oid_;
  static system_parameters *params_;
  static SpaceBase *space_;
  static int n_dim_;
  static double delta_;
  std::string name_;
  species_id sid_;
  obj_type type_ = obj_type::generic;
  comp_type comp_type_ = comp_type::generic;
  shape shape_ = shape::generic;
  graph_struct g_;
  RNG rng_;
  draw_type draw_;
  int n_contact_;
  double position_[3];
  double prev_position_[3];
  double prev_orientation_[3];
  double scaled_position_[3];
  double orientation_[3];
  double force_[3];
  double torque_[3];
  double dr_zero_[3];
  double color_;
  double diameter_;
  double length_;
  double p_energy_;
  double dr_tot_;
  double polar_order_;
  double contact_number_;
  bool interacting_;
  bool is_mesh_;
  bool is_comp_ = false;
  bool has_overlap_;
  bool fixed_ = false;
  int n_anchored_;
  bool interactor_update_;

  std::vector<Object *> interactors_;
  std::vector<object_interaction> ixs_;
  void UpdateKMC();

public:
  Object(unsigned long seed);
  virtual ~Object() = default;
  // kmcx parameter
  int gid;
  double length;
  double radius;
  double pos[3];
  double direction[3];

  // Static functions
  static void SetParams(system_parameters *params);
  static void SetSpace(SpaceBase *space);
  static void SetNDim(int n_dim);
  static void SetDelta(double delta);
  static const double GetDelta();
  static const int GetNextOID();
  static void SetNextOID(const int next_oid);
  // Trivial Set/Get functions
  void SetSID(species_id sid);
  void SetType(obj_type type);
  void SetPosition(const double *const new_pos);
  void SetScaledPosition(const double *const spos);
  void SetOrientation(const double *const u);
  void SetPrevPosition(const double *const ppos);
  void SetPrevOrientation(const double *const pu);
  void SetDiameter(double new_diameter);
  void SetLength(double new_length);
  virtual void AddForce(const double *const f);
  void SubForce(const double *const f);
  void SetForce(const double *const f);
  virtual void AddTorque(const double *const t);
  void SubTorque(const double *const t);
  void SetTorque(const double *const t);
  void AddPotential(const double p);
  void AddPolarOrder(const double po);
  void AddContactNumber(const double cn);
  void SetInteractor(bool ix);
  void IncrementNAnchored();
  void DecrementNAnchored();
  void ToggleIsMesh();
  virtual void CalcPolarOrder();
  virtual void ZeroPolarOrder();
  species_id const GetSID();
  obj_type const GetType();
  comp_type const GetCompType();
  shape const GetShape();
  std::string GetName();
  const int GetOID() const;
  const int GetCompID() const;
  const double *const GetPosition();
  const double *const GetPrevPosition();
  const double *const GetPrevOrientation();
  const double *const GetScaledPosition();
  const double *const GetOrientation();
  virtual void GetAvgPosition(double *ap);
  virtual void GetAvgOrientation(double *au);
  virtual void SetAvgPosition();
  const double GetDiameter();
  const double GetLength();
  const double *const GetForce();
  const double *const GetTorque();
  const double GetPotentialEnergy();
  const double GetPolarOrder();
  const double GetContactNumber();
  const bool IsInteractor();
  const bool IsMesh();
  const bool IsFixed();
  const int GetNAnchored();
  const bool CheckInteractorUpdate();
  void HasOverlap(bool overlap);
  void SetOID(int oid);
  void SetCompID(int cid);
  void SetCompPtr(Object* comp);

  // Virtual functions
  virtual void Init(species_base_parameters *sparams) {}
  virtual void InsertRandom(double buffer = -1);
  virtual void InsertRandomOriented(const double *const u);
  virtual void InsertAt(const double *const new_pos, const double *const u);
  virtual void ZeroForce();
  virtual void UpdatePeriodic();
  virtual void UpdatePosition() {}
  virtual void ResetPreviousPosition();
  virtual void Draw(std::vector<graph_struct *> &graph_array);
  virtual void SetColor(const double c, draw_type dtype);
  virtual void ScalePosition();
  virtual int GetCount();
  virtual void GetInteractors(std::vector<Object *> &ix);
  virtual const double *const GetInteractorPosition();
  virtual const double *const GetInteractorPrevPosition();
  virtual const double *const GetInteractorScaledPosition();
  virtual const double *const GetInteractorOrientation();
  virtual const double GetInteractorDiameter();
  virtual const double GetInteractorLength();
  virtual const double GetVolume();
  virtual const double GetArea();
  virtual void UpdateDrTot();
  virtual const double GetDrTot();
  virtual void ZeroDrTot();
  virtual bool HasNeighbor(int other_id);
  virtual void GiveInteraction(object_interaction ix);
  virtual void ApplyInteractions();
  virtual void FlagDuplicateInteractions();
  virtual void GetInteractions(std::vector<object_interaction> &ixs);
  virtual void CalcPCPosition(double s, double* pos);
  virtual void ClearInteractions();
  virtual void Cleanup();

  // I/O functions
  virtual void Report();
  virtual void WritePosit(std::fstream &oposit);
  virtual void ReadPosit(std::fstream &iposit);
  virtual void WriteSpec(std::fstream &ospec);
  virtual void ReadSpec(std::fstream &ispec);
  virtual void ReadPositFromSpec(std::fstream &ispec);
  virtual void WriteCheckpoint(std::fstream &ocheck);
  virtual void WriteCheckpointHeader(std::fstream &ocheck);
  virtual void ReadCheckpoint(std::fstream &icheck);
  virtual void ReadCheckpointHeader(std::fstream &icheck);
  virtual Object *GetCompPtr() { return comp_ptr_; }

  // Convert binary data to text. Static to avoid needing to istantiate
  // species members.
  static void ConvertPosit(std::fstream &iposit, std::fstream &otext);
  static void ConvertSpec(std::fstream &ispec, std::fstream &otext);
  static void WritePositTextHeader(std::fstream &otext);
  static void WriteSpecTextHeader(std::fstream &otext);
};

// void MinimumDistance(Object* o1, Object* o2, Interaction *ix, SpaceBase
// *space); void BoundaryConditions(Object * o1, SpaceBase *space);

#endif // _CGLASS_OBJECT_H_
