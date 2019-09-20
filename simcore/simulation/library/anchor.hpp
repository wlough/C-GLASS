#ifndef _SIMCORE_ANCHOR_H_
#define _SIMCORE_ANCHOR_H_

#include "mesh.hpp"

class Anchor : public Object {
private:
  bool bound_;
  bool walker_;
  bool diffuse_;
  bool active_;
  bool end_pausing_;

  int step_direction_;

  double bond_length_;
  double bond_lambda_;
  double mesh_length_;
  double mesh_lambda_;
  double velocity_;
  double max_velocity_, diffusion_;
  double f_spring_max_;
  double k_off_;
  double f_stall_;
  double force_dep_vel_flag_;

  Bond *bond_;
  Mesh *mesh_;

  int bond_oid_;
  int mesh_n_bonds_;

  void UpdateAnchorPositionToBond();
  void Diffuse();
  void Walk();
  bool CheckMesh();

public:
  Anchor();
  void Init();
  bool IsBound();
  void UpdatePosition();
  void Activate();
  void Deactivate();
  void ApplyAnchorForces();
  void UpdateAnchorPositionToMesh();
  void SetDiffusion();
  void SetWalker(int dir, double walk_v);
  void AttachObjRandom(Object *o);
  void AttachObjLambda(Object *o, double lambda);
  double const GetMeshLambda();
  double const GetBondLambda();
  void SetBondLambda(double l);
  void SetBound();
  void Clear();
  int const GetBoundOID();
  void Draw(std::vector<graph_struct *> *graph_array);
};

#endif
