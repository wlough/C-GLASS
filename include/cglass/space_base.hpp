#ifndef _CGLASS_SPACE_BASE_H_
#define _CGLASS_SPACE_BASE_H_

#include "auxiliary.hpp"
#include "macros.hpp"
#include "triangle_mesh.hpp"
#include <math.h>

class SpaceBase {
public:
  int n_dim;
  int n_periodic;
  bool bud;
  double radius;
  double pressure_tensor[9]; // pressure tensor
  double pressure;           // isometric pressure
  double volume;
  double bud_radius;
  double bud_height;
  double bud_neck_radius;
  double bud_neck_height;
  double *unit_cell;
  double *unit_cell_inv; // inverse unit cell
  double *a;             // direct lattice vector
  double *b;             // reciprocal lattice vector
  double *a_perp; // perpendicular distance between opposite unit cell faces
  double *mu;     // scaling matrix for constant pressure
  int n_bound;    // number of bound motors
  double concentration; // C of motors
  boundary_type type;
  double BoundaryArea() const;
  // triangle mesh data
  TriMesh mesh_;
};

#endif
