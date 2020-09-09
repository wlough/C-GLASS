#ifndef _CGLASS_DEFINITIONS_H_
#define _CGLASS_DEFINITIONS_H_

#include "macros.hpp"
#include "logger.hpp"
#include <math.h>

#define BETTER_ENUMS_DEFAULT_CONSTRUCTOR(Enum) \
 public:                                       \
  Enum() = default;
#include "enum.hpp"

#if defined(_OPENMP)
#define ENABLE_OPENMP
#endif

BETTER_ENUM(species_id, unsigned char, br_bead, filament, rigid_filament,
            spherocylinder, spindle, crosslink, none);
BETTER_ENUM(draw_type, unsigned char, fixed, orientation, bw, none);
BETTER_ENUM(potential_type, unsigned char, none, wca, soft, lj);
BETTER_ENUM(boundary_type, unsigned char, none = 0, box = 1, sphere = 2,
            budding = 3, wall = 4);
BETTER_ENUM(poly_state, unsigned char, grow, shrink, pause);  // make these 0,
                                                              // 1, 2 explicitly
BETTER_ENUM(bind_state, unsigned char, unbound, singly, doubly);
BETTER_ENUM(obj_type, unsigned char, generic, bond, site);

struct space_struct {
  int n_dim;
  int n_periodic;
  bool bud;
  double radius;
  double pressure_tensor[9];  // pressure tensor
  double pressure;            // isometric pressure
  double volume;
  double bud_radius;
  double bud_height;
  double bud_neck_radius;
  double bud_neck_height;
  double *unit_cell;
  double *unit_cell_inv;  // inverse unit cell
  double *a;              // direct lattice vector
  double *b;              // reciprocal lattice vector
  double *a_perp;  // perpendicular distance between opposite unit cell faces
  double *mu;      // scaling matrix for constant pressure
  int n_bound;     // number of bound motors
  double concentration;  // C of motors
  boundary_type type;
  
  // TO-DO: address Space vs. space question- could we turn into just a class?
  double BoundaryArea() const {
    switch (type) {
      case +boundary_type::none: {
        Logger::Error("Boundary area calculation requires a boundary.");
        break;
      }
      case +boundary_type::sphere: {
        if (n_dim == 2) {
          return 2.0 * M_PI * radius;
        } else {
          return 4 * M_PI * SQR(radius);
        }
      } 
      case +boundary_type::box: {
        if (n_dim == 2) {
          return 8.0 * radius;
        } else {
          return 24.0 * SQR(radius);
        }
      }
      case +boundary_type::budding: {
        double R = radius;
        double r = bud_radius;
        double d = bud_height;
        if (n_dim == 2) {
          // arc length is 2*r*(pi-theta) where theta is intersect angle
          return 2.0 * (r * (M_PI - acos((SQR(d) + SQR(r) - SQR(R)) / (2.0 * d * r))) 
                 + R * (M_PI - acos((SQR(d) - SQR(r) + SQR(R)) / (2.0 * d * R))));
        } else {
          // segment area is 2*pi*r^2*(1+cos(theta)) where theta is intersect angle
          return 2.0 * M_PI * (SQR(r) * (1 + ((SQR(d) + SQR(r) - SQR(R)) / (2.0 * d * r))) 
                 + SQR(R) * (1 + ((SQR(d) - SQR(r) + SQR(R)) / (2.0 * d * R))));
        }
      }
      case +boundary_type::wall: {
        Logger::Error("Cannot calculate area for wall boundary.");
        break;
      }
      default: 
        Logger::Error("Boundary type not recognized in space struct.");
    }
    return -1;
  }
};

struct graph_struct {
  double r[3];
  double u[3];
  double length;
  double diameter;
  double color;
  draw_type draw;
};

/* For unit testing */

// namespace unit_test {
// class Tester;
//}

#ifdef TESTS
#define UNIT_TEST       \
  template <typename T> \
  friend class UnitTest;
template <typename T>
class UnitTest {};
#else
#define UNIT_TEST
#endif

#endif
