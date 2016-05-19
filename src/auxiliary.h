#ifndef _AUXILIARY_H_
#define _AUXILIARY_H_

#include <gsl/gsl_rng.h>
#include <vector>
#include <iostream>
#include "allocate.h"
#include "macros.h"
#include "timetester.h"
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "math.h"
#include "parameters.h"
#include <ios>
#include <fstream>
#include <sstream>

struct graph_struct {
  double r[3];
  double u[3];
  double length;
  double diameter;
};

struct rng_properties { 
  gsl_rng *r;
  const gsl_rng_type *T;
  void init(long seed) {
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
  }
  void clear() {
    gsl_rng_free(r);
  }
};

struct space_struct {
  int n_dim;
  int n_periodic;
  bool bud;
  double radius;
  double bud_radius;
  double bud_height;
  double **unit_cell;
  double **unit_cell_inv;
  double **a;
  double **b;
  double *a_perp;
  std::string type;
};

struct interaction {
  double *pot(double *r); // interaction potential (returns force)
  double *dr; // distance to object
  double rad; // interaction radius of object
};

double cpu();
void generate_random_unit_vector(int n_dim, double *vect, gsl_rng *r);
void rotate_orientation_vector(int n_dim, double *vect1, double *vect2);
double dot_product(int n_dim, double *a, double *b);
double determinant(int n, double **mat);
//double *separation_vector(int n_dim, int n_periodic, double *r1, double *s1, double *r2, double *s2, double **unit_cell);
void separation_vector(int n_dim, int n_periodic, double const * const r1, double const * const s1, double const * const r2, double const * const s2, double **unit_cell, double *dr);
void cross_product(double *a, double *b, double *c, int n_dim);
void normalize_vector(double *a, int n_dim);
void error_exit(const char *error_msg, ...);
void warning(const char *warning_msg);
void tridiagonal_solver(double *a, double *b, double *c, double *d, int n);
void rotate_3d_vector(double theta, double *a, double *b);
void invert_sym_2d_matrix(double **a, double **b); 
void invert_sym_3d_matrix(double **a, double **b); 
void periodic_boundary_conditions(int n_periodic, double **h, double **h_inv,
                                         double *r, double *s);

#endif // _AUXILIARY_H_

