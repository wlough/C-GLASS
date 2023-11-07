#ifndef _CGLASS_TRIANGLE_MESH_
#define _CGLASS_TRIANGLE_MESH_

// #include "common_libs.hpp"
// #include "definitions.hpp"
#include "rng.hpp"
#include "site.hpp"

// TODO add param storage and sync site seeds

struct Triangle;
struct Vertex : public Site {
  int seed{0};
  double pos_[3];

  int n_tris_ = 0;
  Triangle *tris_[10]; // triangles this vertex is a part of
  int n_neighbs_ = 0;
  std::vector<Vertex *> neighbs_;

  // radial tether force jawn
  double fmag_[3]; // per neighb
  double rhat_[3]; // per neighb
  // bend energy jawn -- summed over all dem neighbs
  double sum_lsqT_;         // scalar
  double sum_del_lsqT_[3];  // vec; scalar in each dim
  double sum_rT_[3];        // vec; scalr in each dim
  double sum_del_rT_[3][3]; // tensor; vec. in each dim
  // area conservation energy jawn

  Vertex() : Site(seed) {
    pos_[0] = pos_[1] = pos_[2] = 0;
    Site::SetPositionXYZ(0, 0, 0);
  }
  Vertex(double x, double y, double z) : Site(seed) {
    pos_[0] = x;
    pos_[1] = y;
    pos_[2] = z;
    Site::SetPositionXYZ(x, y, z);
  }
  friend bool operator==(const Vertex &lhs, const Vertex &rhs) {
    return (lhs.pos_[0] == rhs.pos_[0] and lhs.pos_[1] == rhs.pos_[1] and
            lhs.pos_[2] == rhs.pos_[2]);
  }
  friend bool operator!=(const Vertex &lhs, const Vertex &rhs) {
    return !(lhs == rhs);
  }
  void SetPos(const double *const new_pos) {
    pos_[0] = position_[0] = new_pos[0];
    pos_[1] = position_[1] = new_pos[1];
    pos_[2] = position_[2] = new_pos[2];
  }
  void ZeroSums() {
    sum_lsqT_ = 0.0;
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      sum_rT_[i_dim] = 0.0;
      sum_del_lsqT_[i_dim] = 0.0;
      for (int j_dim{0}; j_dim < 3; j_dim++) {
        sum_del_rT_[i_dim][j_dim] = 0.0;
      }
    }
  }
};

struct Triangle {
  double area_;
  double nhat_[3];

  Vertex *vrts_[3]; // verteces that compose this triangle
  Triangle *neighbs_[3];
  Triangle(Vertex *v1, Vertex *v2, Vertex *v3) {
    vrts_[0] = v1;
    vrts_[1] = v2;
    vrts_[2] = v3;
  }
};

class TriMesh {

private:
  double r_sys_{0.0};

  double l_avg_{0.0};

  // params for radial force
  double kappa_B_{0.0};
  double l_max_{0.0};
  double l_min_{0.0};
  double l_c0_{0.0};
  double l_c1_{0.0};
  // params for bending force
  double kappa_{0.0};
  // params for area force
  double kappa_l_{0.0};
  double A_prime_{0.0};
  RNG rng_; // SF TODO link with system RNG

public:
  std::vector<Vertex> vrts_;
  std::vector<Triangle> tris_;

private:
  void MakeIcosphere();
  void MakeIcosahedron();
  void DivideFaces();
  void ProjectToUnitSphere();
  void UpdateNeighbors();
  void UpdatePositions();

public:
  TriMesh();
  void Draw(std::vector<graph_struct *> &graph_array);
};

#endif