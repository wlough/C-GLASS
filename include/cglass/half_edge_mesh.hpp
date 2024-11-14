#pragma once

#include "meshbrane/meshbrane_data_types.hpp"
#include "minimum_distance.hpp"
#include "rng.hpp"
#include "site.hpp"

namespace halfedge {
struct HalfEdge;

struct Vertex : public Site {
  int i_{0}; // index in master vertices_ list
  int seed{0};
  meshbrane::Coords3d xyz_coord_;
  HalfEdge *h_out_{nullptr};

  Vertex() : Site(seed) {
    xyz_coord_.setZero();
    Site::SetPositionXYZ(0, 0, 0);
  }
  Vertex(double x, double y, double z) : Site(seed) {
    xyz_coord_[0] = x;
    xyz_coord_[1] = y;
    xyz_coord_[2] = z;
    Site::SetPositionXYZ(x, y, z);
  }
  friend bool operator==(const Vertex &lhs, const Vertex &rhs) {
    return (lhs.xyz_coord_[0] == rhs.xyz_coord_[0] and
            lhs.xyz_coord_[1] == rhs.xyz_coord_[1] and
            lhs.xyz_coord_[2] == rhs.xyz_coord_[2]);
  }
  friend bool operator!=(const Vertex &lhs, const Vertex &rhs) {
    return !(lhs == rhs);
  }
  void SetPos(const double *const new_pos) {
    xyz_coord_ = Eigen::Map<const Eigen::Matrix<double, 3, 1>>(new_pos);
    pos[0] = position_[0] = new_pos[0];
    pos[1] = position_[1] = new_pos[1];
    pos[2] = position_[2] = new_pos[2];
  }
};

struct Face {
  int i_{0}; // index in master faces_ list
  HalfEdge *h_right_{nullptr};

  Face();

  friend bool operator==(const Face &lhs, const Face &rhs) {
    return lhs.h_right_ == rhs.h_right_;
  }
  friend bool operator!=(const Face &lhs, const Face &rhs) {
    return !(lhs == rhs);
  }
};

struct Edge {
  int i_{0}; // index in master faces_ list
  HalfEdge *h_right_{nullptr};
  Edge();

  friend bool operator==(const Edge &lhs, const Edge &rhs) {
    return (lhs.h_right_ == rhs.h_right_ or
            lhs.h_right_->h_twin_ == rhs.h_right_);
  }
  friend bool operator!=(const Edge &lhs, const Edge &rhs) {
    return !(lhs == rhs);
  }
};

struct HalfEdge {
  int i_{0}; // index in master half_edges_ list
  Vertex *v_origin_{nullptr};
  HalfEdge *h_twin_{nullptr};
  HalfEdge *h_next_{nullptr};
  Face *f_left_{nullptr};

  HalfEdge() {}
  HalfEdge(Vertex *v) : v_origin_(v) {}
  HalfEdge(Vertex *v, HalfEdge *twin, HalfEdge *next, Face *f)
      : v_origin_(v), h_twin_(twin), h_next_(next), f_left_(f) {}
};

class HalfEdgeMesh {
private:
  static const size_t n_edges_min_{3};  // true for any connected graph
  static const size_t n_edges_max_{10}; // arbitrary choice

  bool do_not_pass_go_{false};
  int i_datapoint_{0};

  FILE *forces_{nullptr};    // average force from each type of potential
  FILE *vertices_{nullptr};  // position (3D per vrt per step)
  FILE *adjacency_{nullptr}; // adjacency matrix  (2D per vrt per step)
  system_parameters *params_{nullptr};

  double f_avgs_[4]; // indices 0-4: tether, bend, area, vol

  double l_avg_{0.0};
  double gamma_{0.0};

  // params for radial force
  double kappa_B_{0.0};
  double l_max_{0.0};
  double l_min_{0.0};
  double l_c0_{0.0};
  double l_c1_{0.0};
  // params for bending force
  double kappa_{0.0};
  // params for area conservation force
  double kappa_l_{0.0};
  double A_prime_{0.0};
  // params for volume conservation force
  double kappa_v_{0.0};
  double V_prime_{0.0};
  RNG *rng_; // SF TODO link with system RNG
  MinimumDistance mindist_;

public:
  double r_sys_{0.0};
  double centroid_[3];

  std::vector<Object *> boundary_neighbs_;

  std::vector<Vertex> vrts_;
  std::vector<Triangle> tris_;
  std::vector<Edge> edges_;

  std::vector<graph_struct> f_mem_;
  graph_struct o_;

private:
  void SetParameters();
  void MakeIcosphere();
  void MakeIcosahedron();
  void DivideFaces();
  void ProjectToUnitSphere();
  void InitializeMesh();
  void FlipEdges();
  void UpdateCentroid();
  void UpdateTriangles();
  void UpdateNeighbors();
  void UpdateMesh();
  void ApplyMembraneForces();
  void ApplyBoundaryForces();

public:
  TriMesh() {}
  void Init(system_parameters *params);
  void Draw(std::vector<graph_struct *> &graph_array);
  void UpdatePositions();
  void WriteOutputs();

  ////////////////////////////////////////////////////////////////////////////
  // WLOUGH
public:
  std::string ply_path{"none"};
  void load_ply();
  ////////////////////////////////////////////////////////////////////////////
};

} // namespace halfedge