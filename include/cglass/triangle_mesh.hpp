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

  int i_tri_ = 0;
  Triangle *tris_[6]; // triangles this vertex is a part of
  int i_neighb_ = 0;
  std::vector<Vertex *> neighbs_;
  std::vector<bool> neighb_int_;

  Vertex() : Site(seed) {
    pos_[0] = pos[1] = pos[2] = 0;
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
};

struct Triangle {
  Vertex *vrts_[3]; // verteces that compose this triangle
  Triangle *neighbs_[3];
  graph_struct g_;
  Triangle(Vertex *v1, Vertex *v2, Vertex *v3) {
    vrts_[0] = v1;
    vrts_[1] = v2;
    vrts_[2] = v3;
  }
};

class TriMesh {

private:
  RNG rng_;

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