#ifndef _CGLASS_TRIANGLE_MESH_
#define _CGLASS_TRIANGLE_MESH_

#include "common_libs.hpp"
#include "definitions.hpp"
// #include "site.hpp"

struct Vertex {
  double pos_[3];
  graph_struct g_;
  //   Vertex() : g_({0, 0, 0}, {0, 0, 0}, 0, 0, 0, draw_type::orientation) {
  Vertex() {
    pos_[0] = pos_[1] = pos_[2] = 0.0;
    g_ = {{0, 0, 0}, {0, 0, 0}, 0, 0, 1, draw_type::fixed};
  }
  Vertex(double x, double y, double z) {
    pos_[0] = x;
    pos_[1] = y;
    pos_[2] = z;
    g_ = {{0, 0, 0}, {0, 0, 0}, 0, 0, 0, draw_type::orientation};
  }
};

struct Triangle {
  Vertex *vrts_[3];
  graph_struct g_;
  Triangle(Vertex *v1, Vertex *v2, Vertex *v3) {
    vrts_[0] = v1;
    vrts_[1] = v2;
    vrts_[2] = v3;
  }
};

struct Edge {
  double pos_[3];
  double len_;
  Vertex *vrt_one, *vrt_two;
  graph_struct g_;
};

class TriMesh {

private:
  int numVrt_ = -1;
  int numTri_ = -1;

  std::vector<std::vector<int>> indVert_;

  std::vector<Vertex> verts_;
  //   std::vector<Edge> edges_;

public:
  std::vector<Triangle> tris_;

private:
  void MakeIcosahedron();
  void MakeIcosphere();
  void ProjectToUnitSphere();

public:
  TriMesh();
  void Draw(std::vector<graph_struct *> &graph_array);
};

#endif