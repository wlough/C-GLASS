#include "site.hpp"
#include <meshbrane/half_edge_mesh.hpp>
#include <meshbrane/meshbrane_data_types.hpp>
// #include <meshbrane/vertex_base.hpp>
#include <Eigen/Dense>

struct Vertex : public Site, public meshbrane::HalfEdgeVertex {
  size_t i_{0}; // index in master vertices_ list
  int seed{0};

  Vertex() : Site(seed), meshbrane::HalfEdgeVertex() {
    Site::SetPositionXYZ(0, 0, 0);
  }
  Vertex(double x, double y, double z) : Site(seed) {
    xyz_coord_ = meshbrane::Coords3d(x, y, z);
    Site::SetPositionXYZ(x, y, z);
  }

  void SetPos(const double *const new_pos) {
    xyz_coord_ = meshbrane::Coords3d(new_pos[0], new_pos[1], new_pos[2]);
    position_[0] = new_pos[0];
    position_[1] = new_pos[1];
    position_[2] = new_pos[2];
  }
};

using CMesh = meshbrane::HalfEdgeMesh<Vertex>;
class Brane : public CMesh {
public:
  Brane() : CMesh() {}
};