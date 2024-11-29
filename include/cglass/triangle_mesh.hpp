#ifndef _CGLASS_TRIANGLE_MESH_
#define _CGLASS_TRIANGLE_MESH_

// #include "common_libs.hpp"
// #include "definitions.hpp"
#include "meshbrane/half_edge_primitives.hpp"
#include "meshbrane/matrix_mesh.hpp"
#include "meshbrane/meshbrane_data_types.hpp"
#include "meshbrane/simple_generator.hpp"
#include "minimum_distance.hpp"
#include "rng.hpp"
#include "site.hpp"
#include <memory> // std::shared_ptr

namespace mbrn = meshbrane;
namespace hedge = meshbrane::half_edge;
// TODO add param storage and sync site seeds

struct Triangle;
struct Edge;
struct Vertex;
struct HalfEdge;
using HalfEdgePtr = std::shared_ptr<HalfEdge>;
using TrianglePtr = std::shared_ptr<Triangle>;
using EdgePtr = std::shared_ptr<Edge>;
using VertexPtr = std::shared_ptr<Vertex>;

// Vertex //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Vertex : public meshbrane::MeshBraneObject, public Site {
  ////////////////////////////
  // Fundamental attributes //
  ////////////////////////////
  // /**
  //  * @brief Index in master vrts_ list
  //  */
  // size_t index_{0};
  /**
   * @brief Pointer to an outgoing half-edge
   * 
   */
  HalfEdgePtr h_{nullptr};
  /**
   * @brief Seed for random number generation
   */
  int seed{0};
  /**
   * @brief Position in 3D space
   */
  double pos_[3];

  //////////////////////
  // Precomputed data //
  //////////////////////
  /**
   * @brief Number of triangles incident to this vertex
   */
  int n_tris_ = 0;
  /**
   * @brief Triangles incident to this vertex
   */
  std::vector<Triangle *> tris_;
  /**
   * @brief Number of vertices adjacent to this vertex
   */
  int n_neighbs_ = 0;
  /**
   * @brief Vertices adjacent to this vertex
   */
  std::vector<Vertex *> neighbs_;
  /**
   * @brief Number of edges incident to this vertex
   */
  int n_edges_ = 0; // should be replaced by n_neigbs when all is done
  /**
   * @brief Edges incident to this vertex
   */
  std::vector<Edge *> edges_;

  ////////////////////
  // Initialization //
  ////////////////////
  /**
   * @brief Construct a new Vertex object
   */
  Vertex() : Site(seed) {
    pos_[0] = pos_[1] = pos_[2] = 0;
    Site::SetPositionXYZ(0, 0, 0);
  }
  /**
   * @brief Construct a new Vertex object from x, y, and z coordinates
   */
  Vertex(double x, double y, double z) : Site(seed) {
    pos_[0] = x;
    pos_[1] = y;
    pos_[2] = z;
    Site::SetPositionXYZ(x, y, z);
  }
  /**
   * @brief Construct a new Vertex object from a `meshbrane::Coords3d` object
   */
  Vertex(mbrn::Coords3d xyz_coord) : Site(seed) {
    pos_[0] = xyz_coord(0);
    pos_[1] = xyz_coord(1);
    pos_[2] = xyz_coord(2);
    Site::SetPositionXYZ(xyz_coord(0), xyz_coord(1), xyz_coord(2));
  }

  ///////////////////////
  // Operators/Methods //
  ///////////////////////
  friend bool operator==(const Vertex &lhs, const Vertex &rhs) {
    return (lhs.pos_[0] == rhs.pos_[0] and lhs.pos_[1] == rhs.pos_[1] and
            lhs.pos_[2] == rhs.pos_[2]);
  }
  friend bool operator!=(const Vertex &lhs, const Vertex &rhs) {
    return !(lhs == rhs);
  }
  /**
   * @brief Set the position of the vertex
   */
  void SetPos(const double *const new_pos) {
    pos_[0] = pos[0] = position_[0] = new_pos[0];
    pos_[1] = pos[1] = position_[1] = new_pos[1];
    pos_[2] = pos[2] = position_[2] = new_pos[2];
  }
  /////////////////////
  // Getters/Setters //
  /////////////////////
  /**
   * @brief Get pointer to and outgoing half-edge
   */
  HalfEdgePtr h_out() { return h_; }
};

// Edge //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Edge : public meshbrane::MeshBraneObject {
  ////////////////////////////
  // Fundamental attributes //
  ////////////////////////////
  /**
   * @brief Vertices at each end of the edge
   */
  Vertex *vrts_[2]; // each endpoint
  /**
   * @brief Pointer to a (anti)parallel half-edge
   */
  HalfEdgePtr h_{nullptr};
  /**
   * @brief Flag to indicate if the edge was already flipped
   */
  bool just_flipped{false};

  //////////////////////
  // Precomputed data //
  //////////////////////
  /**
   * @brief Adjacent triangles
   */
  Triangle *tris_[2]; // each adjacent triangle
  /**
   * @brief Length of the edge
   */
  double length_{0.0};
  /**
   * @brief Vector pointing from v0 to v1
   */
  double vector_[3];
  /**
   * @brief Update the edge length `length_` and edge vector `vector_`
   * 
   */
  void Update() {
    length_ = 0.0;
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      vector_[i_dim] = vrts_[1]->pos_[i_dim] - vrts_[0]->pos_[i_dim];
      length_ += SQR(vector_[i_dim]);
    }
    length_ = sqrt(length_);
  }

  ////////////////////
  // Initialization //
  ////////////////////
  Edge() {}
  /**
   * @brief Construct a new Edge object from two vertices
   */
  Edge(Vertex *v0, Vertex *v1) {
    vrts_[0] = v0;
    vrts_[1] = v1;
  }

  ////////////////
  // Predicates //
  ////////////////
  friend bool operator==(const Edge &lhs, const Edge &rhs) {
    return ((lhs.vrts_[0] == rhs.vrts_[0] and lhs.vrts_[1] == rhs.vrts_[1]) or
            (lhs.vrts_[0] == rhs.vrts_[1] and lhs.vrts_[1] == rhs.vrts_[0]));
  }
  friend bool operator!=(const Edge &lhs, const Edge &rhs) {
    return !(lhs == rhs);
  }
  /**
   * @brief Check the edge contains a vertex
   */
  bool Contains(Vertex *vrt) { return vrts_[0] == vrt or vrts_[1] == vrt; }

  /////////////////////
  // Getters/Setters //
  /////////////////////
  /**
   * @brief Get the vertex at other end of the edge from a given vertex
   */
  Vertex *GetOtherEnd(Vertex *vrt) {
    return vrt == vrts_[0] ? vrts_[1] : vrts_[0];
  }
  /**
   * @brief Get Triangle on the other side of the edge from a given triangle
   */
  Triangle *GetOtherTriangle(Triangle *tri) {
    return tri == tris_[0] ? tris_[1] : tris_[0];
  }
  /**
   * @brief Get pointer to a parallel half-edge
   */
  HalfEdgePtr h_parallel() { return h_; }
};

// Triangle //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Triangle : public meshbrane::MeshBraneObject {
  ////////////////////////////
  // Fundamental attributes //
  ////////////////////////////
  /**
   * @brief Vertices at the corners of the triangle
   */
  Vertex *vrts_[3]{{}};
  /**
   * @brief Edges of the triangle 
   */
  Edge *edges_[3]{{}};
  /**
   * @brief Neighboring triangles
   */
  Triangle *neighbs_[3]{{}};
  /**
   * @brief Whether the triangle is flipped?
   */
  bool flipped_{false};
  /**
   * @brief Pointer to a half-edge on the positively oriented boundary of this face
   * 
   */
  HalfEdgePtr h_{nullptr};
  //////////////////////
  // Precomputed data //
  //////////////////////
  /**
   * @brief Triangle area
   */
  double area_;
  /**
   * @brief Volume of tetrahedron formed by triangle and origin
   */
  double volume_;
  /**
   * @brief Unit normal vector to the triangle
   */
  double nhat_[3];
  /**
   * @brief Color of the triangle. RGB values in [0, 255]. Doesn't seem to be used?
   */
  double color_[3];
  // euler angles
  /**
   * @brief cos(Euler angle) about ?-axis
   */
  double cosGamma_;
  /**
   * @brief sin(Euler angle) about ?-axis
   */
  double sinGamma_;
  /**
   * @brief cos(Euler angle) about ?-axis
   */
  double cosBeta_;
  /**
   * @brief sin(Euler angle) about ?-axis
   */
  double sinBeta_;
  /**
   * @brief ???
   */
  double Zrot_;
  /**
   * @brief ???
   */
  double XYrot_[2][3]; // [dim][i_vrt]

  void Update(double origin[]) {
    UpdateArea();
    UpdateVolume(origin);
  }
  void UpdateArea() {
    double a{edges_[0]->length_};
    double b{edges_[1]->length_};
    double c{edges_[2]->length_};
    double s{0.5 * (a + b + c)};
    area_ = sqrt(s * (s - a) * (s - b) * (s - c));
  }
  void UpdateVolume(double origin[]) {
    // update nhat
    // update volume
    double A[3];
    double B[3];
    double C[3];
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      A[i_dim] = vrts_[0]->pos_[i_dim] - origin[i_dim];
      B[i_dim] = vrts_[1]->pos_[i_dim] - origin[i_dim];
      C[i_dim] = vrts_[2]->pos_[i_dim] - origin[i_dim];
    }
    double BxC[3];
    cross_product(B, C, BxC, 3);
    volume_ = std::fabs(dot_product(3, A, BxC) / 6.0);
  }

  ////////////////////
  // Initialization //
  ////////////////////

  Triangle *neighbs_[3]{{}};

  Triangle(Vertex *v1, Vertex *v2, Vertex *v3) {
    vrts_[0] = v1;
    vrts_[1] = v2;
    vrts_[2] = v3;
    color_[0] = rand() % 255;
    color_[1] = rand() % 255;
    color_[2] = rand() % 255;
  }
  ////////////////
  // Predicates //
  ////////////////
  bool Contains(Vertex *vrt) {
    return vrts_[0] == vrt or vrts_[1] == vrt or vrts_[2] == vrt;
  }

  /////////////////////
  // Getters/Setters //
  /////////////////////
  double GetCenterPos(int i_dim) {
    return (vrts_[0]->pos_[i_dim] + vrts_[1]->pos_[i_dim] +
            vrts_[2]->pos_[i_dim]) /
           3.0;
  }
  Vertex *GetOtherVertex(Vertex *vrt1, Vertex *vrt2) {
    if (vrts_[0] != vrt1 and vrts_[0] != vrt2) {
      return vrts_[0];
    } else if (vrts_[1] != vrt1 and vrts_[1] != vrt2) {
      return vrts_[1];
    } else if (vrts_[2] != vrt1 and vrts_[2] != vrt2) {
      return vrts_[2];
    } else {
      printf("Error finding other vertex in triangle %zu\n", index_);
      exit(1);
    }
  }
  Edge *GetEdge(Vertex *vrt1, Vertex *vrt2) {
    for (auto &&edge : edges_) {
      if (edge->Contains(vrt1) and edge->Contains(vrt2)) {
        return edge;
      }
    }
    printf("Error finding edge in triangle %zu\n", index_);
    exit(1);
  }
};

// HalfEdge //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct HalfEdge {
  size_t index_{0}; // index in master half_edges_ list
};

class TriMesh : public mbrn::MatrixMesh {

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

  double l_avg_{0.0};                 // average edge length
  double node_drag_coefficient_{0.0}; // a

  // params for radial force
  double tether_stiffness_{0.0};
  double tether_attractive_singularity_{0.0};
  double tether_repulsive_singularity_{0.0};
  double tether_attractive_onset_{0.0};
  double tether_repulsive_onset_{0.0};
  // params for bending force
  // double bending_modulus_{0.0};
  // params for area conservation force
  double area_reg_stiffness_{0.0};
  double A_prime_{0.0};
  // params for volume conservation force
  double volume_reg_stiffness_{0.0};
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
  std::vector<HalfEdge> half_edges_;

  std::vector<graph_struct> f_mem_;
  graph_struct o_;

private:
  /**
 * @brief Set the Parameters object
 * 
 */
  void SetParameters();
  void MakeIcosphere();
  void MakeIcosahedron();
  void DivideFaces();
  void ProjectToUnitSphere();
  void InitializeMesh();
  void InitializeMeshBrane();
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

  ////////////////////////////////////////
  // WBL /////////////////////////////////
  ////////////////////////////////////////
public:
  //////////
  // Data //
  //////////
  std::string ply_path{"none"};

  double preferred_area_{1.0};
  double preferred_volume_{0.09403159725796};
  double spontaneous_curvature_{0.0};
  double bending_modulus_{0.0};
  double splay_modulus_{0.0};
  // double volume_reg_stiffness_{15.2};
  // double area_reg_stiffness_{3.0};
  double dimensionless_tether_repulsive_onset_{0.8};
  double dimensionless_tether_repulsive_singularity_{0.25};
  double dimensionless_tether_attractive_onset_{1.2};
  double dimensionless_tether_attractive_singularity_{2.5};
  // double tether_stiffness_{80.5};
  // double tether_repulsive_onset_{0.8};
  // double tether_repulsive_singularity_{0.25};
  // double tether_attractive_onset_{1.2};
  // double tether_attractive_singularity_{2.5};
  // double node_drag_coefficient_{0.03};
  double timestep_{1e-5};
  int flip_sweeps_per_step_{1};
  double flipping_probability_{0.3};

  //////////////////////////////////////
  // Constructors and related methods //
  //////////////////////////////////////
  /**
   * @brief Copy constructor
   */
  TriMesh(const TriMesh &other) : mbrn::MatrixMesh(other) {};
  void LoadPly();
  /**
   * @brief Refresh vertex, edge, and face lists from the matrix mesh data
   */
  void RefreshVEF_from_mats();

  void RefreshEdgeParams();

  /////////////////
  // Generators //
  ////////////////
};

#endif