// #ifndef _CGLASS_TRIANGLE_MESH_
// #define _CGLASS_TRIANGLE_MESH_
#pragma once
// #include "common_libs.hpp"
// #include "definitions.hpp"
#include "meshbrane/matrix_mesh.hpp"
#include "meshbrane/meshbrane_data_types.hpp"
#include "meshbrane/meshbrane_object.hpp"
#include "meshbrane/simple_generator.hpp"
#include "minimum_distance.hpp"
#include "rng.hpp"  //
#include "site.hpp" // Site
// #include <memory>   // std::shared_ptr
// #include <optional> // std::optional
// namespace half_edge {
namespace mbrn = meshbrane;
// TODO add param storage and sync site seeds

struct Edge;
struct HalfEdge;
struct Triangle;
struct Vertex;
using EdgePtr = std::shared_ptr<Edge>;
using HalfEdgePtr = std::shared_ptr<HalfEdge>;
using TrianglePtr = std::shared_ptr<Triangle>;
using VertexPtr = std::shared_ptr<Vertex>;
using HalfEdgeGenerator = meshbrane::utils::SimpleGenerator<HalfEdge *>;
using TriangleGenerator = meshbrane::utils::SimpleGenerator<Triangle *>;

// namespace mth {
// template <typename T>
// inline T cross(T u, T v) {
//   return T({u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
//             u[0] * v[1] - u[1] * v[0]});
// }
// } // namespace mth

// Vertex //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief A vertex in a meshed surface. See `
 * 
 */
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
  HalfEdge *h_{nullptr};
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

  /** 
   * @brief Update the list of neighboring vertices, edges, and triangles
   */
  void UpdateNeighbors();
  ////////////////////
  // Initialization //
  ////////////////////
  Vertex(size_t index) : mbrn::MeshBraneObject(index), Site(seed) {
    pos_[0] = pos_[1] = pos_[2] = 0;
    Site::SetPositionXYZ(0, 0, 0);
  }
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
   * @brief Construct a new Vertex object
   * @param index Index in master vrts_ list
   * @param xyz_coord Coordinates of the vertex
   * @param h Pointer to an outgoing half-edge
   */
  Vertex(size_t index, mbrn::Coords3d xyz_coord, HalfEdge *h)
      : mbrn::MeshBraneObject(index), Site(seed), h_(h) {
    pos_[0] = xyz_coord(0);
    pos_[1] = xyz_coord(1);
    pos_[2] = xyz_coord(2);
    Site::SetPositionXYZ(xyz_coord(0), xyz_coord(1), xyz_coord(2));
  }

  ///////////////
  // Operators //
  //////////////
  friend bool operator==(const Vertex &lhs, const Vertex &rhs) {
    return (lhs.pos_[0] == rhs.pos_[0] and lhs.pos_[1] == rhs.pos_[1] and
            lhs.pos_[2] == rhs.pos_[2]);
  }
  friend bool operator!=(const Vertex &lhs, const Vertex &rhs) {
    return !(lhs == rhs);
  }
  /////////////////////
  // Getters/Setters //
  /////////////////////
  /**
   * @brief Set the position of the vertex
   */
  void SetPos(const double *const new_pos) {
    pos_[0] = pos[0] = position_[0] = new_pos[0];
    pos_[1] = pos[1] = position_[1] = new_pos[1];
    pos_[2] = pos[2] = position_[2] = new_pos[2];
  }
  /**
   * @brief Get pointer to and outgoing half-edge
   */
  HalfEdge *h_out() { return h_; }
  void set_outgoing_half_edge(HalfEdge *h) { h_ = h; }

  ////////////////
  // Generators //
  ////////////////
  /**
   * @brief Generate half-edges which originate at this vertex in clockwise order
   */
  HalfEdgeGenerator generate_H_out_clockwise();
  /**
   * @brief Generate half-edges which originate at this vertex in clockwise order
   * @param h_start Starting half-edge
   */
  HalfEdgeGenerator generate_H_out_clockwise(HalfEdge *h_start);
  /**
   * @brief Generate faces incident to this vertex
   */
  TriangleGenerator generate_F_incident_clockwise();
};

// Edge //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Edge : public meshbrane::MeshBraneObject {
  ////////////////////////////
  // Fundamental attributes //
  ////////////////////////////
  /**
   * @brief Vertices at each end of the edge
   */
  Vertex *vrts_[2]{{}}; // each endpoint
  /**
   * @brief Pointer to a (anti)parallel half-edge
   */
  HalfEdge *h_{nullptr};
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
  Triangle *tris_[2]{{}}; // each adjacent triangle
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

  void UpdateNeighborTris();

  ////////////////////
  // Initialization //
  ////////////////////
  Edge() = default;
  /**
   * @brief Construct a new Edge object from two vertices
   */
  Edge(Vertex *v0, Vertex *v1) {
    vrts_[0] = v0;
    vrts_[1] = v1;
  }

  Edge(size_t index, Vertex *v0, Vertex *v1, HalfEdge *h)
      : mbrn::MeshBraneObject(index), h_(h) {
    vrts_[0] = v0;
    vrts_[1] = v1;
  }
  ////////////////
  // Predicates //
  ////////////////
  bool is_in_some_boundary() const;
  bool is_flippable() const;        // TODO
  bool is_locally_delaunay() const; // TODO

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
  HalfEdge *h_parallel() const { return h_; }

  /**
   * @brief Set the parallel half-edge
   */
  void set_parallel_half_edge(HalfEdge *h) { h_ = h; }
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
  HalfEdge *h_{nullptr};
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

  /**
   * @brief Update face area volume of tetrahedron formed by triangle and origin
   * 
   * @param origin 
   */
  void Update(double origin[]) {
    UpdateArea();
    UpdateVolume(origin);
  }
  /**
   * @brief Update stored triangle area
   * 
   * @param origin
   */
  void UpdateArea() {
    double a{edges_[0]->length_};
    double b{edges_[1]->length_};
    double c{edges_[2]->length_};
    double s{0.5 * (a + b + c)};
    area_ = sqrt(s * (s - a) * (s - b) * (s - c));
  }
  /**
   * @brief Update unit normal vector and volume of tetrahedron formed by triangle and origin
   * 
   * @param origin 
   */
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
  Triangle() = default;
  Triangle(size_t index) : mbrn::MeshBraneObject(index) {}
  /**
   * @brief Construct a new Triangle object from three vertices
   */
  Triangle(Vertex *v1, Vertex *v2, Vertex *v3) {
    vrts_[0] = v1;
    vrts_[1] = v2;
    vrts_[2] = v3;
    color_[0] = rand() % 255;
    color_[1] = rand() % 255;
    color_[2] = rand() % 255;
  }

  Triangle(size_t index, HalfEdge *h) : mbrn::MeshBraneObject(index), h_(h) {
    // color_[0] = rand() % 255;
    // color_[1] = rand() % 255;
    // color_[2] = rand() % 255;
    // h_ = std::make_shared<HalfEdge>(h);
  }
  Triangle(size_t index, Vertex *v1, Vertex *v2, Vertex *v3, HalfEdge *h)
      : mbrn::MeshBraneObject(index), h_(h) {
    vrts_[0] = v1;
    vrts_[1] = v2;
    vrts_[2] = v3;
    // color_[0] = rand() % 255;
    // color_[1] = rand() % 255;
    // color_[2] = rand() % 255;
  }

  // /**
  //  * @brief Copy constructor
  //  */
  // Triangle(const Triangle &other) : mbrn::MeshBraneObject(other) {};
  ////////////////
  // Predicates //
  ////////////////
  /**
   * @brief Check if the triangle contains a vertex
   */
  bool Contains(Vertex *vrt) {
    return vrts_[0] == vrt or vrts_[1] == vrt or vrts_[2] == vrt;
  }

  friend bool operator==(const Triangle &lhs, const Triangle &rhs) {
    return std::is_permutation(std::begin(lhs.vrts_), std::end(lhs.vrts_),
                               std::begin(rhs.vrts_));
  }

  friend bool operator!=(const Triangle &lhs, const Triangle &rhs) {
    return !(lhs == rhs);
  }

  /////////////////////
  // Getters/Setters //
  /////////////////////
  /**
   * @brief Get centroid of the triangle
   */
  double GetCenterPos(int i_dim) {
    return (vrts_[0]->pos_[i_dim] + vrts_[1]->pos_[i_dim] +
            vrts_[2]->pos_[i_dim]) /
           3.0;
  }
  /**
   * @brief Get the vertex at the other end of the edge from two vertices
   */
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
  /**
   * @brief Get the edge between two vertices
   */
  Edge *GetEdge(Vertex *vrt1, Vertex *vrt2) {
    for (auto &&edge : edges_) {
      if (edge->Contains(vrt1) and edge->Contains(vrt2)) {
        return edge;
      }
    }
    printf("Error finding edge in triangle %zu\n", index_);
    exit(1);
  }

  void set_half_edge(HalfEdge *h) { h_ = h; }
  HalfEdge *h_right() const { return h_; }
  void set_vertices(Vertex *v0, Vertex *v1, Vertex *v2) {
    vrts_[0] = v0;
    vrts_[1] = v1;
    vrts_[2] = v2;
  }
};

// HalfEdge //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct HalfEdge : public meshbrane::MeshBraneObject {

  ////////////////////////////
  // Fundamental attributes //
  ////////////////////////////
  /**
   * @brief Pointer to origin vertex
   */
  Vertex *v_{nullptr};
  /**
   * @brief Pointer to parallel edge
   */
  Edge *e_{nullptr};
  /**
   * @brief Pointer to face containing this half-edge
   */
  Triangle *f_{nullptr};
  /**
   * @brief Pointer to the twin half-edge
   */
  HalfEdge *h_twin_{nullptr};
  /**
   * @brief Pointer to the next half-edge in the face cycle
   */
  HalfEdge *h_next_{nullptr};

  ////////////////////////
  // Combinatorial maps //
  ////////////////////////
  /**
   * @brief Get pointer to twin half-edge
   */
  HalfEdge *h_twin() { return h_twin_; }
  /**
   * @brief Get pointer to next half-edge in the face cycle
   */
  HalfEdge *h_next() { return h_next_; }
  /**
   * @brief Get pointer to origin vertex
   */
  Vertex *v_origin() { return v_; }
  /**
   * @brief Get pointer to parallel edge
   */
  Edge *e_parallel() { return e_; }
  /**
   * @brief Get pointer to left face
   */
  Triangle *f_left() { return f_; }
  // Derived maps
  /**
   * @brief Get pointer to half-edge rotated clockwise about the origin vertex
   */
  HalfEdge *h_rotcw() { return this->h_twin_->h_next_; }
  /**
   * @brief Get pointer to the vertex at the head of the half-edge
   */
  Vertex *v_head() { return this->h_twin_->v_; }

  //////////////////////
  // Precomputed data //
  //////////////////////

  ////////////////////
  // Initialization //
  ////////////////////
  HalfEdge(size_t index) : mbrn::MeshBraneObject(index) {}

  ////////////////
  // Predicates //
  ////////////////
  bool is_in_some_negative_boundary() const { return f_->is_ghost(); }
  bool is_in_some_positive_boundary() const { return h_twin_->f_->is_ghost(); }

  friend bool operator==(const HalfEdge &lhs, const HalfEdge &rhs) {
    return (*lhs.v_ == *rhs.v_ and *lhs.e_ == *rhs.e_ and *lhs.f_ == *rhs.f_);
  }
  friend bool operator!=(const HalfEdge &lhs, const HalfEdge &rhs) {
    return !(lhs == rhs);
  }

  /////////////////////
  // Getters/Setters //
  /////////////////////
  void set_origin_vertex(Vertex *v) { v_ = v; }

  void set_parallel_edge(Edge *e) { e_ = e; }
  void set_left_face(Triangle *f) { f_ = f; }
  void set_twin_half_edge(HalfEdge *h) { h_twin_ = h; }
  void set_next_half_edge(HalfEdge *h) { h_next_ = h; }

  ////////////////
  // Generators //
  ////////////////
  HalfEdgeGenerator generate_H_rotcw();
};

// } // namespace half_edge
// #endif