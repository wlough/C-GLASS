// #ifndef _CGLASS_TRIANGLE_MESH_
// #define _CGLASS_TRIANGLE_MESH_
#pragma once
// #include "common_libs.hpp"
// #include "definitions.hpp"
#include "half_edge_primitives.hpp"
#include "meshbrane/matrix_mesh.hpp"
#include "meshbrane/meshbrane_data_types.hpp"
#include "meshbrane/meshbrane_object.hpp"
#include "meshbrane/simple_generator.hpp"
#include "minimum_distance.hpp"
#include "rng.hpp"
#include "site.hpp"
#include <memory>   // std::shared_ptr
#include <optional> // std::optional

// TODO add param storage and sync site seeds

/**
 * @brief A dynamically triangulated surface
 * 
 */
class TriMesh : public meshbrane::MatrixMesh {
  ////////////////////////////
  // Fundamental attributes //
  ////////////////////////////
public:
  std::vector<Vertex> vrts_;
  std::vector<Triangle> tris_;
  std::vector<Edge> edges_;
  std::vector<HalfEdge> half_edges_;
  std::vector<Triangle> boundaries_;
  std::vector<Object *> boundary_neighbs_;
  system_parameters *params_{nullptr};
  std::string ply_path{"none"};
  double r_sys_{0.0};
  // Bending force
  double bending_modulus_{0.0};
  double spontaneous_curvature_{0.0};
  double splay_modulus_{0.0};
  // Tether force
  double dimensionless_tether_repulsive_onset_{0.8};
  double dimensionless_tether_repulsive_singularity_{0.25};
  double dimensionless_tether_attractive_onset_{1.2};
  double dimensionless_tether_attractive_singularity_{2.5};
  double tether_stiffness_{0.0};
  double tether_attractive_singularity_{0.0};
  double tether_repulsive_singularity_{0.0};
  double tether_attractive_onset_{0.0};
  double tether_repulsive_onset_{0.0};
  // Area conservation force
  double area_reg_stiffness_{0.0};
  double preferred_area_{1.0};
  // Volume conservation force
  double preferred_volume_{0.09403159725796};
  double volume_reg_stiffness_{0.0};
  // Drag force
  double node_drag_coefficient_{0.0};
  //  Edge flipping
  int flip_sweeps_per_step_{1};
  double flipping_probability_{0.3};
  //  Misc
  double timestep_{1e-5};
  static const size_t n_edges_min_{3};  // true for any connected graph
  static const size_t n_edges_max_{10}; // arbitrary choice
  bool do_not_pass_go_{false};
  int i_datapoint_{0};

  ////////////////////
  // Initialization //
  ////////////////////
  TriMesh() {}
  /**
   * @brief calls SetParameters, LoadPly/MakeIcosphere, InitializeMesh
   */
  void Init(system_parameters *params);
  /**
   * @brief Set parameters from system_parameters
   */
  void SetParameters();
  /**
   * @brief Initialize MatrixMesh from a ply file.
   */
  void InitializeMesh();
  /**
   * @brief Copy constructor
   */
  TriMesh(const TriMesh &other) : meshbrane::MatrixMesh(other) {};
  void LoadPly();
  void MakeIcosphere();
  void MakeIcosahedron();
  void DivideFaces();
  void ProjectToUnitSphere();
  /**
   * @brief Refresh vertex, edge, and face lists from the matrix mesh data
   */
  void RefreshFromMats();
  void RefreshMats();
  void TestFun();
  //////////////////////
  // Precomputed data //
  //////////////////////
  double average_face_area_{0.0};
  double average_face_volume_{0.0};
  double centroid_[3];
  double f_avgs_[4];                // indices 0-4: tether, bend, area, vol
  double average_edge_length_{0.0}; // average edge length

  //////////////////////////////
  // Getters/Setters/updaters //
  //////////////////////////////
  void RefreshEdgeParams();
  void UpdatePositions();

  ///////////////////////////
  // Output methods/params //
  //////////////////////////
  FILE *forces_{nullptr};    // average force from each type of potential
  FILE *vertices_{nullptr};  // position (3D per vrt per step)
  FILE *adjacency_{nullptr}; // adjacency matrix  (2D per vrt per step)
  void WriteOutputs();

  ///////////////////
  // Other methods //
  ///////////////////
private:
  void FlipEdges();
  void UpdateCentroid();
  void UpdateTriangles();
  void UpdateNeighbors();
  void UpdateMesh();
  void ApplyMembraneForces();
  void ApplyBoundaryForces();

public:
  void Draw(std::vector<graph_struct *> &graph_array);

  //////////////////////
  // ???????????????? //
  //////////////////////
public:
  std::vector<graph_struct> f_mem_;
  graph_struct o_;
  RNG *rng_; // SF TODO link with system RNG
  MinimumDistance mindist_;

  //////////////////////
  // To be deprecated //
  //////////////////////
private:
  void InitializeMeshOG();
};

// #endif