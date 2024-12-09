// #ifndef _CGLASS_TRIANGLE_MESH_
// #define _CGLASS_TRIANGLE_MESH_
#pragma once
// #include "common_libs.hpp"
// #include "definitions.hpp"
#include "half_edge_primitives.hpp"
#include "meshbrane/matrix_mesh.hpp"
#include "meshbrane/meshbrane_data_types.hpp"
#include "meshbrane/meshbrane_object.hpp"
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
  bool scale_to_system_radius{false};
  bool make_a_movie{false};
  double target_face_area_{0.0};
  double target_volume_{0.0};

  std::string get_new_output_ply_path() {
    // printf("  getting new output ply path\n");
    std::string directory = "output/ply_files/";
    std::string filename = params_->run_name;
    // output path is directory + filename + zero_pad_to_6_digits + i_datapoint_ + ".ply"
    // number of digits in i_datapoint_
    size_t n_digits = 1;
    if (i_datapoint_ > 0)
      n_digits = (size_t)std::log10(i_datapoint_) + 1;
    size_t num_zeros = 6 - n_digits;
    std::string output_ply_path = directory + filename + "_";
    for (size_t i = 0; i < num_zeros; ++i) {
      output_ply_path += "0";
    }
    output_ply_path += std::to_string(i_datapoint_) + ".ply";
    return output_ply_path;
  }
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
  void InitializeHalfEdgeMats();
  /**
   * @brief Copy constructor
   */
  TriMesh(const TriMesh &other) : meshbrane::MatrixMesh(other) {};
  /**
   * @brief Load half-edge data from a ply file. Does not update vrts_, edges_, tris_, half_edges_, boundaries_ lists.
   */
  void LoadPly();
  void MakeIcosphere();
  void MakeIcosahedron();
  void DivideFaces();
  void ProjectToUnitSphere();
  /**
 * @brief Sync vrts_, edges_, tris_, half_edges_, boundaries_ lists with data in half-edge matrices. Does NOT update cached incidence/adjacency/geometric data.
 */
  void SyncWithMats();
  void SyncIndices();
  void VEFdataToMats();
  void ScaleToSystemRadius();

  void RefreshFromMatsBack();
  void RefreshMats();

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
  /**
   * @brief Update incidence/adjacency data for vrts_, edges_, tris_
   */
  void UpdateIncidenceData();
  /**
   * @brief Update centroid and geometric data (lengths, areas, volumes, tangent/normal vectors) for edges_ and tris_.
   */
  void UpdateGeometricData();
  void RefreshTetherParams();
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
  // private:
  void DrawVerts();
  void FlipEdges();
  void UpdateCentroid();
  void UpdateTriangles();
  void UpdateMesh();
  void ApplyMembraneForces();
  void ApplyBoundaryForces();
  void ApplyBendingForces();
  // public:
  void Draw(std::vector<graph_struct *> &graph_array);

  //////////////////////
  // ???????????????? //
  //////////////////////
  // public:
  std::vector<graph_struct> f_mem_;
  graph_struct o_;
  RNG *rng_; // SF TODO link with system RNG
  MinimumDistance mindist_;
  /**
   * @brief Updates Euler angles???
   * 
   */
  void UpdateEulerAngles();

  ///////////////
  // Debugging //
  ///////////////
  void TestFun();
  int CheckPtrs();
  int CheckMats();

private:
  //////////////////////
  // To be deprecated //
  //////////////////////
  void InitializeMeshOG();
  void UpdateNeighborsOG();
};

// #endif