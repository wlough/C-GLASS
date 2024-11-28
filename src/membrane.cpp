#include <cglass/filament.hpp>
// #include <cglass/ply_tools.hpp>
// #include "cglass/triangle_mesh.hpp"
#include <cglass/membrane.hpp>
#include <filesystem>
#include <meshbrane/matrix_mesh.hpp>
#include <meshbrane/meshbrane_data_types.hpp>
#include <unistd.h>
#include <vector>

// void Membrane::LoadPly() {
//   // printf("Loading ply file\n");
//   printf("Loading ply file %s\n", ply_path.c_str());
//   std::filesystem::path path(ply_path);
//   std::string directory = path.parent_path().string();
//   std::string filename = path.filename().string();
//   meshbrane::MatrixMesh m = meshbrane::MatrixMesh::from_he_ply(ply_path);
//   int euler_characteristic = m.get_euler_characteristic();
//   int num_vertices = m.get_num_vertices();
//   int num_faces = m.get_num_faces();
//   int num_edges = m.get_num_edges();

//   auto [xyz_coord_V, V_cycle_E, V_cycle_F] = m.vef_samples();
//   meshbrane::Samples3d xyz_coord_V2 = xyz_coord_V;
//   // for (int i = 0; i < num_vertices; i++) {
//   //   xyz_coord_V(i, 1) = xyz_coord_V2(i, 2);
//   //   xyz_coord_V(i, 2) = xyz_coord_V2(i, 0);
//   //   xyz_coord_V(i, 0) = xyz_coord_V2(i, 1);
//   // }
//   // m.set_xyz_coord_V(xyz_coord_V);
//   m.write_he_ply("output/" + filename);

//   tris_.reserve(num_faces);
//   edges_.reserve(num_edges);
//   vrts_.reserve(num_vertices);
//   printf("%zu faces, %zu edges, %zu verts\n", num_faces, num_edges,
//          num_vertices);
//   // MakeIcosahedron();
//   for (int v = 0; v < num_vertices; v++) {
//     // double xyz_coord[3] = {xyz_coord_V(v, 0), xyz_coord_V(v, 1), xyz_coord_V(v, 2)};
//     vrts_.emplace_back(xyz_coord_V(v, 0), xyz_coord_V(v, 1), xyz_coord_V(v, 2));
//   }
//   for (int e = 0; e < num_edges; e++) {
//     int v0 = V_cycle_E(e, 0);
//     int v1 = V_cycle_E(e, 1);
//     edges_.emplace_back(&vrts_[v0], &vrts_[v1]);
//   }
//   for (int f = 0; f < num_faces; f++) {
//     int v0 = V_cycle_F(f, 0);
//     int v1 = V_cycle_F(f, 1);
//     int v2 = V_cycle_F(f, 2);
//     tris_.emplace_back(&vrts_[v0], &vrts_[v1], &vrts_[v2]);
//   }
// }