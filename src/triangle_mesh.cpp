#include <cglass/triangle_mesh.hpp>

TriMesh::TriMesh() {

  double R_SYS{15};
  MakeIcosahedron();

  /*
  verts_.resize(6);
  for (int i_vrt{0}; i_vrt < verts_.size(); i_vrt++) {
    verts_[i_vrt].pos_[0] = 0.0;
    verts_[i_vrt].pos_[1] = 0.0;
    verts_[i_vrt].pos_[2] = 2 * i_vrt;
  }
  */
  //   verts_.emplace_back(10, 0, 0);
  //   verts_.emplace_back(0, 10, 0);
  //   verts_.emplace_back(0, 0, 10);
  //   tris_.emplace_back(&verts_[0], &verts_[1], &verts_[2]);
}

void TriMesh::ProjectToUnitSphere() {
  // p much just normalization i believe
  for (Vertex vrt : verts_) {
    double norm{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      norm += SQR(vrt.pos_[i_dim]);
    }
    norm = sqrt(norm);
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      double val{vrt.pos_[i_dim]};
      vrt.pos_[i_dim] = val / norm;
    }
  }
}

void TriMesh::MakeIcosahedron() {

  double r_sys{10.0};

  float phi = (1.0f + sqrt(5.0f)) * 0.5f; // golden ratio
  float a = r_sys * 1.0f;
  float b = r_sys * 1.0f / phi;
  // add vertices
  //   verts_.push_back({{0, b, -a}, base});
  verts_.emplace_back(0, b, -a);
  verts_.emplace_back(b, a, 0);
  verts_.emplace_back(-b, a, 0);
  verts_.emplace_back(0, b, a);
  verts_.emplace_back(0, -b, a);
  verts_.emplace_back(-a, 0, b);
  verts_.emplace_back(0, -b, -a);
  verts_.emplace_back(a, 0, -b);
  verts_.emplace_back(a, 0, b);
  verts_.emplace_back(-a, 0, -b);
  verts_.emplace_back(b, -a, 0);
  verts_.emplace_back(-b, -a, 0);
  // make ptrs to vertices to pass to triangle structures
  std::vector<Vertex *> vrt_ptrs_;
  for (int i_vrt{0}; i_vrt < verts_.size(); i_vrt++) {
    vrt_ptrs_.emplace_back(&verts_[i_vrt]);
  }
  // normalize radii of vertices
  //   ProjectToUnitSphere();
  // add triangles
  printf("%zu\n", vrt_ptrs_.size());
  tris_.emplace_back(vrt_ptrs_[2], vrt_ptrs_[1], vrt_ptrs_[0]);
  tris_.emplace_back(vrt_ptrs_[1], vrt_ptrs_[2], vrt_ptrs_[3]);
  tris_.emplace_back(vrt_ptrs_[5], vrt_ptrs_[4], vrt_ptrs_[3]);
  tris_.emplace_back(vrt_ptrs_[4], vrt_ptrs_[8], vrt_ptrs_[3]);
  tris_.emplace_back(vrt_ptrs_[7], vrt_ptrs_[6], vrt_ptrs_[0]);
  tris_.emplace_back(vrt_ptrs_[6], vrt_ptrs_[9], vrt_ptrs_[0]);
  tris_.emplace_back(vrt_ptrs_[11], vrt_ptrs_[10], vrt_ptrs_[4]);
  tris_.emplace_back(vrt_ptrs_[10], vrt_ptrs_[11], vrt_ptrs_[6]);
  tris_.emplace_back(vrt_ptrs_[9], vrt_ptrs_[5], vrt_ptrs_[2]);
  tris_.emplace_back(vrt_ptrs_[5], vrt_ptrs_[9], vrt_ptrs_[11]);
  tris_.emplace_back(vrt_ptrs_[8], vrt_ptrs_[7], vrt_ptrs_[1]);
  tris_.emplace_back(vrt_ptrs_[7], vrt_ptrs_[8], vrt_ptrs_[10]);
  tris_.emplace_back(vrt_ptrs_[2], vrt_ptrs_[5], vrt_ptrs_[3]);
  tris_.emplace_back(vrt_ptrs_[8], vrt_ptrs_[1], vrt_ptrs_[3]);
  tris_.emplace_back(vrt_ptrs_[9], vrt_ptrs_[2], vrt_ptrs_[0]);
  tris_.emplace_back(vrt_ptrs_[1], vrt_ptrs_[7], vrt_ptrs_[0]);
  tris_.emplace_back(vrt_ptrs_[11], vrt_ptrs_[9], vrt_ptrs_[6]);
  tris_.emplace_back(vrt_ptrs_[7], vrt_ptrs_[10], vrt_ptrs_[6]);
  tris_.emplace_back(vrt_ptrs_[5], vrt_ptrs_[11], vrt_ptrs_[4]);
  tris_.emplace_back(vrt_ptrs_[10], vrt_ptrs_[8], vrt_ptrs_[4]);
  /*
  mesh.add_triangle(v3, v2, v1);
  mesh.add_triangle(v2, v3, v4);
  mesh.add_triangle(v6, v5, v4);
  mesh.add_triangle(v5, v9, v4);
  mesh.add_triangle(v8, v7, v1);
  mesh.add_triangle(v7, v10, v1);
  mesh.add_triangle(v12, v11, v5);
  mesh.add_triangle(v11, v12, v7);
  mesh.add_triangle(v10, v6, v3);
  mesh.add_triangle(v6, v10, v12);
  mesh.add_triangle(v9, v8, v2);
  mesh.add_triangle(v8, v9, v11);
  mesh.add_triangle(v3, v6, v4);
  mesh.add_triangle(v9, v2, v4);
  mesh.add_triangle(v10, v3, v1);
  mesh.add_triangle(v2, v8, v1);
  mesh.add_triangle(v12, v10, v7);
  mesh.add_triangle(v8, v11, v7);
  mesh.add_triangle(v6, v12, v5);
  mesh.add_triangle(v11, v9, v5);
  */
}

void TriMesh::Draw(std::vector<graph_struct *> &graph_array) {
  for (int i_vrt{0}; i_vrt < verts_.size(); i_vrt++) {
    for (int i{0}; i < 3; i++) {
      verts_[i_vrt].g_.r[i] = verts_[i_vrt].pos_[i];
    }
    verts_[i_vrt].g_.color = M_PI;
    verts_[i_vrt].g_.diameter = 1.0;
    verts_[i_vrt].g_.length = 0.0;
    graph_array.push_back(&verts_[i_vrt].g_);
  }
}