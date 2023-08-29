#include <cglass/triangle_mesh.hpp>

TriMesh::TriMesh() {

  double R_SYS{15};
  verts_.reserve(1024);
  tris_.reserve(1024);
  MakeIcosahedron();
  DivideFaces();
  ProjectToUnitSphere();
  //   DivideFaces();
  //   ProjectToUnitSphere();
  //   DivideFaces();
  //   ProjectToUnitSphere();
}

void TriMesh::ProjectToUnitSphere() {
  double R_SYS{15};
  //   return;
  // p much just normalization i believe
  for (auto &&vrt : verts_) {
    double norm{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      norm += SQR(vrt.pos_[i_dim]);
    }
    norm = sqrt(norm);
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      double val{vrt.pos_[i_dim]};
      vrt.pos_[i_dim] = R_SYS * val / norm;
    }
  }
}

void TriMesh::MakeIcosahedron() {

  double r_sys{15.0};

  float phi = (1.0f + sqrt(5.0f)) * 0.5f; // golden ratio
  float a = r_sys * 1.0f;
  float b = r_sys * 1.0f / phi;
  // add vertices
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
  //   std::vector<Vertex *> vrt_ptrs_;
  for (int i_vrt{0}; i_vrt < verts_.size(); i_vrt++) {
    // vrt_ptrs_.emplace_back(&verts_[i_vrt]);
    verts_[i_vrt].g_.color = M_PI;
    verts_[i_vrt].g_.draw = draw_type::fixed;
    verts_[i_vrt].g_.diameter = 1.0;
  }
  // normalize radii of vertices
  ProjectToUnitSphere();
  // add triangles
  tris_.emplace_back(&verts_[2], &verts_[1], &verts_[0]);
  tris_.emplace_back(&verts_[1], &verts_[2], &verts_[3]);
  tris_.emplace_back(&verts_[5], &verts_[4], &verts_[3]);
  tris_.emplace_back(&verts_[4], &verts_[8], &verts_[3]);
  tris_.emplace_back(&verts_[7], &verts_[6], &verts_[0]);
  tris_.emplace_back(&verts_[6], &verts_[9], &verts_[0]);
  tris_.emplace_back(&verts_[11], &verts_[10], &verts_[4]);
  tris_.emplace_back(&verts_[10], &verts_[11], &verts_[6]);
  tris_.emplace_back(&verts_[9], &verts_[5], &verts_[2]);
  tris_.emplace_back(&verts_[5], &verts_[9], &verts_[11]);
  tris_.emplace_back(&verts_[8], &verts_[7], &verts_[1]);
  tris_.emplace_back(&verts_[7], &verts_[8], &verts_[10]);
  tris_.emplace_back(&verts_[2], &verts_[5], &verts_[3]);
  tris_.emplace_back(&verts_[8], &verts_[1], &verts_[3]);
  tris_.emplace_back(&verts_[9], &verts_[2], &verts_[0]);
  tris_.emplace_back(&verts_[1], &verts_[7], &verts_[0]);
  tris_.emplace_back(&verts_[11], &verts_[9], &verts_[6]);
  tris_.emplace_back(&verts_[7], &verts_[10], &verts_[6]);
  tris_.emplace_back(&verts_[5], &verts_[11], &verts_[4]);
  tris_.emplace_back(&verts_[10], &verts_[8], &verts_[4]);
}

void TriMesh::DivideFaces() {

  std::vector<Triangle> new_faces;

  for (auto &&face : tris_) {
    Vertex *vrt0{face.vrts_[0]};
    Vertex *vrt1{face.vrts_[1]};
    Vertex *vrt2{face.vrts_[2]};
    double pos10[3] = {(vrt0->pos_[0] + vrt1->pos_[0]) / 2.0,
                       (vrt0->pos_[1] + vrt1->pos_[1]) / 2.0,
                       (vrt0->pos_[2] + vrt1->pos_[2]) / 2.0};
    double pos12[3] = {(vrt2->pos_[0] + vrt1->pos_[0]) / 2.0,
                       (vrt2->pos_[1] + vrt1->pos_[1]) / 2.0,
                       (vrt2->pos_[2] + vrt1->pos_[2]) / 2.0};
    double pos20[3] = {(vrt0->pos_[0] + vrt2->pos_[0]) / 2.0,
                       (vrt0->pos_[1] + vrt2->pos_[1]) / 2.0,
                       (vrt0->pos_[2] + vrt2->pos_[2]) / 2.0};
    Vertex *vrt10{nullptr}, *vrt12{nullptr}, *vrt20{nullptr};
    bool dup10{false}, dup12{false}, dup20{false};
    for (auto &&vrt : verts_) {
      if (vrt.pos_[0] == pos10[0] and vrt.pos_[1] == pos10[1] and
          vrt.pos_[2] == pos10[2]) {
        dup10 = true;
        vrt10 = &vrt;
      }
      if (vrt.pos_[0] == pos12[0] and vrt.pos_[1] == pos12[1] and
          vrt.pos_[2] == pos12[2]) {
        dup12 = true;
        vrt12 = &vrt;
      }
      if (vrt.pos_[0] == pos20[0] and vrt.pos_[1] == pos20[1] and
          vrt.pos_[2] == pos20[2]) {
        dup20 = true;
        vrt20 = &vrt;
      }
    }
    if (!dup10) {
      verts_.emplace_back(pos10[0], pos10[1], pos10[2]);
      vrt10 = &verts_.back();
      verts_.back().g_.color = 1 * M_PI;
      //   verts_.back().g_.diameter = 2;
    }
    if (!dup12) {
      verts_.emplace_back(pos12[0], pos12[1], pos12[2]);
      vrt12 = &verts_.back();
      verts_.back().g_.color = 2 * M_PI;
      //   verts_.back().g_.diameter = 3;
    }
    if (!dup20) {
      verts_.emplace_back(pos20[0], pos20[1], pos20[2]);
      vrt20 = &verts_.back();
      verts_.back().g_.color = 3 * M_PI;
      //   verts_.back().g_.diameter = 4;
    }
    if (vrt10 == nullptr or vrt12 == nullptr or vrt20 == nullptr) {
      printf("nope\n");
      exit(1);
    }
    new_faces.emplace_back(vrt0, vrt10, vrt20);
    new_faces.emplace_back(vrt1, vrt10, vrt12);
    new_faces.emplace_back(vrt2, vrt12, vrt20);
    new_faces.emplace_back(vrt10, vrt12, vrt20);
  }
  tris_ = new_faces;
}

void TriMesh::Draw(std::vector<graph_struct *> &graph_array) {

  for (int i_vrt{0}; i_vrt < verts_.size(); i_vrt++) {
    for (int i{0}; i < 3; i++) {
      verts_[i_vrt].g_.r[i] = verts_[i_vrt].pos_[i];
    }
    // verts_[i_vrt].g_.color = 1.5 * M_PI;
    // verts_[i_vrt].g_.diameter = 1;
    // verts_[i_vrt].g_.length = 0.0;
    graph_array.push_back(&verts_[i_vrt].g_);
  }
  //   DivideFaces();
}