#include <cglass/triangle_mesh.hpp>

// todo lol fix seed
TriMesh::TriMesh() : rng_(0) {

  return; // need to ensure this isn't initialized without proper boundary type

  double R_SYS{15};
  vrts_.reserve(1280);
  tris_.reserve(1280);
  MakeIcosahedron(); // 20 triangle faces initially
  DivideFaces();     // 80 triangles
  DivideFaces();     // 320
  // DivideFaces();     // 1280
  ProjectToUnitSphere();
  for (auto &&vrt : vrts_) {
    vrt.neighbs_.resize(6);
    vrt.neighb_int_.resize(6);
  }
  UpdateNeighbors();
  size_t n_flawed{0};
  size_t n_gucci{0};
  for (auto &&vrt : vrts_) {
    if (vrt.i_neighb_ != 6) {
      n_flawed++;
    } else if (vrt.i_neighb_ == 6) {
      n_gucci++;
    } else {
      printf("ummm?? in TRIANGLE MESH INIT\n");
      exit(1);
    }
  }
  printf("%zu vrts (%zu flawed; %zu good)\n", vrts_.size(), n_flawed, n_gucci);
}

void TriMesh::UpdateNeighbors() {

  // Reset storage of neighb and triangle ptrs in each vertex
  for (auto &&vrt : vrts_) {
    vrt.i_tri_ = 0;
    vrt.i_neighb_ = 0;
  }
  // Update triangle ptrs stored by each vertex
  for (auto &&tri : tris_) {
    tri.vrts_[0]->tris_[tri.vrts_[0]->i_tri_++] = &tri;
    tri.vrts_[1]->tris_[tri.vrts_[1]->i_tri_++] = &tri;
    tri.vrts_[2]->tris_[tri.vrts_[2]->i_tri_++] = &tri;
  }
  // verify
  for (auto &&vrt : vrts_) {
    // for a defect-free mesh, all points should incorporate 6 triangles
    // mAYBE? i dont really know
    printf("I compose %i triangles!\n", vrt.i_tri_);
  }
  // Update vertex neighbors held by each vertex
  for (auto &&vrt : vrts_) {
    // Each triangle we are a part of will contain 2 neighbors
    for (int i_tri{0}; i_tri < vrt.i_tri_; i_tri++) {
      Triangle *tri{vrt.tris_[i_tri]};
      // Add all points that are not ourselves or duplicated vertices
      for (int i{0}; i < 3; i++) {
        if (*tri->vrts_[i] != vrt) {
          printf("NOT EQUAL\n");
          bool duplicate{false};
          for (int i_neighb{0}; i_neighb < vrt.i_neighb_; i_neighb++) {
            if (tri->vrts_[i] == vrt.neighbs_[i_neighb])
              duplicate = true;
          }
          if (!duplicate) {
            vrt.neighbs_[vrt.i_neighb_++] = tri->vrts_[i];
            printf("ADDED\n");
          } else {
            printf("DUPLICATE IGNORED\n");
          }
        } else if (*tri->vrts_[i] == vrt) {
          printf("EQUAL\n");
        }
      }
    }
    printf("%i neighbors from %i triangles\n\n", vrt.i_neighb_, vrt.i_tri_);
  }

  /*
  return;
  for (auto &&tri : tris_) {
    double l1{0.0};
    double l2{0.0};
    double l3{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      l1 += SQR(tri.vrts_[0]->pos_[i_dim] - tri.vrts_[1]->pos_[i_dim]);
      l2 += SQR(tri.vrts_[1]->pos_[i_dim] - tri.vrts_[2]->pos_[i_dim]);
      l3 += SQR(tri.vrts_[0]->pos_[i_dim] - tri.vrts_[2]->pos_[i_dim]);
    }
    l1 = sqrt(l1);
    l2 = sqrt(l2);
    l3 = sqrt(l3);
  }
  for (auto &&vrt : vrts_) {
    int i_neighb{0};
    for (auto &&neighb : vrts_) {
      double dist{0.0};
      for (int i{0}; i < 3; i++) {
        dist += SQR(vrt.pos_[i] - neighb.pos_[i]);
      }
      dist = sqrt(dist);
      if (dist == 0.0) {
        continue;
      }
      if (dist < 5.0) {
        if (i_neighb > 6) {
          printf("MESH IRREGULARITY\n");
          // exit(1);
        }
        vrt.neighbs_[i_neighb++] = &neighb;
        // printf("%i\n", i_neighb);
      }
    }
  }
  */
}

void TriMesh::ProjectToUnitSphere() {
  double R_SYS{15.0};
  // p much just normalization i believe
  for (auto &&vrt : vrts_) {
    // const double *old_pos{vrt.GetPosition()};
    double norm{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      norm += SQR(vrt.pos_[i_dim]);
    }
    norm = sqrt(norm);
    double new_pos[3];
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      double val{vrt.pos_[i_dim]};
      new_pos[i_dim] = R_SYS * val / norm;
      // printf("%g = %g * %g / %g\n", new_pos[i_dim], R_SYS, val, norm);
    }
    vrt.SetPos(new_pos);
  }
}

void TriMesh::MakeIcosahedron() {

  // sf todo update system radius from params
  double r_sys{15.0};
  double phi = (1.0f + sqrt(5.0f)) * 0.5f; // golden ratio
  double a = r_sys * 1.0f;
  double b = r_sys * 1.0f / phi;
  // add vertices
  vrts_.emplace_back(0.0, b, -a);
  vrts_.emplace_back(b, a, 0.0);
  vrts_.emplace_back(-b, a, 0.0);
  vrts_.emplace_back(0.0, b, a);
  vrts_.emplace_back(0.0, -b, a);
  vrts_.emplace_back(-a, 0.0, b);
  vrts_.emplace_back(0.0, -b, -a);
  vrts_.emplace_back(a, 0.0, -b);
  vrts_.emplace_back(a, 0.0, b);
  vrts_.emplace_back(-a, 0.0, -b);
  vrts_.emplace_back(b, -a, 0.0);
  vrts_.emplace_back(-b, -a, 0.0);

  tris_.emplace_back(&vrts_[2], &vrts_[1], &vrts_[0]);
  tris_.emplace_back(&vrts_[1], &vrts_[2], &vrts_[3]);
  tris_.emplace_back(&vrts_[5], &vrts_[4], &vrts_[3]);
  tris_.emplace_back(&vrts_[4], &vrts_[8], &vrts_[3]);
  tris_.emplace_back(&vrts_[7], &vrts_[6], &vrts_[0]);
  tris_.emplace_back(&vrts_[6], &vrts_[9], &vrts_[0]);
  tris_.emplace_back(&vrts_[11], &vrts_[10], &vrts_[4]);
  tris_.emplace_back(&vrts_[10], &vrts_[11], &vrts_[6]);
  tris_.emplace_back(&vrts_[9], &vrts_[5], &vrts_[2]);
  tris_.emplace_back(&vrts_[5], &vrts_[9], &vrts_[11]);
  tris_.emplace_back(&vrts_[8], &vrts_[7], &vrts_[1]);
  tris_.emplace_back(&vrts_[7], &vrts_[8], &vrts_[10]);
  tris_.emplace_back(&vrts_[2], &vrts_[5], &vrts_[3]);
  tris_.emplace_back(&vrts_[8], &vrts_[1], &vrts_[3]);
  tris_.emplace_back(&vrts_[9], &vrts_[2], &vrts_[0]);
  tris_.emplace_back(&vrts_[1], &vrts_[7], &vrts_[0]);
  tris_.emplace_back(&vrts_[11], &vrts_[9], &vrts_[6]);
  tris_.emplace_back(&vrts_[7], &vrts_[10], &vrts_[6]);
  tris_.emplace_back(&vrts_[5], &vrts_[11], &vrts_[4]);
  tris_.emplace_back(&vrts_[10], &vrts_[8], &vrts_[4]);
}

void TriMesh::DivideFaces() {

  std::vector<Triangle> new_faces;
  for (auto &&face : tris_) {
    Vertex *vrt0{face.vrts_[0]};
    Vertex *vrt1{face.vrts_[1]};
    Vertex *vrt2{face.vrts_[2]};
    // SF todo validate this
    double pos10[3] = {(vrt0->pos_[0] + vrt1->pos_[0]) / 2.0,
                       (vrt0->pos_[1] + vrt1->pos_[1]) / 2.0,
                       (vrt0->pos_[2] + vrt1->pos_[2]) / 2.0};
    double pos12[3] = {(vrt2->pos_[0] + vrt1->pos_[0]) / 2.0,
                       (vrt2->pos_[1] + vrt1->pos_[1]) / 2.0,
                       (vrt2->pos_[2] + vrt1->pos_[2]) / 2.0};
    double pos20[3] = {(vrt0->pos_[0] + vrt2->pos_[0]) / 2.0,
                       (vrt0->pos_[1] + vrt2->pos_[1]) / 2.0,
                       (vrt0->pos_[2] + vrt2->pos_[2]) / 2.0};

    double norm1{0.0}, norm2{0.0}, norm3{0.0};
    for (int i{0}; i < 3; i++) {
      norm1 += SQR(pos10[i]);
      norm2 += SQR(pos12[i]);
      norm3 += SQR(pos20[i]);
    }
    norm1 = sqrt(norm1);
    norm2 = sqrt(norm2);
    norm3 = sqrt(norm3);
    Vertex *vrt10{nullptr}, *vrt12{nullptr}, *vrt20{nullptr};
    bool dup10{false}, dup12{false}, dup20{false};
    for (auto &&vrt : vrts_) {
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
      vrts_.emplace_back(pos10[0], pos10[1], pos10[2]);
      vrt10 = &vrts_.back();
    }
    if (!dup12) {
      vrts_.emplace_back(pos12[0], pos12[1], pos12[2]);
      vrt12 = &vrts_.back();
    }
    if (!dup20) {
      vrts_.emplace_back(pos20[0], pos20[1], pos20[2]);
      vrt20 = &vrts_.back();
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
  printf("%zu -> %zu triangles\n", tris_.size(), new_faces.size());
  tris_ = new_faces;
}

void TriMesh::Draw(std::vector<graph_struct *> &graph_array) {

  UpdatePositions();
  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    vrts_[i_vrt].Draw(graph_array);
  }
}

void TriMesh::UpdatePositions() {

  // sum forces over all vertices
  for (auto &&vrt : vrts_) {
    vrt.ZeroForce();
    for (int i_neighb{0}; i_neighb < vrt.i_neighb_; i_neighb++) {
      vrt.neighb_int_[i_neighb] = true;
      Vertex *neighb{vrt.neighbs_[i_neighb]};
      if (neighb == nullptr) {
        printf("what in TriMesh::UpdatePosition()\n");
        exit(1);
      }

      double rhat[3];
      double rmag{0.0};
      for (int i{0}; i < 3; i++) {
        rhat[i] = vrt.GetPosition()[i] - neighb->GetPosition()[i];
        rmag += SQR(rhat[i]);
      }
      rmag = sqrt(rmag);
      double fmag{0.01 * SQR(rmag - 4.5)};
      double f[3];
      for (int i{0}; i < 3; i++) {
        rhat[i] = rhat[i] / rmag;
        f[i] = -1 * fmag * rhat[i];
      }
      vrt.AddForce(f);
      neighb->SubForce(f);
    }
  }
  // apply displacements
  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    //Expected diffusion length of the crosslink in the solution in delta
    double sigma_d = sqrt(2.0 * 0.0005 * 0.0131);
    //Distance actually diffused
    double dr[3] = {0.0, 0.0, 0.0};
    //Previous and final location of the crosslink
    double const *const r_prev = vrts_[i_vrt].GetPosition();
    double r_final[3];
    //Update position of the crosslink
    for (int i = 0; i < 3; ++i) {
      dr[i] = rng_.RandomNormal(sigma_d);
      double f = vrts_[i_vrt].GetForce()[i];
      r_final[i] = r_prev[i] + dr[i] + 0.00001 * f;
    }
    vrts_[i_vrt].SetPos(r_final);
    // vrts_[i_vrt].SetPositionXYZ(r_final[0], r_final[1], r_final[2]);
  }
}