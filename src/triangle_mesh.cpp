#include <cglass/triangle_mesh.hpp>

// todo lol fix seed
TriMesh::TriMesh() : rng_(0) {

  // return; // need to ensure this isn't initialized without proper boundary type

  // fix these got dang parameters and put em in the yaml file
  double R_SYS{15};
  kappa_ = 10.0;
  kappa_B_ = 10.0;
  // l_max_ = 20;
  // l_min_ = 9;
  // l_c0_ = 19;
  // l_c1_ = 18;
  l_max_ = 6;
  l_min_ = 3;
  l_c0_ = 4.7;
  l_c1_ = 4.3;

  vrts_.reserve(1280);
  tris_.reserve(1280);
  MakeIcosahedron(); // 20 triangle faces initially
  DivideFaces();     // 80 triangles
  DivideFaces();     // 320
  // DivideFaces();     // 1280
  ProjectToUnitSphere();
  for (auto &&vrt : vrts_) {
    vrt.neighbs_.resize(6);
    // vrt.neighb_int_.resize(6);
  }
  UpdateNeighbors();
  size_t n_flawed{0};
  size_t n_gucci{0};
  for (auto &&vrt : vrts_) {
    if (vrt.n_neighbs_ != 6) {
      n_flawed++;
    } else if (vrt.n_neighbs_ == 6) {
      n_gucci++;
    } else {
      printf("ummm?? in TRIANGLE MESH INIT\n");
      exit(1);
    }
  }
  printf("%zu vrts (%zu flawed; %zu ideal)\n", vrts_.size(), n_flawed, n_gucci);
}

void TriMesh::UpdateNeighbors() {

  // This function is purposefully coded to be explicit yet inefficient
  // (it does not get called often if at all beyond initialization)
  // Reset storage of neighb and triangle ptrs in each vertex
  for (auto &&vrt : vrts_) {
    vrt.n_tris_ = 0;
    vrt.n_neighbs_ = 0;
  }
  // Update triangle ptrs stored by each vertex
  for (auto &&tri : tris_) {
    tri.vrts_[0]->tris_[tri.vrts_[0]->n_tris_++] = &tri;
    tri.vrts_[1]->tris_[tri.vrts_[1]->n_tris_++] = &tri;
    tri.vrts_[2]->tris_[tri.vrts_[2]->n_tris_++] = &tri;
  }
  // Update vertex neighbors held by each vertex
  for (auto &&vrt : vrts_) {
    // Each triangle we are a part of will contain 2 neighbors
    for (int i_tri{0}; i_tri < vrt.n_tris_; i_tri++) {
      Triangle *tri{vrt.tris_[i_tri]};
      // Add all points that are not ourselves or duplicated vertices
      for (int i{0}; i < 3; i++) {
        if (*tri->vrts_[i] != vrt) {
          // printf("NOT EQUAL\n");
          bool duplicate{false};
          for (int i_neighb{0}; i_neighb < vrt.n_neighbs_; i_neighb++) {
            if (tri->vrts_[i] == vrt.neighbs_[i_neighb])
              duplicate = true;
          }
          if (!duplicate) {
            vrt.neighbs_[vrt.n_neighbs_++] = tri->vrts_[i];
            // printf("ADDED\n");
          }
          // else {
          // printf("DUPLICATE IGNORED\n");
          // }
        }
        // else if (*tri->vrts_[i] == vrt) {
        // printf("EQUAL\n");
        // }
      }
    }
    // printf("%i neighbors from %i triangles\n\n", vrt.n_neighbs_, vrt.n_tris_);
  }
  // rearrange all neighbor lists to order them in a counterclockwise ring
  // SF TODO: need to arrange in counter-clockwise ring order
  for (auto &&vrt : vrts_) {
    Vertex *neighbs_ordered[vrt.n_neighbs_]; // list to be ordered in ccw ring
    neighbs_ordered[0] = vrt.neighbs_[0];    // use 1st entry as starting point
    // scan over all neighbors to find which is next in the ring
    for (int i_neighb{0}; i_neighb < vrt.n_neighbs_ - 1; i_neighb++) {
      Vertex *current_entry{neighbs_ordered[i_neighb]};
      double r_ij[3]; // points from j (current_entry) to i (vrt)
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        r_ij[i_dim] = vrt.pos_[i_dim] - current_entry->pos_[i_dim];
      }
      bool found{false};
      // the next node in the ring order will be a neighbor of current entry
      for (int j_neighb{0}; j_neighb < current_entry->n_neighbs_; j_neighb++) {
        if (found) {
          break;
        }
        Vertex *entry_neighb{current_entry->neighbs_[j_neighb]};
        for (int k_neighb{0}; k_neighb < vrt.n_neighbs_; k_neighb++) {
          Vertex *this_entry{vrt.neighbs_[k_neighb]};
          if (entry_neighb == this_entry) {
            double r_ik[3]; // points from k (entry_neighb) to i (vrt)
            for (int i_dim{0}; i_dim < 3; i_dim++) {
              r_ik[i_dim] = vrt.pos_[i_dim] - this_entry->pos_[i_dim];
            }
            // to ensure ccw ordering, need to make sure n_tri points away from origin
            double nhat[3];
            cross_product(r_ij, r_ik, nhat, 3);
            // if nhat is antiparallel with vertex position, its pointing inwards
            if (dot_product(3, vrt.pos_, nhat) == 0.0) {
            }
            if (dot_product(3, vrt.pos_, nhat) < 0.0) {
              continue;
            }
            if (this_entry != current_entry and this_entry != &vrt) {
              if (i_neighb == 0) {
                neighbs_ordered[i_neighb + 1] = this_entry;
                found = true;
                break;
              } else if (this_entry != neighbs_ordered[i_neighb - 1]) {
                neighbs_ordered[i_neighb + 1] = this_entry;
                found = true;
                break;
              }
            }
          }
        }
      }
      if (!found) {
        printf("tri mesh neighbor not found ?? \n");
        exit(1);
      }
    }
    for (int i_neighb{0}; i_neighb < vrt.n_neighbs_; i_neighb++) {
      vrt.neighbs_[i_neighb] = neighbs_ordered[i_neighb];
    }
  }
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

  printf("initializing icosahedron (20 triangles; 12 points)\n");
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

  double edge_length{0.0};
  for (int i{0}; i < 3; i++) {
    edge_length += SQR(vrts_[0].pos_[i] - vrts_[1].pos_[i]);
  }
  edge_length = sqrt(edge_length);
  printf("edge length is %g\n", edge_length);

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
  double avg_edge_length{0.0};
  size_t n_edges{0};
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
      norm1 += SQR(vrt1->pos_[i] - pos10[i]);
      norm2 += SQR(vrt2->pos_[i] - pos12[i]);
      norm3 += SQR(vrt2->pos_[i] - pos20[i]);
      // norm1 += SQR(pos10[i]);
      // norm2 += SQR(pos12[i]);
      // norm3 += SQR(pos20[i]);
    }
    norm1 = sqrt(norm1);
    norm2 = sqrt(norm2);
    norm3 = sqrt(norm3);
    avg_edge_length += norm1;
    avg_edge_length += norm2;
    avg_edge_length += norm3;
    n_edges += 3;
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
  printf("avg edge length is %g\n", avg_edge_length / n_edges);
  tris_ = new_faces;
}

void TriMesh::Draw(std::vector<graph_struct *> &graph_array) {

  UpdatePositions();
  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    vrts_[i_vrt].Draw(graph_array);
  }
}

void TriMesh::UpdatePositions() {

  // SF TODO all calculations are currently done redundantly
  // SF TODO i.e., we do not take advtange of newton's 3rd law
  // SF TODO once validated, increase computational efficiency by addressing this

  // zero forces out first
  for (auto &&vrt : vrts_) {
    vrt.ZeroForce();
  }
  // sum forces over all vertices -- radial tether attraction/repulsion
  for (auto &&vrt : vrts_) {
    for (int i_neighb{0}; i_neighb < vrt.n_neighbs_; i_neighb++) {
      // vrt.neighb_int_[i_neighb] = true;
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
      // free movement within a certain range
      if (rmag <= l_c0_ and rmag >= l_c1_) {
        // printf("%g < %g < %g\n", l_c1_, rmag, l_c0_);
        continue;
      }
      double fmag{0.0};
      // attraction
      if (rmag > l_c0_) {
        fmag = (1 / (l_max_ - rmag)) *
               (SQR(1.0 / (rmag - l_c1_)) - (1.0 / (l_max_ - rmag))) *
               exp(1.0 / (rmag - l_c1_));
      }
      // repulsion
      if (rmag < l_c1_) {
        fmag = (1.0 / (rmag - l_min_)) *
               ((1.0 / (rmag - l_min_)) - SQR(1.0 / (l_c0_ - rmag))) *
               exp(1.0 / (l_c0_ - rmag));
      }
      fmag *= kappa_B_;
      double f[3];
      for (int i{0}; i < 3; i++) {
        rhat[i] = rhat[i] / rmag;
        f[i] = fmag * rhat[i];
      }
      vrt.AddForce(f);
      // neighb->SubForce(f);
    }
  }
  // sum forces over all vertices -- discrete bending forces
  std::vector<double> observed_thetas;
  std::vector<int> observed_counts;
  for (auto &&vrt : vrts_) {
    //scan over all neighbs
    // first we need to construct weights (sums over all neighbs)
    double sum_num{0.0};
    double sum_den{0.0};
    for (int i_neighb{0}; i_neighb < vrt.n_neighbs_; i_neighb++) {
      // this assumes neighbors are ordered in a ring
      // (this assumption is incorrect lol)
      int i_plus{i_neighb == (vrt.n_neighbs_ - 1) ? 0 : i_neighb + 1};
      int i_minus{i_neighb == 0 ? vrt.n_neighbs_ - 1 : i_neighb - 1};
      Vertex *neighb{vrt.neighbs_[i_neighb]};
      Vertex *neighb_plus{vrt.neighbs_[i_plus]};
      Vertex *neighb_minus{vrt.neighbs_[i_minus]};
      double r_ij[3];
      double r_ij_plus[3];
      double r_ij_minus[3];
      double r_jj_plus[3];
      double r_jj_minus[3];
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        r_ij[i_dim] = vrt.pos_[i_dim] - neighb->pos_[i_dim];
        r_ij_plus[i_dim] = vrt.pos_[i_dim] - neighb_plus->pos_[i_dim];
        r_ij_minus[i_dim] = vrt.pos_[i_dim] - neighb_minus->pos_[i_dim];
        r_jj_plus[i_dim] = neighb->pos_[i_dim] - neighb_plus->pos_[i_dim];
        r_jj_minus[i_dim] = neighb->pos_[i_dim] - neighb_minus->pos_[i_dim];
      }
      // the edge l_ij = r_ij connects two triangles, get unit normal of each
      double n_a[3];
      double n_b[3];
      cross_product(r_ij_minus, r_ij, n_a, 3);
      cross_product(r_ij, r_ij_plus, n_b, 3);
      // normalize these jonnies
      double l_ij{0.0};
      double l_ij_plus{0.0};
      double l_ij_minus{0.0};
      double l_jj_plus{0.0};
      double l_jj_minus{0.0};
      double norm_a{0.0};
      double norm_b{0.0};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        l_ij += SQR(r_ij[i_dim]);
        l_ij_plus += SQR(r_ij_plus[i_dim]);
        l_ij_minus += SQR(r_ij_minus[i_dim]);
        l_jj_plus += SQR(r_jj_plus[i_dim]);
        l_jj_minus += SQR(r_jj_minus[i_dim]);
        norm_a += SQR(n_a[i_dim]);
        norm_b += SQR(n_b[i_dim]);
      }
      l_ij = sqrt(l_ij);
      l_ij_plus = sqrt(l_ij_plus);
      l_ij_minus = sqrt(l_ij_minus);
      l_jj_plus = sqrt(l_jj_plus);
      l_jj_minus = sqrt(l_jj_minus);
      norm_a = sqrt(norm_a);
      norm_b = sqrt(norm_b);
      double rhat[3];
      double nhat_a[3];
      double nhat_b[3];
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        rhat[i_dim] = r_ij[i_dim] / l_ij;
        nhat_a[i_dim] = n_a[i_dim] / norm_a;
        nhat_b[i_dim] = n_b[i_dim] / norm_b;
      }
      // calculate individual variables used in force expression
      // printf("%g\n", l_ij);
      double n_cross[3];
      cross_product(n_a, n_b, n_cross, 3);
      // double theta{acos(dot_product(3, nhat_a, nhat_b) / (norm_a * norm_b))};
      // double phi_ij{M_PI - atan2(dot_product(3, r_ij, nhat_cross),
      //                            dot_product(3, nhat_a, nhat_b))};
      // SF TODO validate this; should be good tho
      double phi_ij{atan2(dot_product(3, r_ij, n_cross),
                          l_ij * dot_product(3, n_a, n_b))};
      double alpha{acos((SQR(l_ij) + SQR(l_jj_minus) - SQR(l_ij_minus)) /
                        (2 * l_ij * l_jj_minus))};
      double beta{acos((SQR(l_ij) + SQR(l_jj_plus) - SQR(l_ij_plus)) /
                       (2 * l_ij * l_jj_plus))};
      // printf("alpha: %g | beta %g\n", alpha * 180 / M_PI, beta * 180 / M_PI);
      double theta_1{acos((SQR(l_ij_minus) + SQR(l_jj_minus) - SQR(l_ij)) /
                          (2 * l_ij_minus * l_jj_minus))};
      double theta_2{acos((SQR(l_ij_plus) + SQR(l_jj_plus) - SQR(l_ij)) /
                          (2 * l_ij_plus * l_jj_plus))};

      sum_num += l_ij * phi_ij;
      sum_den += SQR(l_ij) *
                 (cos(theta_1) / sin(theta_1) + cos(theta_2) / sin(theta_2));
      // printf("theta 1: %g | theta 2: %g\n", theta_1 * 180 / M_PI,
      //        theta_2 * 180 / M_PI);
      // double leftover_1{acos((SQR(l_ij_minus) + SQR(l_ij) - SQR(l_jj_minus)) /
      //                        (2 * l_ij_minus * l_ij))};
      // double leftover_2{acos((SQR(l_ij_plus) + SQR(l_ij) - SQR(l_jj_plus)) /
      //                        (2 * l_ij_plus * l_ij))};
      // printf("%g vs %g | %g vs %g\n", leftover_1, M_PI - alpha - theta_1,
      //        leftover_2, M_PI - beta - theta_2);
      // printf("\n");
      // if (observed_thetas.empty()) {
      //   observed_thetas.push_back(theta);
      //   observed_counts.push_back(1);
      // } else {
      //   bool already_seen{false};
      //   double epsilon{1e-6}; // precision level we care about
      //   for (int i_entry{0}; i_entry < observed_thetas.size(); i_entry++) {
      //     double entry{observed_thetas[i_entry]};
      //     if (std::fabs(theta - entry) < epsilon) {
      //       already_seen = true;
      //       observed_counts[i_entry]++;
      //     }
      //   }
      //   if (!already_seen) {
      //     observed_thetas.push_back(theta);
      //     observed_counts.push_back(1);
      //   }
      // }
      // printf("(%i- vs %i-faced vertex)\n", vrt.n_neighbs_, neighb->n_neighbs_);
      // printf("theta = %g\n", theta * 180.0 / M_PI);
      // printf("phi_ij = %g deg\n\n", phi_ij * 180.0 / M_PI);
    }
    // then we can use these weights to find vector forces
    double coeff{sum_num / sum_den};
    double f_bend[3] = {0.0, 0.0, 0.0};
    for (int i_neighb{0}; i_neighb < vrt.n_neighbs_; i_neighb++) {
      // this assumes neighbors are ordered in a ring
      // (this assumption is incorrect lol)
      int i_plus{i_neighb == (vrt.n_neighbs_ - 1) ? 0 : i_neighb + 1};
      int i_minus{i_neighb == 0 ? vrt.n_neighbs_ - 1 : i_neighb - 1};
      Vertex *neighb{vrt.neighbs_[i_neighb]};
      Vertex *neighb_plus{vrt.neighbs_[i_plus]};
      Vertex *neighb_minus{vrt.neighbs_[i_minus]};
      double r_ij[3];
      double r_ij_plus[3];
      double r_ij_minus[3];
      double r_jj_plus[3];
      double r_jj_minus[3];
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        r_ij[i_dim] = vrt.pos_[i_dim] - neighb->pos_[i_dim];
        r_ij_plus[i_dim] = vrt.pos_[i_dim] - neighb_plus->pos_[i_dim];
        r_ij_minus[i_dim] = vrt.pos_[i_dim] - neighb_minus->pos_[i_dim];
        r_jj_plus[i_dim] = neighb->pos_[i_dim] - neighb_plus->pos_[i_dim];
        r_jj_minus[i_dim] = neighb->pos_[i_dim] - neighb_minus->pos_[i_dim];
      }
      // the edge l_ij = r_ij connects two triangles, get unit normal of each
      double n_a[3];
      double n_b[3];
      cross_product(r_ij_minus, r_ij, n_a, 3);
      cross_product(r_ij, r_ij_plus, n_b, 3);
      // normalize these jonnies
      double l_ij{0.0};
      double l_ij_plus{0.0};
      double l_ij_minus{0.0};
      double l_jj_plus{0.0};
      double l_jj_minus{0.0};
      double norm_a{0.0};
      double norm_b{0.0};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        l_ij += SQR(r_ij[i_dim]);
        l_ij_plus += SQR(r_ij_plus[i_dim]);
        l_ij_minus += SQR(r_ij_minus[i_dim]);
        l_jj_plus += SQR(r_jj_plus[i_dim]);
        l_jj_minus += SQR(r_jj_minus[i_dim]);
        norm_a += SQR(n_a[i_dim]);
        norm_b += SQR(n_b[i_dim]);
      }
      l_ij = sqrt(l_ij);
      l_ij_plus = sqrt(l_ij_plus);
      l_ij_minus = sqrt(l_ij_minus);
      l_jj_plus = sqrt(l_jj_plus);
      l_jj_minus = sqrt(l_jj_minus);
      norm_a = sqrt(norm_a);
      norm_b = sqrt(norm_b);
      double rhat[3];
      double nhat_a[3];
      double nhat_b[3];
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        rhat[i_dim] = r_ij[i_dim] / l_ij;
        nhat_a[i_dim] = n_a[i_dim] / norm_a;
        nhat_b[i_dim] = n_b[i_dim] / norm_b;
      }
      // calculate individual variables used in force expression
      // printf("%g\n", l_ij);
      double n_cross[3];
      cross_product(n_a, n_b, n_cross, 3);
      // double theta{acos(dot_product(3, nhat_a, nhat_b) / (norm_a * norm_b))};
      // double phi_ij{M_PI - atan2(dot_product(3, r_ij, nhat_cross),
      //                            dot_product(3, nhat_a, nhat_b))};
      // SF TODO validate this; should be good tho
      double phi_ij{atan2(dot_product(3, r_ij, n_cross),
                          l_ij * dot_product(3, n_a, n_b))};
      double alpha{acos((SQR(l_ij) + SQR(l_jj_minus) - SQR(l_ij_minus)) /
                        (2 * l_ij * l_jj_minus))};
      double beta{acos((SQR(l_ij) + SQR(l_jj_plus) - SQR(l_ij_plus)) /
                       (2 * l_ij * l_jj_plus))};
      // printf("alpha: %g | beta %g\n", alpha * 180 / M_PI, beta * 180 / M_PI);
      double theta_1{acos((SQR(l_ij_minus) + SQR(l_jj_minus) - SQR(l_ij)) /
                          (2 * l_ij_minus * l_jj_minus))};
      double theta_2{acos((SQR(l_ij_plus) + SQR(l_jj_plus) - SQR(l_ij)) /
                          (2 * l_ij_plus * l_jj_plus))};
      double T{(cos(theta_1) / sin(theta_1)) + (cos(theta_2) / sin(theta_2))};
      // square of sums, not sum of
      // radial component of force first
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        f_bend[i_dim] += 2 * coeff * phi_ij * rhat[i_dim];
        f_bend[i_dim] -= 2 * SQR(coeff) * l_ij * T;
      }
      // normal components
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        f_bend[i_dim] += 2 * coeff * cos(alpha) / sin(alpha) * nhat_a[i_dim];
        f_bend[i_dim] += 2 * coeff * cos(beta) / sin(beta) * nhat_b[i_dim];
      }
      // perp components
      double perp_a[3];
      double perp_b[3];
      cross_product(nhat_a, r_ij_minus, perp_a, 3);
      cross_product(nhat_b, r_ij_plus, perp_b, 3);
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        f_bend[i_dim] += SQR(coeff) * SQR(l_ij) /
                         (SQR(l_ij_minus) * SQR(sin(theta_1))) * perp_a[i_dim];
        f_bend[i_dim] -= SQR(coeff) * SQR(l_ij) /
                         (SQR(l_ij_plus) * SQR(sin(theta_2))) * perp_b[i_dim];
      }
    }
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      f_bend[i_dim] *= kappa_;
    }
    // finally, add the force
    vrt.AddForce(f_bend);
  }
  // for (int i_entry{0}; i_entry < observed_thetas.size(); i_entry++) {
  //   printf("angle %.4g observed %i times\n",
  //          observed_thetas[i_entry] * 180 / M_PI, observed_counts[i_entry]);
  // }
  // printf("\n");
  // apply displacements
  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    //Expected diffusion length of the crosslink in the solution in delta
    double sigma_d = sqrt(2.0 * 0.005 * 0.0131);
    //Distance actually diffused
    double dr[3] = {0.0, 0.0, 0.0};
    //Previous and final location of the crosslink
    double const *const r_prev = vrts_[i_vrt].GetPosition();
    double r_final[3];
    //Update position of the crosslink
    for (int i = 0; i < 3; ++i) {
      dr[i] = rng_.RandomNormal(sigma_d);
      double f = vrts_[i_vrt].GetForce()[i];
      // if (f > 100) {
      //   printf("%g\n", f);
      // }
      // SF TODO add gamma and timestep for them proper fiziks
      r_final[i] = r_prev[i] + dr[i] + 0.00005 * f;
    }
    vrts_[i_vrt].SetPos(r_final);
  }
}