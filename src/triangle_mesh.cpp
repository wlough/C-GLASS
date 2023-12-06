#include <cglass/triangle_mesh.hpp>

// todo lol fix seed

void TriMesh::Init(system_parameters *params) {
  // TriMesh::TriMesh() : rng_(0) {

  // SF temp before integrating with output_manager
  std::string filename{params->run_name + "_membrane_forces.file"};
  forces_ = fopen(filename.c_str(), "w");
  if (forces_ == nullptr) {
    printf("aw hell naw\n");
    exit(1);
  } else {
    file_open_ = true;
  }

  params_ = params;
  long seed{params == nullptr ? 0 : params_->seed};
  rng_ = new RNG(seed);

  r_sys_ = params->system_radius;
  kappa_B_ = params->mesh_kB;
  kappa_ = params->mesh_k;
  kappa_l_ = params->mesh_kl;
  gamma_ = params->node_gamma;

  size_t n_faces{size_t(20 * std::pow(4, params_->n_subdivisions))};
  printf("%zu faces expected\n", n_faces);

  tris_.reserve(n_faces);
  vrts_.reserve(n_faces);
  MakeIcosahedron(); // 20 triangle faces initially
  // 80 -> 320 -> 1,280 -> 5,120 -> 20,480 -> 81,920 -> 327,680 -> 1,310,720 faces
  printf("%zu VRTS TOTAL\n", vrts_.size());
  for (int i_divide{0}; i_divide < params_->n_subdivisions; i_divide++) {
    DivideFaces(); // 80 triangles
    printf("%zu VRTS TOTAL\n", vrts_.size());
  }
  ProjectToUnitSphere();
  for (auto &&vrt : vrts_) {
    vrt.neighbs_.resize(6);
  }
  // apparently the above shuffle destroys triangle ptrs??
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

  // SF TODO sloppy hack; fix
  double l_sum{0.0};
  double A_sum{0.0};
  size_t n_entries{0};
  size_t n_good{0};
  size_t n_bad{0};
  for (auto const &tri : tris_) {
    double l1{0.0};
    double l2{0.0};
    double l3{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      l1 += SQR(tri.vrts_[0]->pos_[i_dim] - tri.vrts_[1]->pos_[i_dim]);
      l2 += SQR(tri.vrts_[0]->pos_[i_dim] - tri.vrts_[2]->pos_[i_dim]);
      l3 += SQR(tri.vrts_[1]->pos_[i_dim] - tri.vrts_[2]->pos_[i_dim]);
    }
    l1 = sqrt(l1);
    l2 = sqrt(l2);
    l3 = sqrt(l3);
    if (l1 < 1) {
      n_good++;
    } else {
      n_bad++;
    }
    if (l2 < 1) {
      n_good++;
    } else {
      n_bad++;
    }
    if (l3 < 1) {
      n_good++;
    } else {
      n_bad++;
    }
    l_sum += l1;
    l_sum += l2;
    l_sum += l3;
    double s{0.5 * (l1 + l2 + l3)};
    A_sum += sqrt(s * (s - l1) * (s - l2) * (s - l3));
    n_entries += 3;
  }
  l_avg_ = l_sum / n_entries;
  A_prime_ = A_sum / (n_entries / 3);
  printf("%zu good lengths; %zu bad\n", n_good, n_bad);
  // also calculate std_err of mean to check for bugs
  double var{0.0};
  for (auto const &tri : tris_) {
    double l1{0.0};
    double l2{0.0};
    double l3{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      l1 += SQR(tri.vrts_[0]->pos_[i_dim] - tri.vrts_[1]->pos_[i_dim]);
      l2 += SQR(tri.vrts_[0]->pos_[i_dim] - tri.vrts_[2]->pos_[i_dim]);
      l3 += SQR(tri.vrts_[1]->pos_[i_dim] - tri.vrts_[2]->pos_[i_dim]);
    }
    var += SQR(l_avg_ - sqrt(l1));
    var += SQR(l_avg_ - sqrt(l2));
    var += SQR(l_avg_ - sqrt(l3));
    // n_entries += 3;
  }
  var = sqrt(var / n_entries);
  printf("l_avg = %g +/- %g\n", l_avg_, var);
  printf("A_prime = %g\n", A_prime_);
  printf("A_calc = %g\n", 4 * M_PI * SQR(r_sys_) / tris_.size());
  if (var < 0.05 * l_avg_) {
    var = 0.05 * l_avg_;
  }
  // l_c0_ = l_avg_ + var;
  // l_c1_ = l_avg_ - var;
  // l_max_ = 1.25 * l_avg_;
  // l_min_ = 0.75 * l_avg_;
  // l_avg_ *= 0.8;

  l_avg_ *= 0.9;
  l_c0_ = 1.2 * l_avg_;
  l_c1_ = 0.8 * l_avg_;
  l_max_ = 1.4 * l_avg_;
  l_min_ = 0.6 * l_avg_;

  for (auto &&vrt : vrts_) {
    vrt.SetDiameter(params_->node_diameter);
  }
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
        // check if entry_neighb is also part of vrt's neighbor list
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
              printf("ya dun goofed\n");
              exit(1);
            }
            if (dot_product(3, vrt.pos_, nhat) < 0.0) {
              // printf("this aint it cuz\n");
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
              } else {
                printf("this should never happen\n");
                exit(1);
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
      new_pos[i_dim] = r_sys_ * val / norm;
    }
    vrt.SetPos(new_pos);
  }
}

void TriMesh::MakeIcosahedron() {

  printf("initializing icosahedron (20 triangles; 12 points)\n");
  // sf todo update system radius from params
  double phi = (1.0f + sqrt(5.0f)) * 0.5f; // golden ratio
  double a = r_sys_ * 1.0f;
  double b = r_sys_ * 1.0f / phi;
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
  double l_sum{0.0};
  size_t n_entries{0};

  for (auto const &tri : tris_) {
    double l1{0.0};
    double l2{0.0};
    double l3{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      l1 += SQR(tri.vrts_[0]->pos_[i_dim] - tri.vrts_[1]->pos_[i_dim]);
      l2 += SQR(tri.vrts_[0]->pos_[i_dim] - tri.vrts_[2]->pos_[i_dim]);
      l3 += SQR(tri.vrts_[1]->pos_[i_dim] - tri.vrts_[2]->pos_[i_dim]);
    }
    l1 = sqrt(l1);
    l2 = sqrt(l2);
    l3 = sqrt(l3);
    l_sum += l1;
    l_sum += l2;
    l_sum += l3;
    n_entries += 3;
  }
  printf("edge length is %g (alt)\n", l_sum / n_entries);
}

// Minimum distance from segment to a polygon
// INPUT
// (x, y, z) line segment
// tri triangle mesh
// RETURN: triangle number in mesh
// OUTPUT
// t fractional distance along segment
// (xRot, yRot) rotated coordinates of point on triangle mesh
// dist distance of closest approach
// rcontact lab frame coordinates of closest point of approach on
//     polygon
int TriMesh::SegmentToPolygon(double xStart, double yStart, double zStart,
                              double xEnd, double yEnd, double zEnd, double *t,
                              double *xRot, double *yRot, double *dist,
                              double *rcontact) {

  // Set dist to something large
  int itriang = 0;
  *dist = 1e10;
  // for (int itri = 0; itri < tri->numTriang; ++itri) {
  for (int itri{0}; itri < tris_.size(); itri++) {
    // printf("itri = %i\n", itri);
    double tt, txRot, tyRot, tdist;
    tris_[itri].MinDist_Segment(xStart, yStart, zStart, xEnd, yEnd, zEnd, &tt,
                                &txRot, &tyRot, &tdist);
    // segment_to_triangle(xStart, yStart, zStart, xEnd, yEnd, zEnd, tri, itri,
    //                     &tt, &txRot, &tyRot, &tdist);
    if (tdist < *dist) {
      *dist = tdist;
      *t = tt;
      *xRot = txRot;
      *yRot = tyRot;
      itriang = itri;
    }
  }

  // sf todo
  // undo the rotation
  double xPrime = *xRot * tris_[itriang].cosBeta_ +
                  tris_[itriang].Zrot_ * tris_[itriang].sinBeta_;
  rcontact[0] =
      xPrime * tris_[itriang].cosGamma_ + *yRot * tris_[itriang].sinGamma_;
  rcontact[1] =
      -xPrime * tris_[itriang].sinGamma_ + *yRot * tris_[itriang].cosGamma_;
  rcontact[2] = -(*xRot) * tris_[itriang].sinBeta_ +
                tris_[itriang].Zrot_ * tris_[itriang].cosBeta_;
  //std::cout << "Contact (" << rcontact[0] << ", "
  //                         << rcontact[1] << ", "
  //                         << rcontact[2] << ")\n";

  return itriang;
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

// Finds the minimum distance from a rod to a polygon described in the triangle mesh
// OUTPUTS
// rmin: minimum r vector from point on rod to point on polygon
// rminmag2: squared minimum distance
// rcontact: lab coordinate of contact point on polygon
// mu: distance along rod for contact
int TriMesh::MinDist_Sphero(double *r_1, double *s_1, double *u_1,
                            double length_1, double *rmin, double *rminmag2,
                            double *rcontact, double *mu) {

  size_t n_dim{3};
  double r1[3] = {0.0};
  double r2[3] = {0.0};
  double t, dist;
  double xRot, yRot;

  // Convert to the tirangle mesh representation
  // r1 and r2 are simply the endpoints of each filament
  for (int i = 0; i < n_dim; ++i) {
    r1[i] = r_1[i] - 0.5 * length_1 * u_1[i];
    r2[i] = r_1[i] + 0.5 * length_1 * u_1[i];
  }

  // printf("\nr1 = <%g, %g, %g>\n", r1[0], r1[1], r1[2]);
  // printf("r2 = <%g, %g, %g>\n", r2[0], r2[1], r2[2]);

  // Run the other version
  int itri = SegmentToPolygon(r1[0], r1[1], r1[2], r2[0], r2[1], r2[2], &t,
                              &xRot, &yRot, &dist, rcontact);
  // Convert to what we need
  *rminmag2 = dist * dist;
  *mu = (t - 0.5) * length_1;
  for (int i = 0; i < n_dim; ++i) {
    rmin[i] = rcontact[i] - (r_1[i] + *mu * u_1[i]);
  }
  return itri;
}

void TriMesh::Draw(std::vector<graph_struct *> &graph_array) {

  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    vrts_[i_vrt].Draw(graph_array);
  }
  graph_array.push_back(&f_mem_);
}

void TriMesh::UpdateMesh() {
  // for (auto &&tri : tris_) {
  for (int itri{0}; itri < tris_.size(); itri++) {
    Triangle *tri{&tris_[itri]};
    double aVector[3] = {0.0};
    double bVector[3] = {0.0};
    double cVector[3] = {0.0};

    // Compute rotation to a plane and rotated coordinates
    for (int i = 0; i < 3; ++i) {
      aVector[i] =
          tri->vrts_[1]->GetPosition()[i] - tri->vrts_[0]->GetPosition()[i];
      bVector[i] =
          tri->vrts_[2]->GetPosition()[i] - tri->vrts_[1]->GetPosition()[i];
    }
    cross_product(aVector, bVector, cVector, 3);
    double cLenInXY = sqrt(cVector[0] * cVector[0] + cVector[1] * cVector[1]);
    double cLength = sqrt(cVector[0] * cVector[0] + cVector[1] * cVector[1] +
                          cVector[2] * cVector[2]);

    double sinGam, cosGam, sinBet, cosBet;
    if (cLenInXY < 1.0e-10) {
      sinGam = 0.0;
      cosGam = 1.0;
    } else {
      sinGam = -cVector[1] / cLenInXY;
      cosGam = cVector[0] / cLenInXY;
    }
    sinBet = cLenInXY / cLength;
    cosBet = cVector[2] / cLength;

    tri->cosGamma_ = cosGam;
    tri->sinGamma_ = sinGam;
    tri->cosBeta_ = cosBet;
    tri->sinBeta_ = sinBet;

    // Calculate and save rotated coordiantes
    double zRotSum = 0;
    for (int ivert = 0; ivert < 3; ++ivert) {
      int iv = ivert; //tri->indVert[ivert][itri];
      tri->XYrot_[0][ivert] = (tri->vrts_[iv]->GetPosition()[0] * cosGam -
                               tri->vrts_[iv]->GetPosition()[1] * sinGam) *
                                  cosBet -
                              tri->vrts_[iv]->GetPosition()[2] * sinBet;
      tri->XYrot_[1][ivert] = tri->vrts_[iv]->GetPosition()[0] * sinGam +
                              tri->vrts_[iv]->GetPosition()[1] * cosGam;
      zRotSum = zRotSum +
                (tri->vrts_[iv]->GetPosition()[0] * cosGam -
                 tri->vrts_[iv]->GetPosition()[1] * sinGam) *
                    sinBet +
                tri->vrts_[iv]->GetPosition()[2] * cosBet;
    }
    tri->Zrot_ = zRotSum / 3.;
  }
}

void TriMesh::ApplyBoundaryForces() {

  double r_cutoff2 = 4.0;
  double sigma2 = 4.0;
  double four_epsilon = 1.0;

  for (int i_neighb{0}; i_neighb < neighbs_.size(); i_neighb++) {
    // printf("neighb %i\n", i_neighb);
    double rmin[3], r_min_mag2, rcontact[3], mu;
    // double pos[3] = neighbs_[i_neighb].GetPosition();
    Object *neighb{neighbs_[i_neighb]};
    double r[3] = {neighb->GetPosition()[0], neighb->GetPosition()[1],
                   neighb->GetPosition()[2]};
    double u[3] = {neighb->GetOrientation()[0], neighb->GetOrientation()[1],
                   neighb->GetOrientation()[2]};
    // printf("r = <%g, %g, %g>\n", r[0], r[1], r[2]);
    // printf("u = <%g, %g, %g>\n", u[0], u[1], u[2]);
    // printf("l = %g\n", neighb->GetLength());
    // SF TODO incorporate periodic space (incorporate s and h)
    int i_tri = MinDist_Sphero(r, r, u, neighb->GetLength(), rmin, &r_min_mag2,
                               rcontact, &mu);
    double rmagcalc{0.0};
    for (int i{0}; i < 3; i++) {
      f_mem_.r[i] = tris_[i_tri].GetCenterPos(i);
      rmagcalc += SQR(rmin[i]);
    }
    // printf("%g vs %g\n", rmagcalc, r_min_mag2);
    rmagcalc = sqrt(rmagcalc);
    for (int i{0}; i < 3; i++) {
      f_mem_.u[i] = rmin[i] / rmagcalc;
    }
    // printf("u_f = <%g, %g, %g>\n", f_mem_.u[0], f_mem_.u[1], f_mem_.u[2]);
    f_mem_.length = 0.0;
    f_mem_.color = 1.8 * M_PI;
    f_mem_.diameter = 1.0;
    f_mem_.draw = draw_type::fixed;
    // printf("%g @ triangle %i\n", r_min_mag2, i_tri);
    if (r_min_mag2 < r_cutoff2) {
      //std::cout << "Firing off WCA calc\n";
      // Calculate WCA potential and forces
      double rho2 = sigma2 / r_min_mag2;
      double rho6 = CUBE(rho2);
      double rho12 = SQR(rho6);

      // printf("r_min_mag = %g\n", sqrt(r_min_mag2));
      // u += four_epsilon * (rho12 - rho6) + u_shift;
      double factor = 6.0 * four_epsilon * (2.0 * rho12 - rho6) / r_min_mag2;
      // printf("factor = %g\n", factor);

      for (int i{0}; i < 3; i++) {
        f_mem_.u[i] = rmin[i] / sqrt(r_min_mag2);
      }
      f_mem_.length = factor;
      // sf todo incorporate gamma properly
      // double f_cutoff =
      //     0.1 / params_->delta *
      //     MIN(properties->bonds.gamma_par[ibond], chromosomes->gamma_t_);
      double f_cutoff = 10000;
      // Truncate the forces if necessary
      double r_min_mag = sqrt(r_min_mag2);
      if (factor * r_min_mag > f_cutoff) {
        //std::cout << "NOTE tripping fcutoff in kcmt, resetting force factor, original " << factor << " to ";
        factor = f_cutoff / r_min_mag;
        //std::cout << factor << std::endl;
        printf(" *** Force exceeded f_cutoff "
               "kinetochoremesh_mt_wca_potential_neighbors ***\n");
      }
      double f_lj[3] = {0.0};
      for (int i = 0; i < params_->n_dim; ++i) {
        f_lj[i] = factor * rmin[i];
      }
      //std::cout << "KC MT Current forces Steric:\n";
      //std::cout << "  KC[" << ikc << "] bond[" << ibond << "]\n";
      // std::cout << "  f_lj: (" << f_lj[0] << ", " << f_lj[1] << ", " << f_lj[2]
      //           << ")\n";
      //std::cout << "  fkc:   (" << f_kc[ikc][0] << ", " << f_kc[ikc][1] << ", " << f_kc[ikc][2] << ")\n";
      //std::cout << "  fbond: (" << f_bond[ibond][0] << ", " << f_bond[ibond][1] << ", " << f_bond[ibond][2] << ")\n";

      Triangle *tri{&tris_[i_tri]};
      for (auto &&vrt : tri->vrts_) {
        vrt->AddForce(f_lj);
      }
      // tri->AddForce(f_lj);
      neighb->SubForce(f_lj);
      // Add to accumulators
      // for (int i = 0; i < ndim; ++i) {
      //   f_kc[ikc][i] += f_lj[i];
      //   f_bond[ibond][i] -= f_lj[i];
      //   //ftip[ibond][i] -= f_lj[i]; // we are definitely talking to the tip
      // }

      // SF TODO incorporate this into kinetochores
      /*
      // We need to do a special check to see if we're within the tip distance for the KC-MT interaction
      // for force dependent catastrophe
      if (properties->bonds.length[ibond] -
              (mu - 0.5) * properties->bonds.length[ibond] <
          chromosomes->chromatid_mt_fc_distance_) {
        for (int i = 0; i < ndim; ++i) {
          ftip[ibond][i] -= f_lj[i];
        }
      }
      */

      // // Calculate the virial contribution
      // if (properties->control.virial_flag == 1) {
      //   for (int i = 0; i < parameters->n_dim; ++i) {
      //     for (int j = 0; j < parameters->n_dim; ++j) {
      //       virial[i][j] += rmin[i] * f_lj[j];
      //     }
      //   }
      // }

      // SF TODO incorporate torques
      /*
      // Calculate torques
      double rcontact_kc[3] = {0.0};
      double rcontact_mt[3] = {0.0};
      for (int i = 0; i < ndim; ++i) {
        //rcontact_kc[i] = kc_iter->r_[i] - rcontact[i];
        rcontact_kc[i] = chromosomes->r_[ikc][i] - rcontact[i];
        rcontact_mt[i] = mu * u_bond[ibond][i];
      }
      //std::cout << "rcontact_kc(" << rcontact_kc[0] << ", " << rcontact_kc[1] << ", " << rcontact_kc[2] << ")\n";
      //std::cout << "rcontact_mt(" << rcontact_mt[0] << ", " << rcontact_mt[1] << ", " << rcontact_mt[2] << ")\n";
      double tau[3] = {0.0};
      cross_product(rcontact_kc, f_lj, tau, 3);
      for (int i = 0; i < 3; ++i) {
        t_kc[ikc][i] -= tau[i];
      }
      cross_product(rcontact_mt, f_lj, tau, 3);
      for (int i = 0; i < 3; ++i) {
        t_bond[ibond][i] -= tau[i];
      }
      */
    }
  }
}

void TriMesh::UpdatePositions() {

  // shrink this jonny
  l_avg_ *= 0.99999;
  l_c0_ = 1.2 * l_avg_;
  l_c1_ = 0.8 * l_avg_;
  l_max_ = 1.4 * l_avg_;
  l_min_ = 0.6 * l_avg_;

  // SF TODO all calculations are currently done redundantly
  // SF TODO i.e., we do not take advtange of newton's 3rd law
  // SF TODO once validated, increase computational efficiency by addressing this

  // for (auto &&neighb : neighbs_) {
  //   printf("position is <%g, %g, %g>\n", neighb->GetPosition()[0],
  //          neighb->GetPosition()[1], neighb->GetPosition()[2]);
  // }

  // zero forces out first
  for (auto &&vrt : vrts_) {
    vrt.ZeroForce();
  }
  double f_teth_flaw{0.0};
  size_t n_teth_flaw{0};
  double f_bend_flaw{0.0};
  size_t n_bend_flaw{0};
  double f_area_flaw{0.0};
  size_t n_area_flaw{0};
  double f_teth_good{0.0};
  size_t n_teth_good{0};
  double f_bend_good{0.0};
  size_t n_bend_good{0};
  double f_area_good{0.0};
  size_t n_area_good{0};

  double l_edge_avg{0.0};
  size_t n_edges{0};
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
      l_edge_avg += rmag;
      n_edges++;
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
      if (vrt.n_neighbs_ == 5) {
        f_teth_flaw += fmag;
        n_teth_flaw++;
      } else if (vrt.n_neighbs_ == 6) {
        f_teth_good += fmag;
        n_teth_good++;
      } else {
        printf("what [1]\n");
        exit(1);
      }
    }
  }
  // sum forces over all vertices -- discrete bending forces
  // printf("avg edge length is %g\n", avg_edge_length / n_edges);
  for (auto &&vrt : vrts_) {
    vrt.ZeroSums();
    for (int i_neighb{0}; i_neighb < vrt.n_neighbs_; i_neighb++) {
      // this assumes neighbors are ordered in a ring with + oriented ccw
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
      // the edge l_ij = |r_ij| connects two triangles, get length of each
      double l_ij{0.0};
      double l_ij_plus{0.0};
      double l_ij_minus{0.0};
      double l_jj_plus{0.0};
      double l_jj_minus{0.0};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        l_ij += SQR(r_ij[i_dim]);
        l_ij_plus += SQR(r_ij_plus[i_dim]);
        l_ij_minus += SQR(r_ij_minus[i_dim]);
        l_jj_plus += SQR(r_jj_plus[i_dim]);
        l_jj_minus += SQR(r_jj_minus[i_dim]);
      }
      l_ij = sqrt(l_ij);
      l_ij_plus = sqrt(l_ij_plus);
      l_ij_minus = sqrt(l_ij_minus);
      l_jj_plus = sqrt(l_jj_plus);
      l_jj_minus = sqrt(l_jj_minus);
      double chi_minus{dot_product(3, r_ij_minus, r_jj_minus) /
                       (l_ij_minus * l_jj_minus)};
      double chi_plus{dot_product(3, r_ij_plus, r_jj_plus) /
                      (l_ij_plus * l_jj_plus)};
      double T_ij{(chi_minus / sqrt(1.0 - SQR(chi_minus))) +
                  (chi_plus / sqrt(1.0 - SQR(chi_plus)))};
      double grad_lsq[3];
      double grad_chi_plus[3];
      double grad_chi_minus[3];
      double grad_T[3];
      //some bullshit
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        grad_lsq[i_dim] = 2 * r_ij[i_dim];
        grad_chi_plus[i_dim] =
            (1.0 / (l_ij_plus * l_jj_plus)) *
            (r_jj_plus[i_dim] -
             (l_jj_plus / l_ij_plus) * chi_plus * r_ij_plus[i_dim]);
        grad_chi_minus[i_dim] =
            (1.0 / (l_ij_minus * l_jj_minus)) *
            (r_jj_minus[i_dim] -
             (l_jj_minus / l_ij_minus) * chi_minus * r_ij_minus[i_dim]);
        grad_T[i_dim] =
            grad_chi_plus[i_dim] / std::pow((1.0 - SQR(chi_plus)), 1.5) +
            grad_chi_minus[i_dim] / std::pow((1.0 - SQR(chi_minus)), 1.5);
      }
      vrt.sum_lsqT_ += SQR(l_ij) * T_ij;
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        vrt.sum_rT_[i_dim] += r_ij[i_dim] * T_ij;
        vrt.sum_del_lsqT_[i_dim] +=
            (T_ij * grad_lsq[i_dim] + SQR(l_ij) * grad_T[i_dim]);
        for (int j_dim{0}; j_dim < 3; j_dim++) {
          vrt.sum_del_rT_[i_dim][j_dim] += r_ij[j_dim] * grad_T[i_dim];
          if (i_dim == j_dim) {
            vrt.sum_del_rT_[i_dim][j_dim] += T_ij;
          }
        }
      }
    }
    // then we can use these weights to find vector forces
    double f_bend[3] = {0.0, 0.0, 0.0};
    double f_area[3] = {0.0, 0.0, 0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      f_bend[i_dim] += 2 * dot_product(3, vrt.sum_rT_, vrt.sum_rT_) /
                       SQR(vrt.sum_lsqT_) * vrt.sum_del_lsqT_[i_dim];
      f_bend[i_dim] -= 4 * dot_product(3, vrt.sum_rT_, vrt.sum_del_rT_[i_dim]) /
                       vrt.sum_lsqT_;
      f_bend[i_dim] *= kappa_;
    }
    for (int i_neighb{0}; i_neighb < vrt.n_neighbs_; i_neighb++) {
      // this assumes neighbors are ordered in a ring with + oriented ccw
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
      // the edge l_ij = |r_ij| connects two triangles, get length of each
      double n_b[3];
      cross_product(r_ij, r_ij_plus, n_b, 3);
      double l_ij{0.0};
      double l_ij_plus{0.0};
      double l_ij_minus{0.0};
      double l_jj_plus{0.0};
      double l_jj_minus{0.0};
      double norm_b{0.0};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        l_ij += SQR(r_ij[i_dim]);
        l_ij_plus += SQR(r_ij_plus[i_dim]);
        l_ij_minus += SQR(r_ij_minus[i_dim]);
        l_jj_plus += SQR(r_jj_plus[i_dim]);
        l_jj_minus += SQR(r_jj_minus[i_dim]);
        norm_b += SQR(n_b[i_dim]);
      }
      l_ij = sqrt(l_ij);
      l_ij_plus = sqrt(l_ij_plus);
      l_ij_minus = sqrt(l_ij_minus);
      l_jj_plus = sqrt(l_jj_plus);
      l_jj_minus = sqrt(l_jj_minus);
      norm_b = sqrt(norm_b);
      // force from area conservation
      // (only add contribution from fwd neighbor)
      double area_force_vec[3];
      // cross product order intentionally reversed cuz dir of r_jj+ is wrong
      double nhat_b[3];
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        nhat_b[i_dim] = n_b[i_dim] / norm_b;
      }
      // the triangle with this nhat will be composed of vrt, neighb, and neighb_fwd
      // SF TODO fix this monstrosity with improved data structure
      bool found_tri{false};
      for (int i_tri{0}; i_tri < vrt.n_tris_; i_tri++) {
        if (found_tri) {
          break;
        }
        for (int j_tri{0}; j_tri < neighb->n_tris_; j_tri++) {
          if (found_tri) {
            break;
          }
          for (int k_tri{0}; k_tri < neighb_plus->n_tris_; k_tri++) {
            if (vrt.tris_[i_tri] == neighb->tris_[j_tri] and
                neighb->tris_[j_tri] == neighb_plus->tris_[k_tri]) {
              for (int i_dim{0}; i_dim < 3; i_dim++) {
                vrt.tris_[i_tri]->nhat_[i_dim] = nhat_b[i_dim];
              }
              found_tri = true;
              break;
            }
          }
        }
      }

      cross_product(r_jj_plus, nhat_b, area_force_vec, 3);
      // get current area of triangle
      double s{0.5 * (l_ij + l_jj_plus + l_ij_plus)};
      double A{sqrt(s * (s - l_ij) * (s - l_jj_plus) * (s - l_ij_plus))};
      double f_area_mag{-0.25 * kappa_l_ * (A - A_prime_) / A_prime_};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        f_area[i_dim] += f_area_mag * area_force_vec[i_dim];
      }
    }
    // finally, add the force
    // printf("f_bend = <%g, %g, %g>\n", f_bend[0], f_bend[1], f_bend[2]);
    // printf("f_area = <%g, %g, %g>\n\n", f_area[0], f_area[1], f_area[2]);
    vrt.AddForce(f_bend);
    vrt.AddForce(f_area);
    double fmag_bend{0.0};
    double fmag_area{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      fmag_bend += SQR(f_bend[i_dim]);
      fmag_area += SQR(f_area[i_dim]);
    }
    if (vrt.n_neighbs_ == 5) {
      f_bend_flaw += sqrt(fmag_bend);
      f_area_flaw += sqrt(fmag_area);
      n_bend_flaw++;
      n_area_flaw++;
    } else if (vrt.n_neighbs_ == 6) {
      f_bend_good += sqrt(fmag_bend);
      f_area_good += sqrt(fmag_area);
      n_bend_good++;
      n_area_good++;
    } else {
      printf("what [1]\n");
      exit(1);
    }
  }
  // for (int i_entry{0}; i_entry < observed_thetas.size(); i_entry++) {
  //   printf("angle %.4g observed %i times\n",
  //          observed_thetas[i_entry] * 180 / M_PI, observed_counts[i_entry]);
  // }
  // printf("\n");
  // apply displacements
  UpdateMesh();
  ApplyBoundaryForces();

  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    //Expected diffusion length of the crosslink in the solution in delta
    double sigma{sqrt(2 * params_->delta / gamma_)}; // kbT = 1??
    //Distance actually diffused
    double dr[3] = {0.0, 0.0, 0.0};
    //Previous and final location of the crosslink
    double const *const r_prev = vrts_[i_vrt].GetPosition();
    double r_final[3];
    //Update position of the crosslink
    for (int i = 0; i < 3; ++i) {
      // you forgot the timestep you IDJEET
      double vel{vrts_[i_vrt].GetForce()[i] / gamma_};
      double noise{rng_->RandomNormal(sigma)};
      dr[i] = vel * params_->delta + noise;
      r_final[i] = r_prev[i] + dr[i];
    }
    // if (i_vrt == 1) {
    //   r_final[1] += l_avg_ / 100;
    // printf("r_final[1] = %g\n", r_final[1]);
    // }
    vrts_[i_vrt].SetPos(r_final);
  }
  double f_avg_good[3] = {f_teth_good / n_teth_good, f_bend_good / n_bend_good,
                          f_area_good / n_area_good};
  double f_avg_flaw[3] = {f_teth_flaw / n_teth_flaw, f_bend_flaw / n_bend_flaw,
                          f_area_flaw / n_area_flaw};

  // SF temp before integrating into output_manager
  if (n_datapoints < 10000) {
    int n_written_good = fwrite(f_avg_good, sizeof(double), 3, forces_);
    int n_written_flaw = fwrite(f_avg_flaw, sizeof(double), 3, forces_);
    // if (n_written != 1) {
    //   printf("%i\n", n_written);
    //   exit(1);
    // }
    // printf("counter is %i\n", n_datapoints);
    n_datapoints++;
    // exit(1);
  } else if (file_open_) {
    file_open_ = false;
    printf("done!\n");
    fclose(forces_);
    // exit(1);
  }
  // printf("writing forces\n");
}