#include <cglass/triangle_mesh.hpp>
#include <unistd.h>

void TriMesh::Init(system_parameters *params) {

  // SF temp before integrating with output_manager
  std::string force_filename{params->run_name + "_membrane_forces.file"};
  std::string vrt_filename{params->run_name + "_membrane_vrt_positions.file"};
  std::string adj_filename{params->run_name + "_membrane_vrt_adjacency.file"};
  forces_ = fopen(force_filename.c_str(), "w");
  vertices_ = fopen(vrt_filename.c_str(), "w");
  adjacency_ = fopen(adj_filename.c_str(), "w");
  if (forces_ == nullptr or vertices_ == nullptr or adjacency_ == nullptr) {
    printf("Error generating mesh data files\n");
    exit(1);
  }
  params_ = params;
  SetParameters();
  MakeIcosphere();
  UpdateNeighbors();
  CheckVertices();
  // SF TODO relocate thi
  o_.diameter = 5.0;
  o_.length = 0.0;
  o_.draw = draw_type::fixed;
}

void TriMesh::SetParameters() {

  long seed{params_ == nullptr ? 0 : params_->seed};
  rng_ = new RNG(seed);
  r_sys_ = params_->system_radius;
  kappa_B_ = params_->mesh_kB;
  kappa_ = params_->mesh_k;
  kappa_l_ = params_->mesh_kl;
  kappa_v_ = params_->mesh_kV;
  gamma_ = params_->node_gamma;
}

void TriMesh::MakeIcosphere() {

  size_t n_faces{size_t(20 * std::pow(4, params_->n_subdivisions))};
  tris_.reserve(n_faces);
  vrts_.reserve(n_faces);
  MakeIcosahedron(); // 20 triangle faces initially
  // 80 -> 320 -> 1,280 -> 5,120 -> 20,480 -> 81,920 -> 327,680 -> 1,310,720 faces
  for (int i_divide{0}; i_divide < params_->n_subdivisions; i_divide++) {
    printf("Subdivision iteration %i:\n", i_divide + 1);
    DivideFaces();
  }
  ProjectToUnitSphere();
  for (size_t i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    vrts_[i_vrt].SetID(i_vrt);
  }
}

void TriMesh::MakeIcosahedron() {

  printf("Initializing icosahedron (20 triangles; 12 vertices)\n");
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
  // add triangles
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
  // add edges
  edges_.emplace_back(&vrts_[0], &vrts_[1]);   // 1: 1
  edges_.emplace_back(&vrts_[0], &vrts_[2]);   // 2: 1
  edges_.emplace_back(&vrts_[0], &vrts_[6]);   // 6: 1
  edges_.emplace_back(&vrts_[0], &vrts_[7]);   // 7: 1
  edges_.emplace_back(&vrts_[0], &vrts_[9]);   // 9: 1
  edges_.emplace_back(&vrts_[1], &vrts_[2]);   // 2: 2
  edges_.emplace_back(&vrts_[1], &vrts_[3]);   // 3: 1
  edges_.emplace_back(&vrts_[1], &vrts_[7]);   // 7: 2
  edges_.emplace_back(&vrts_[1], &vrts_[8]);   // 8: 1
  edges_.emplace_back(&vrts_[2], &vrts_[3]);   // 3: 2
  edges_.emplace_back(&vrts_[2], &vrts_[5]);   // 5: 1
  edges_.emplace_back(&vrts_[2], &vrts_[9]);   // 9: 2
  edges_.emplace_back(&vrts_[3], &vrts_[4]);   // 4: 1
  edges_.emplace_back(&vrts_[3], &vrts_[5]);   // 5: 2
  edges_.emplace_back(&vrts_[3], &vrts_[8]);   // 8: 2
  edges_.emplace_back(&vrts_[4], &vrts_[5]);   // 5: 3
  edges_.emplace_back(&vrts_[4], &vrts_[8]);   // 8: 3
  edges_.emplace_back(&vrts_[4], &vrts_[10]);  // 10: 1
  edges_.emplace_back(&vrts_[4], &vrts_[11]);  // 11: 1
  edges_.emplace_back(&vrts_[5], &vrts_[9]);   // 9: 3
  edges_.emplace_back(&vrts_[5], &vrts_[11]);  // 11: 2
  edges_.emplace_back(&vrts_[6], &vrts_[7]);   // 7: 3
  edges_.emplace_back(&vrts_[6], &vrts_[9]);   // 9: 4
  edges_.emplace_back(&vrts_[6], &vrts_[10]);  // 10: 2
  edges_.emplace_back(&vrts_[6], &vrts_[11]);  // 11: 3
  edges_.emplace_back(&vrts_[7], &vrts_[8]);   // 8: 4
  edges_.emplace_back(&vrts_[7], &vrts_[10]);  // 10: 3
  edges_.emplace_back(&vrts_[8], &vrts_[10]);  // 10: 4
  edges_.emplace_back(&vrts_[9], &vrts_[11]);  // 11: 4
  edges_.emplace_back(&vrts_[10], &vrts_[11]); // 11: 5
}

void TriMesh::DivideFaces() {

  std::vector<Triangle> new_faces;
  std::vector<Edge> new_edges;
  size_t n_vrts_pre{vrts_.size()};
  // loop over all triangles in mesh
  for (auto &&face : tris_) {
    Vertex *vrt0{face.vrts_[0]};
    Vertex *vrt1{face.vrts_[1]};
    Vertex *vrt2{face.vrts_[2]};
    // each edge is bisected to generate a finer mesh
    double pos10[3] = {(vrt0->pos_[0] + vrt1->pos_[0]) / 2.0,
                       (vrt0->pos_[1] + vrt1->pos_[1]) / 2.0,
                       (vrt0->pos_[2] + vrt1->pos_[2]) / 2.0};
    double pos12[3] = {(vrt2->pos_[0] + vrt1->pos_[0]) / 2.0,
                       (vrt2->pos_[1] + vrt1->pos_[1]) / 2.0,
                       (vrt2->pos_[2] + vrt1->pos_[2]) / 2.0};
    double pos20[3] = {(vrt0->pos_[0] + vrt2->pos_[0]) / 2.0,
                       (vrt0->pos_[1] + vrt2->pos_[1]) / 2.0,
                       (vrt0->pos_[2] + vrt2->pos_[2]) / 2.0};
    // ensure we do not add new vertices more than once
    Vertex *vrt10{nullptr}, *vrt12{nullptr}, *vrt20{nullptr};
    bool dup10{false}, dup12{false}, dup20{false};
    bool error{false};
    for (auto &&vrt : vrts_) {
      if (vrt.pos_[0] == pos10[0] and vrt.pos_[1] == pos10[1] and
          vrt.pos_[2] == pos10[2]) {
        if (vrt10 != nullptr) {
          error = true;
        }
        dup10 = true;
        vrt10 = &vrt;
      }
      if (vrt.pos_[0] == pos12[0] and vrt.pos_[1] == pos12[1] and
          vrt.pos_[2] == pos12[2]) {
        if (vrt12 != nullptr) {
          error = true;
        }
        dup12 = true;
        vrt12 = &vrt;
      }
      if (vrt.pos_[0] == pos20[0] and vrt.pos_[1] == pos20[1] and
          vrt.pos_[2] == pos20[2]) {
        if (vrt20 != nullptr) {
          error = true;
        }
        dup20 = true;
        vrt20 = &vrt;
      }
    }
    if (error) {
      printf("Critical error in TriMesh::DivideFaces()\n");
      exit(1);
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
    // 1 face turns into 4
    new_faces.emplace_back(vrt0, vrt10, vrt20);
    new_faces.emplace_back(vrt1, vrt10, vrt12);
    new_faces.emplace_back(vrt2, vrt12, vrt20);
    new_faces.emplace_back(vrt10, vrt12, vrt20);
    // 3 edges turns into 9
    // only add outer edges if vertex is new -- otherwise will double count
    if (!dup10) {
      new_edges.emplace_back(vrt0, vrt10);
      new_edges.emplace_back(vrt1, vrt10);
    }
    if (!dup20) {
      new_edges.emplace_back(vrt0, vrt20);
      new_edges.emplace_back(vrt2, vrt20);
    }
    if (!dup12) {
      new_edges.emplace_back(vrt1, vrt12);
      new_edges.emplace_back(vrt2, vrt12);
    }
    // inner edges are always added since they cannot be double counted
    new_edges.emplace_back(vrt10, vrt20);
    new_edges.emplace_back(vrt10, vrt12);
    new_edges.emplace_back(vrt12, vrt20);
  }
  printf("  %zu -> %zu triangles\n", tris_.size(), new_faces.size());
  printf("  %zu -> %zu vertices\n", n_vrts_pre, vrts_.size());
  printf("  %zu -> %zu edges\n", edges_.size(), new_edges.size());
  // printf("  (avg edge length is now %g)\n", avg_edge_length / n_edges);
  tris_ = new_faces;
  edges_ = new_edges;
}

void TriMesh::UpdateOrigin() {

  // SF TODO include area weights
  for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
    origin_[i_dim] = 0.0;
  }
  for (auto &&tri : tris_) {
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      origin_[i_dim] += tri.GetCenterPos(i_dim) / tris_.size();
    }
  }
  o_.r[0] = origin_[0];
  o_.r[1] = origin_[1];
  o_.r[2] = origin_[2];
}

void TriMesh::UpdateNeighbors() {

  UpdateTriangles();
  // NOTE: currently only works with convex objects
  // NOTE: an update to ensure neighboring normals are aligned would fix this
  // This function is purposefully coded to be explicit yet inefficient
  // (it does not get called often if at all beyond initialization)
  for (auto &&vrt : vrts_) {
    vrt.neighbs_.resize(10);
    vrt.tris_.resize(10);
    vrt.edges_.resize(10);
  }
  // Reset storage of neighb and triangle ptrs in each vertex
  size_t i_vrt{0};
  for (auto &&vrt : vrts_) {
    vrt.i_ = i_vrt++;
    vrt.n_tris_ = 0;
    vrt.n_neighbs_ = 0;
    vrt.n_edges_ = 0;
  }
  // Update triangle ptrs stored by each vertex
  size_t i_tri{0};
  for (auto &&tri : tris_) {
    tri.i_ = i_tri++;
    tri.vrts_[0]->tris_[tri.vrts_[0]->n_tris_++] = &tri;
    tri.vrts_[1]->tris_[tri.vrts_[1]->n_tris_++] = &tri;
    tri.vrts_[2]->tris_[tri.vrts_[2]->n_tris_++] = &tri;
    for (int i_edge{0}; i_edge < 3; i_edge++) {
      tri.edges_[i_edge] = nullptr;
    }
  }
  // Update edge + neighb ptrs stored by each vertex
  size_t i_edge{0};
  for (auto &&edge : edges_) {
    edge.i_ = i_edge++;
    edge.vrts_[0]->edges_[edge.vrts_[0]->n_edges_++] = &edge;
    edge.vrts_[1]->edges_[edge.vrts_[1]->n_edges_++] = &edge;
    edge.vrts_[0]->neighbs_[edge.vrts_[0]->n_neighbs_++] = edge.vrts_[1];
    edge.vrts_[1]->neighbs_[edge.vrts_[1]->n_neighbs_++] = edge.vrts_[0];
  }

  // Update triangles and edges so that they track each other
  for (auto &&edge : edges_) {
    int n_found{0};
    // the two triangles we are a part of will be contained by each vertex
    for (int i_tri{0}; i_tri < edge.vrts_[0]->n_tris_; i_tri++) {
      Triangle *tri{edge.vrts_[0]->tris_[i_tri]};
      if (tri->Contains(edge.vrts_[1])) {
        // store triangle pointer in edge
        edge.tris_[n_found++] = tri;
        // add edge pointer to said triangle
        for (int i_edge{0}; i_edge < 3; i_edge++) {
          // printf("query %i\n", i_edge);
          if (tri->edges_[i_edge] == nullptr) {
            // printf("edge %i added\n", i_edge);
            tri->edges_[i_edge] = &edge;
            break;
          } else if (i_edge == 2) {
            printf("error? i think?\n");
            exit(1);
          }
        }
      }
    }
    if (n_found != 2) {
      printf("Error: found %i triangles that contain one edge\n", n_found);
      exit(1);
    }
  }
  // for (auto &&tri : tris_) {
  //   printf("tri %zu composed of vrts %zu, %zu, %zu\n", tri.i_, tri.vrts_[0]->i_,
  //          tri.vrts_[1]->i_, tri.vrts_[2]->i_);
  //   printf("tri %zu composed of edges %zu, %zu, %zu\n", tri.i_,
  //          tri.edges_[0]->i_, tri.edges_[1]->i_, tri.edges_[2]->i_);
  // }
  // Order triangle ptrs so that they are in consecutive order (share edges)
  for (auto &&vrt : vrts_) {
    Triangle *tris_ordered[vrt.n_tris_]{{}};
    Edge *edges_ordered[vrt.n_edges_]{{}};
    Vertex *neighbs_ordered[vrt.n_neighbs_]{{}};
    if (vrt.n_tris_ != vrt.n_edges_ or vrt.n_tris_ != vrt.n_neighbs_) {
      printf("Error in neighbor lists\n");
      exit(1);
    }
    // starting point doesn't matter, so just use whichever entry is first
    neighbs_ordered[0] = vrt.neighbs_[0];
    for (int i_edge{0}; i_edge < vrt.n_edges_; i_edge++) {
      Edge *edge{vrt.edges_[i_edge]};
      if (edge->Contains(&vrt) and edge->Contains(vrt.neighbs_[0])) {
        edges_ordered[0] = edge;
        break;
      }
    }
    if (edges_ordered[0] == nullptr) {
      printf("Error: first edge not found in UpdateNeighbors\n");
      exit(1);
    }
    // once starting point is established, loop over other edges to sort
    for (int i_entry{1}; i_entry < vrt.n_edges_; i_entry++) {
      Edge *prev_edge{edges_ordered[i_entry - 1]};
      Vertex *prev_neighb{neighbs_ordered[i_entry - 1]};
      for (int i_edge{0}; i_edge < vrt.n_edges_; i_edge++) {
        Edge *next_edge{vrt.edges_[i_edge]};
        // disregard self-comparisons (free life advice 4 ya)
        if (prev_edge == next_edge) {
          // printf("doink\n");
          continue;
        }
        Triangle *next_tri{nullptr};
        // SF TODO remove chicanery
        if (prev_edge == nullptr) {
          printf("no thx @ %i\n", i_entry);
          printf("(%i entries total)\n", vrt.n_edges_);
          vrt.SetColor(2 * M_PI, draw_type::fixed);
          vrt.SetDiameter(2);
          do_not_pass_go_ = true;
          return;
          // exit(1);
        }
        for (auto &&tri_i : prev_edge->tris_) {
          for (auto &&tri_j : next_edge->tris_) {
            if (tri_i == tri_j) {
              if (next_tri != nullptr) {
                printf("can't have more than 1 match my guy\n");
                exit(1);
              }
              next_tri = tri_i;
            }
          }
        }
        // if triangle match wasn't found, these edges are not adjacent
        if (next_tri == nullptr) {
          continue;
        }
        // if triangle was found, need to check its ring-wise direction
        Vertex *next_neighb{next_edge->GetOtherEnd(&vrt)};
        // points from prev_neighb (j_prev) to vrt (i)
        double r_ij_prev[3];
        // ponts from next_neighb (j_next) to vrt (i);
        double r_ij_next[3];
        for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
          r_ij_prev[i_dim] = vrt.pos_[i_dim] - prev_neighb->pos_[i_dim];
          r_ij_next[i_dim] = vrt.pos_[i_dim] - next_neighb->pos_[i_dim];
        }
        // for proper ccw order, r_ij_prev x r_ij_next should align with nhat
        double n_tri[3];
        cross_product(r_ij_prev, r_ij_next, n_tri, 3);
        // printf("n_tri = <%.2f, %.2f, %.2f>\n", n_tri[0], n_tri[1], n_tri[2]);
        // printf("nhat = <%.2f, %.2f, %.2f>\n", next_tri->nhat_[0],
        //        next_tri->nhat_[1], next_tri->nhat_[2]);
        if (dot_product(3, n_tri, next_tri->nhat_) > 0.0) {
          // printf("success\n");
          tris_ordered[i_entry - 1] = next_tri;
          edges_ordered[i_entry] = next_edge;
          neighbs_ordered[i_entry] = next_neighb;
          // since we are adding tris "behind" edges/neighbs, last tri must be added manually
          if (i_entry == vrt.n_edges_ - 1) {
            tris_ordered[i_entry] = next_edge->GetOtherTriangle(next_tri);
          }
        }
      }
    }
    for (int i_entry{0}; i_entry < vrt.n_neighbs_; i_entry++) {
      vrt.neighbs_[i_entry] = neighbs_ordered[i_entry];
      vrt.edges_[i_entry] = edges_ordered[i_entry];
      vrt.tris_[i_entry] = tris_ordered[i_entry];
    }
  }
}

void TriMesh::CheckVertices() {

  size_t n_flawed{0};
  size_t n_gucci{0};
  for (auto &&vrt : vrts_) {
    vrt.SetDiameter(params_->node_diameter);
    vrt.SetColor(1.5, draw_type::fixed);
    if (vrt.n_neighbs_ != 6) {
      n_flawed++;
    } else if (vrt.n_neighbs_ == 6) {
      n_gucci++;
    } else {
      printf("ummm?? in TRIANGLE MESH INIT\n");
      exit(1);
    }
  }
  printf("Final mesh statistics:\n");
  printf("  %zu vrts (%zu flawed; %zu ideal)\n", vrts_.size(), n_flawed,
         n_gucci);

  // SF TODO sloppy hack; fix
  double l_sum{0.0};
  double A_sum{0.0};
  double V_sum{0.0};
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
    double s{0.5 * (l1 + l2 + l3)};
    A_sum += sqrt(s * (s - l1) * (s - l2) * (s - l3));
    double A[3];
    double B[3];
    double C[3];
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      A[i_dim] = tri.vrts_[0]->pos_[i_dim] - origin_[i_dim];
      B[i_dim] = tri.vrts_[1]->pos_[i_dim] - origin_[i_dim];
      C[i_dim] = tri.vrts_[2]->pos_[i_dim] - origin_[i_dim];
    }
    double BxC[3];
    cross_product(B, C, BxC, 3);
    double vol_prism{std::fabs(dot_product(3, A, BxC) / 6.0)};
    V_sum += vol_prism;
    n_entries += 3;
  }
  l_avg_ = l_sum / n_entries;
  A_prime_ = A_sum / (n_entries / 3);
  V_prime_ = V_sum / (n_entries / 3);
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
  }
  var = sqrt(var / n_entries);
  printf("  l_avg = %g +/- %g\n", l_avg_, var);
  printf("  A_prime = %g\n", A_prime_);
  printf("  A_calc = %g\n", 4.0 * M_PI * SQR(r_sys_) / tris_.size());
  printf("  V_prime = %g\n", V_prime_);
  printf("  V_calc = %g\n", (4.0 / 3.0) * M_PI * CUBE(r_sys_) / tris_.size());
  printf("  V_calc_alt = %g\n", (1.0 / 3.0) * A_prime_ * r_sys_);
  if (var < 0.05 * l_avg_) {
    var = 0.05 * l_avg_;
  }
  l_c0_ = 1.2 * l_avg_;
  l_c1_ = 0.8 * l_avg_;
  // l_max_ = 1.4 * l_avg_;
  // l_min_ = 0.6 * l_avg_;
  l_max_ = 1.67 * l_avg_;
  l_min_ = 0.33 * l_avg_;
}

void TriMesh::FlipEdges() {

  for (auto &&edge : edges_) {
    edge.just_flipped = false;
  }
  bool flipparino{false};
  for (auto &&edge : edges_) {
    if (edge.just_flipped) {
      continue;
    }
    Vertex *vrt1 = edge.vrts_[0];
    Vertex *vrt2 = edge.vrts_[1];
    // Find adjacent triangles
    Triangle *left{edge.tris_[0]};
    Triangle *right{edge.tris_[1]};
    if (left == nullptr or right == nullptr) {
      printf("issue at edge %zu\n", edge.i_);
      printf("ABORT\n");
      exit(1);
    }
    Vertex *vrt3{left->GetOtherVertex(vrt1, vrt2)};
    Vertex *vrt4{right->GetOtherVertex(vrt1, vrt2)};
    if (vrt3 == nullptr or vrt4 == nullptr) {
      printf("issue 2 at edge %zu\n", edge.i_);
      printf("ABORT\n");
      exit(1);
    }
    // original triangles: 1->2->3->1 and 1->2->4->1
    // post-flip triangles: 1->3->4->1 and 2->3->4->2
    double l12{0.0}; // original edge
    double l34{0.0}; // flipped edge
    // shared edges
    double l23{0.0};
    double l31{0.0};
    double l24{0.0};
    double l41{0.0};
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      l12 += SQR(vrt1->pos_[i_dim] - vrt2->pos_[i_dim]);
      l34 += SQR(vrt3->pos_[i_dim] - vrt4->pos_[i_dim]);
      l23 += SQR(vrt2->pos_[i_dim] - vrt3->pos_[i_dim]);
      l31 += SQR(vrt3->pos_[i_dim] - vrt1->pos_[i_dim]);
      l24 += SQR(vrt2->pos_[i_dim] - vrt4->pos_[i_dim]);
      l41 += SQR(vrt4->pos_[i_dim] - vrt1->pos_[i_dim]);
    }
    l12 = sqrt(l12);
    l34 = sqrt(l34);
    l23 = sqrt(l23);
    l31 = sqrt(l31);
    l24 = sqrt(l24);
    l41 = sqrt(l41);
    // calculate energy pre-flip
    double energy_i{0.0};
    if (l12 > l_c0_) {
      energy_i += kappa_B_ * exp(1.0 / (l_c0_ - l12)) / (l_max_ - l12);
    } else if (l12 < l_c1_) {
      energy_i += kappa_B_ * exp(1.0 / (l12 - l_c1_)) / (l12 - l_min_);
    }
    double s1_i{0.5 * (l12 + l23 + l31)};
    double area1_i{sqrt(s1_i * (s1_i - l12) * (s1_i - l23) * (s1_i - l31))};
    energy_i += 0.5 * kappa_l_ * SQR(A_prime_ - area1_i) / A_prime_;
    double s2_i{0.5 * (l12 + l24 + l41)};
    double area2_i{sqrt(s2_i * (s2_i - l12) * (s2_i - l24) * (s2_i - l41))};
    energy_i += 0.5 * kappa_l_ * SQR(A_prime_ - area2_i) / A_prime_;
    // calculate energy post-flip
    double energy_f{0.0};
    bool impossible{false};
    if (l34 >= l_max_ or l34 <= l_min_) {
      impossible = true;
    } else if (l34 > l_c0_) {
      energy_f += kappa_B_ * exp(1.0 / (l_c0_ - l34)) / (l_max_ - l34);
    } else if (l34 < l_c1_) {
      energy_f += kappa_B_ * exp(1.0 / (l34 - l_c1_)) / (l34 - l_min_);
    }
    double s1_f{0.5 * (l34 + l41 + l31)};
    double area1_f{sqrt(s1_f * (s1_f - l34) * (s1_f - l41) * (s1_f - l31))};
    energy_f += 0.5 * kappa_l_ * SQR(A_prime_ - area1_f) / A_prime_;
    double s2_f{0.5 * (l34 + l24 + l23)};
    double area2_f{sqrt(s2_f * (s2_f - l34) * (s2_f - l41) * (s2_f - l31))};
    energy_f += 0.5 * kappa_l_ * SQR(A_prime_ - area2_f) / A_prime_;
    // area can become NaN if three vertices are in a line
    if (area1_f != area1_f or area2_f != area2_f) {
      impossible = true;
    }
    bool do_flip{false};
    if (not impossible) {
      if (energy_f < energy_i) {
        // printf("AUTO PASS \n");
        do_flip = true;
      } else {
        // printf("pre: %.4f | post %.4f\n", energy_i, energy_f);
        double delta{energy_f - energy_i};
        double p_flip{exp(-delta)};
        // printf("  p = %.4f\n", p_flip);
        double ran = rng_->RandomUniform();
        if (ran < p_flip) {
          do_flip = true;
        }
      }
    }
    if (do_flip) {
      // printf("FLIP edge %i\n", i_edge);
      flipparino = true;
      // flipped triangles: 3->2->4 and 3->1->4
      edges_[edge.i_].vrts_[0] = vrt3;
      edges_[edge.i_].vrts_[1] = vrt4;
      tris_[left->i_].vrts_[0] = vrt2;
      tris_[left->i_].vrts_[1] = vrt3;
      tris_[left->i_].vrts_[2] = vrt4;
      tris_[right->i_].vrts_[0] = vrt1;
      tris_[right->i_].vrts_[1] = vrt3;
      tris_[right->i_].vrts_[2] = vrt4;
      // need to mark all edges as just flipped; otherwise can get wonky triangles
      for (int i_edge{0}; i_edge < 3; i_edge++) {
        left->edges_[i_edge]->just_flipped = true;
        right->edges_[i_edge]->just_flipped = true;
      }
      vrt1->SetColor(0.0, draw_type::fixed);
      vrt2->SetColor(0.0, draw_type::fixed);
      vrt3->SetColor(0.0, draw_type::fixed);
      vrt4->SetColor(0.0, draw_type::fixed);
    }
  }
  if (flipparino) {
    UpdateTriangles();
    UpdateNeighbors();
  }
}

void TriMesh::UpdateTriangles() {

  UpdateOrigin();
  for (auto &&tri : tris_) {
    double r1[3], r2[3]; // two edges from this triangle
    double r_origin[3];  // points from origin to center of triangle
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      r1[i_dim] = tri.vrts_[0]->pos_[i_dim] - tri.vrts_[1]->pos_[i_dim];
      r2[i_dim] = tri.vrts_[1]->pos_[i_dim] - tri.vrts_[2]->pos_[i_dim];
      r_origin[i_dim] = tri.GetCenterPos(i_dim) - origin_[i_dim];
    }
    double n[3];
    cross_product(r1, r2, n, 3);
    if (dot_product(3, n, r_origin) < 0.0) {
      // cross_product(r2, r1, n, 3);
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        n[i_dim] *= -1;
      }
      if (dot_product(3, n, r_origin) < 0.0) {
        printf("what the fuc\n");
      }
    }
    double n_mag{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      n_mag += SQR(n[i_dim]);
    }
    n_mag = sqrt(n_mag);
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      tri.nhat_[i_dim] = n[i_dim] / n_mag;
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

void TriMesh::Draw(std::vector<graph_struct *> &graph_array) {

  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    vrts_[i_vrt].Draw(graph_array);
  }
  for (int i{0}; i < f_mem_.size(); i++) {
    graph_array.push_back(&f_mem_[i]);
  }
  graph_array.push_back(&o_);
}

void TriMesh::UpdateMesh() {
  for (auto &&vrt : vrts_) {
    vrt.ZeroForce();
  }
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

void TriMesh::ApplyMembraneForces() {

  // SF TODO all calculations are currently done redundantly
  // SF TODO i.e., we do not take advtange of newton's 3rd law
  // SF TODO once validated, increase computational efficiency by addressing this

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

  // SAF TODO -- make loop over edges (with pre-updated r vec+mag)
  // TETHER FORCE
  for (auto &&vrt : vrts_) {
    for (int i_neighb{0}; i_neighb < vrt.n_neighbs_; i_neighb++) {
      Vertex *neighb{vrt.neighbs_[i_neighb]};
      if (neighb == nullptr) {
        printf("what in TriMesh::UpdatePosition()\n");
        exit(1);
      }
      double r_ij[3];
      double rmag{0.0};
      for (int i{0}; i < 3; i++) {
        r_ij[i] = vrt.GetPosition()[i] - neighb->GetPosition()[i];
        rmag += SQR(r_ij[i]);
      }
      rmag = sqrt(rmag);
      // free movement within a certain range
      if (rmag <= l_c0_ and rmag >= l_c1_) {
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
        f[i] = fmag * r_ij[i] / rmag;
      }
      vrt.AddForce(f);
      if (vrt.n_neighbs_ == 5) {
        f_teth_flaw += fmag;
        n_teth_flaw++;
      } else if (vrt.n_neighbs_ == 6) {
        f_teth_good += fmag;
        n_teth_good++;
      } else {
        // printf("what [1] - %i\n", vrt.n_neighbs_);
        // exit(1);
      }
    }
  }
  // SF TODO think about better data structure -- vrts seems best tho
  // DISCRETE BENDING FORCES
  for (auto &&vrt : vrts_) {
    double sum_lsqT{0.0};        // scalar
    double sum_del_lsqT[3]{{}};  // vec; scalar in each dim
    double sum_rT[3]{{}};        // vec; scalr in each dim
    double sum_del_rT[3][3]{{}}; // tensor; vec. in each dim
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
      // c.f. Appendix A of Guckenberger et al. Comp. Phys. Comm (2016)
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
      sum_lsqT += SQR(l_ij) * T_ij;
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        sum_rT[i_dim] += r_ij[i_dim] * T_ij;
        sum_del_lsqT[i_dim] +=
            (T_ij * grad_lsq[i_dim] + SQR(l_ij) * grad_T[i_dim]);
        for (int j_dim{0}; j_dim < 3; j_dim++) {
          sum_del_rT[i_dim][j_dim] += r_ij[j_dim] * grad_T[i_dim];
          if (i_dim == j_dim) {
            sum_del_rT[i_dim][j_dim] += T_ij;
          }
        }
      }
    }
    // then we can use these weights to find vector forces
    double f_bend[3]{{}};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      f_bend[i_dim] += 2 * dot_product(3, sum_rT, sum_rT) / SQR(sum_lsqT) *
                       sum_del_lsqT[i_dim];
      f_bend[i_dim] -= 4 * dot_product(3, sum_rT, sum_del_rT[i_dim]) / sum_lsqT;
      f_bend[i_dim] *= kappa_;
    }
    vrt.AddForce(f_bend);
    double fmag_bend{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      fmag_bend += SQR(f_bend[i_dim]);
    }
    if (vrt.n_neighbs_ == 5) {
      f_bend_flaw += sqrt(fmag_bend);
      n_bend_flaw++;
    } else if (vrt.n_neighbs_ == 6) {
      f_bend_good += sqrt(fmag_bend);
      n_bend_good++;
    } else {
      // printf("what [2] - %i\n", vrt.n_neighbs_);
      // exit(1);
    }
  }
  // SF TODO loop over pre-updated triangles?
  // AREA AND VOLUME CONSERVATION FORCE
  for (auto &&vrt : vrts_) {
    double f_area[3]{{}};
    double f_vol[3]{{}};
    for (int i_neighb{0}; i_neighb < vrt.n_neighbs_; i_neighb++) {
      // this assumes neighbors are ordered in a ring with + oriented ccw
      int i_plus{i_neighb == (vrt.n_neighbs_ - 1) ? 0 : i_neighb + 1};
      Vertex *neighb{vrt.neighbs_[i_neighb]};
      Vertex *neighb_plus{vrt.neighbs_[i_plus]};
      Triangle *tri_plus{vrt.tris_[i_neighb]}; // from neighb list construction
      double r_ij[3];
      double r_ij_plus[3];
      double r_jj_plus[3];
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        r_ij[i_dim] = vrt.pos_[i_dim] - neighb->pos_[i_dim];
        r_ij_plus[i_dim] = vrt.pos_[i_dim] - neighb_plus->pos_[i_dim];
        r_jj_plus[i_dim] = neighb->pos_[i_dim] - neighb_plus->pos_[i_dim];
      }
      // the edge l_ij = |r_ij| connects two triangles, get length of each
      // double n_b[3];
      // cross_product(r_ij, r_ij_plus, n_b, 3);
      double l_ij{0.0};
      double l_ij_plus{0.0};
      double l_jj_plus{0.0};
      // double norm_b{0.0};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        l_ij += SQR(r_ij[i_dim]);
        l_ij_plus += SQR(r_ij_plus[i_dim]);
        l_jj_plus += SQR(r_jj_plus[i_dim]);
        // norm_b += SQR(n_b[i_dim]);
      }
      l_ij = sqrt(l_ij);
      l_ij_plus = sqrt(l_ij_plus);
      l_jj_plus = sqrt(l_jj_plus);
      // norm_b = sqrt(norm_b);
      // force from area conservation
      // (only add contribution from fwd neighbor)
      double f_area_vec[3];
      // cross product order intentionally reversed cuz dir of r_jj+ is wrong
      // double nhat_b[3];
      // for (int i_dim{0}; i_dim < 3; i_dim++) {
      //   nhat_b[i_dim] = n_b[i_dim] / norm_b;
      // }
      cross_product(r_jj_plus, tri_plus->nhat_, f_area_vec, 3);
      // get current area of triangle
      double s{0.5 * (l_ij + l_jj_plus + l_ij_plus)};
      double area{sqrt(s * (s - l_ij) * (s - l_jj_plus) * (s - l_ij_plus))};
      double f_area_mag{-0.25 * kappa_l_ * (area - A_prime_) / A_prime_};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        f_area[i_dim] += f_area_mag * f_area_vec[i_dim];
      }
      // neighb vertices for volume force will be same as area (i.e., fwd only)
      // need area of triangle that connects these two vertices w/ polygon center
      // (opposite edge is l_jj_plus)
      // to keep consistent naming, r_oj_plus points from j+ to origin
      // likewise, r_oj points from j to origin
      double r_oj[3];
      double r_oj_plus[3];
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        r_oj[i_dim] = origin_[i_dim] - neighb->GetPosition()[i_dim];
        r_oj_plus[i_dim] = origin_[i_dim] - neighb_plus->GetPosition()[i_dim];
      }
      double l_oj{0.0};
      double l_oj_plus{0.0};
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        l_oj += SQR(r_oj[i_dim]);
        l_oj_plus += SQR(r_oj_plus[i_dim]);
      }
      l_oj = sqrt(l_oj);
      l_oj_plus = sqrt(l_oj_plus);
      double s_opp{0.5 * (l_jj_plus + l_oj + l_oj_plus)};
      double area_opp{sqrt(s_opp * (s_opp - l_jj_plus) * (s_opp - l_oj) *
                           (s_opp - l_oj_plus))};
      double A[3];
      double B[3];
      double C[3];
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        A[i_dim] = vrt.pos_[i_dim] - origin_[i_dim];
        B[i_dim] = neighb->pos_[i_dim] - origin_[i_dim];
        C[i_dim] = neighb_plus->pos_[i_dim] - origin_[i_dim];
      }
      double BxC[3];
      cross_product(B, C, BxC, 3);
      double V{std::fabs(dot_product(3, A, BxC) / 6.0)};
      double f_vol_mag{(-1.0 / 3.0) * kappa_v_ * (V - V_prime_) * area_opp /
                       V_prime_};
      if (f_vol_mag != f_vol_mag) {
        vrt.SetColor(0.0, draw_type::fixed);
        neighb_plus->SetColor(1.0, draw_type::fixed);
        printf("V: %g\n", V);
        printf("V_prime: %g\n", V_prime_);
        printf("area_opp: %g\n", area_opp);
        do_not_pass_go_ = true;
        double zero_force[3] = {0.0, 0.0, 0.0};
        for (auto &&vrt : vrts_) {
          vrt.SetForce(zero_force);
        }
        printf("holla\n");
        return;
      }
      // just the inward pointing normal of the face opposite to vrt in prism
      double f_vol_vec[3];
      cross_product(r_oj_plus, r_jj_plus, f_vol_vec, 3);
      double f_vol_vec_norm{0.0};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        f_vol_vec_norm += SQR(f_vol_vec[i_dim]);
      }
      f_vol_vec_norm = sqrt(f_vol_vec_norm);
      if (f_vol_vec_norm == 0.0) {
        printf("r_oj+ = <%g, %g, %g>\n", r_oj_plus[0], r_oj_plus[1],
               r_oj_plus[2]);
        printf("r_jj+ = <%g, %g, %g>\n", r_jj_plus[0], r_jj_plus[1],
               r_jj_plus[2]);
        printf("NOPE\n");
        exit(1);
      }
      // vec is going -NaN, meaning the norm is zero?? WHY
      // should be in-line with r_ij
      // if (dot_product(3, r_ij, f_vol_vec) < 0.0) {
      //   printf("no way bro\n");
      //   vrt.SetColor(0.0, draw_type::fixed);
      //   neighb_plus->SetColor(1.0, draw_type::fixed);
      //   do_not_pass_go_ = true;
      //   return;
      //   exit(1);
      // }
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        f_vol[i_dim] += f_vol_mag * f_vol_vec[i_dim] / f_vol_vec_norm;
        // f_vol[i_dim] += f_vol_vec[i_dim] / f_vol_vec_norm;
        if (f_vol[i_dim] != f_vol[i_dim]) {
          printf("why - %g\n", f_vol[i_dim]);
          exit(1);
        }
      }
    }
    vrt.AddForce(f_area);
    vrt.AddForce(f_vol);
    /*
    double fmag_area{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      fmag_area += SQR(f_area[i_dim]);
    }
    if (vrt.n_neighbs_ == 5) {
      f_area_flaw += sqrt(fmag_area);
      n_area_flaw++;
    } else if (vrt.n_neighbs_ == 6) {
      f_area_good += sqrt(fmag_area);
      n_area_good++;
    } else {
      // printf("what [3] - %i\n", vrt.n_neighbs_);
      // exit(1);
    }
    */
  }
  // SF TODO add volume forces?
  // f_avg_ideal_[0] = f_teth_good / n_teth_good;
  // f_avg_ideal_[1] = f_bend_good / n_bend_good;
  // f_avg_ideal_[2] = f_area_good / n_area_good;
  // f_avg_flawed_[0] = f_teth_flaw / n_teth_flaw;
  // f_avg_flawed_[1] = f_bend_flaw / n_bend_flaw;
  // f_avg_flawed_[2] = f_area_flaw / n_area_flaw;
}

void TriMesh::ApplyBoundaryForces() {

  // double r_cutoff2 = SQR(pow(2, 1.0 / 6.0) * 0.5);
  // double sigma2 = 0.25;
  double r_cutoff2 = 1.5;
  double sigma2 = 1.5;
  double four_epsilon = 4.0;

  f_mem_.resize(neighbs_.size());
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
    int i_tri = mindist_.SpheroPolygon(this, r, r, u, neighb->GetLength(), rmin,
                                       &r_min_mag2, rcontact, &mu);
    double rmagcalc{0.0};
    for (int i{0}; i < 3; i++) {
      f_mem_[i_neighb].r[i] = tris_[i_tri].GetCenterPos(i);
      rmagcalc += SQR(rmin[i]);
    }
    // printf("%g vs %g\n", rmagcalc, r_min_mag2);
    rmagcalc = sqrt(rmagcalc);
    for (int i{0}; i < 3; i++) {
      f_mem_[i_neighb].u[i] = rmin[i] / rmagcalc;
    }
    // printf("u_f = <%g, %g, %g>\n", f_mem_.u[0], f_mem_.u[1], f_mem_.u[2]);
    f_mem_[i_neighb].length = 0.0;
    f_mem_[i_neighb].color = 1.8 * M_PI;
    f_mem_[i_neighb].diameter = 1.0;
    f_mem_[i_neighb].draw = draw_type::fixed;
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
        f_mem_[i_neighb].u[i] = rmin[i] / sqrt(r_min_mag2);
      }
      f_mem_[i_neighb].length = factor;
      // sf todo incorporate gamma properly
      // double f_cutoff =
      //     0.1 / params_->delta *
      //     MIN(properties->bonds.gamma_par[ibond], chromosomes->gamma_t_);
      double f_cutoff = 250;
      // Truncate the forces if necessary
      double r_min_mag = sqrt(r_min_mag2);
      if (factor * r_min_mag > f_cutoff) {
        //std::cout << "NOTE tripping fcutoff in kcmt, resetting force factor, original " << factor << " to ";
        std::cout << factor << std::endl;
        printf(" *** Force exceeded f_cutoff "
               "kinetochoremesh_mt_wca_potential_neighbors ***\n");
        factor = f_cutoff / r_min_mag;
      }
      double f_lj[3] = {0.0};
      for (int i = 0; i < params_->n_dim; ++i) {
        f_lj[i] = factor * rmin[i];
      }
      Triangle *tri{&tris_[i_tri]};
      for (auto &&vrt : tri->vrts_) {
        vrt->AddForce(f_lj);
      }
      // tri->AddForce(f_lj);
      neighb->SubForce(f_lj);

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

void TriMesh::WriteOutputs() {
  // SF temp before integrating into output_manager
  if (params_->i_step % params_->mesh_steps_per_datapoint != 0) {
    return;
  }
  if (i_datapoint_ < params_->mesh_datapoints) {
    // write different membrane forces
    // int n_written_good = fwrite(f_avg_ideal_, sizeof(double), 3, forces_);
    // int n_written_flaw = fwrite(f_avg_flawed_, sizeof(double), 3, forces_);
    i_datapoint_++;
    for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
      double pos[3];
      int adj[6];
      Vertex *vrt{&vrts_[i_vrt]};
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        pos[i_dim] = vrt->GetPosition()[i_dim];
      }
      for (int i_neighb{0}; i_neighb < vrt->n_neighbs_; i_neighb++) {
        adj[i_neighb] = (int)vrt->neighbs_[i_neighb]->vid_;
      }
      // padding for 5-point vertices
      for (int i_neighb{vrt->n_neighbs_}; i_neighb < 6; i_neighb++) {
        adj[i_neighb] = -1;
      }
      fwrite(pos, sizeof(double), 3, vertices_);
      fwrite(adj, sizeof(int), 6, adjacency_);
    }
  } else {
    printf("data collection done! (%i datapoints written)\n", i_datapoint_);
    fclose(forces_);
    fclose(vertices_);
    fclose(adjacency_);
    early_exit = true;
  }
}

void TriMesh::UpdatePositions() {

  if (do_not_pass_go_) {
    return;
  }
  // shrink this jonny
  if (params_->mesh_shrink_rate > 0.0) {
    l_avg_ *= (1.0 - params_->mesh_shrink_rate);
    l_c0_ = 1.2 * l_avg_;
    l_c1_ = 0.8 * l_avg_;
    l_max_ = 1.4 * l_avg_;
    l_min_ = 0.6 * l_avg_;
  }
  UpdateTriangles();
  FlipEdges();
  if (do_not_pass_go_) {
    return;
  }
  UpdateMesh(); // update edge lengths, triangle area/vol, etc.
  ApplyMembraneForces();
  // ApplyBoundaryForces();
  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    double sigma{sqrt(2 * params_->delta / gamma_)}; // kbT = 1??
    double dr[3] = {0.0, 0.0, 0.0};
    double const *const r_prev = vrts_[i_vrt].GetPosition();
    double r_final[3];
    for (int i = 0; i < 3; ++i) {
      double vel{vrts_[i_vrt].GetForce()[i] / gamma_};
      double noise{rng_->RandomNormal(sigma)};
      // double noise{0.0};
      dr[i] = vel * params_->delta + noise;
      r_final[i] = r_prev[i] + dr[i];
    }
    vrts_[i_vrt].SetPos(r_final);
  }
  // WriteOutputs();
}