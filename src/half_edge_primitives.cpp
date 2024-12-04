#include <cglass/filament.hpp>
// #include <cglass/ply_tools.hpp>
#include "cglass/half_edge_primitives.hpp"
#include "meshbrane/simple_generator.hpp"
#include <coroutine>
#include <filesystem>
#include <iostream>
#include <meshbrane/matrix_mesh.hpp>
#include <meshbrane/meshbrane_data_types.hpp>
#include <unistd.h>
#include <vector>

#include <unordered_set> // std::unordered_set

// namespace half_edge {

///////////////////////////////////////////////////
// Vertex /////////////////////////////////////////
///////////////////////////////////////////////////
HalfEdgeGenerator Vertex::generate_H_out_clockwise() {
  return h_->generate_H_rotcw();
}
HalfEdgeGenerator Vertex::generate_H_out_clockwise(HalfEdge *h_start) {
  if (*h_start->v_origin() != *this) {
    printf("Error in Vertex::generate_H_out_clockwise\n");
    printf("h_start does not originate at this vertex\n");
    exit(1);
  }
  return h_start->generate_H_rotcw();
};
TriangleGenerator Vertex::generate_F_incident_clockwise() {
  HalfEdgeGenerator h = generate_H_out_clockwise();
  for (HalfEdge *h : h) {
    if (h->is_in_some_negative_boundary()) {
      continue;
    }
    co_yield h->f_left();
  }
}

void Vertex::UpdateNeighbors() {
  neighbs_.clear();
  edges_.clear();
  tris_.clear();
  for (HalfEdge *h : generate_H_out_clockwise()) {
    neighbs_.push_back(h->v_head());
    edges_.push_back(h->e_parallel());
    if (h->is_in_some_negative_boundary()) {
      continue;
    }
    tris_.push_back(h->f_left());
  }
  n_neighbs_ = neighbs_.size();
  n_edges_ = edges_.size();
  n_tris_ = tris_.size();
}

///////////////////////////////////////////////////
// HalfEdge ///////////////////////////////////////
///////////////////////////////////////////////////
HalfEdgeGenerator HalfEdge::generate_H_rotcw() {
  HalfEdge *h = this;
  HalfEdge *h_start = h;
  do {
    co_yield h;
    h = h->h_rotcw();
  } while (h != h_start);
};
///////////////////////////////////////////////////
// Edge ///////////////////////////////////////////
///////////////////////////////////////////////////
void Edge::UpdateNeighborTris() {
  tris_[0] = h_parallel()->f_left();
  tris_[1] = h_parallel()->h_twin()->f_left();
}
bool Edge::is_in_some_boundary() const {
  return h_parallel()->is_in_some_negative_boundary() ||
         h_parallel()->is_in_some_positive_boundary();
}
/**
 * @brief 
 * 
 * @return true 
 * @return false 
 */
bool Edge::is_flippable() const {
  // if self.boundary_contains_h(h):
  //     return False
  // hlj = h
  // hjk = self.h_next_h(hlj)
  // # hjl = self.h_twin_h(hlj)
  // hli = self.h_next_h(self.h_twin_h(hlj))
  // vi = self.v_head_h(hli)
  // vk = self.v_head_h(hjk)

  // for him in self.generate_H_out_v_clockwise(vi):
  //     if self.v_head_h(him) == vk:
  //         return False
  if (is_in_some_boundary()) {
    return false;
  }
  HalfEdge *hlj = h_;
  HalfEdge *hjk = hlj->h_next();
  HalfEdge *hli = hlj->h_twin()->h_next();
  Vertex *vi = hli->v_head();
  Vertex *vk = hjk->v_head();

  // for (HalfEdgePtr him : hlj->generate_H_out_v_clockwise(vi)) {
  //   if (him->v_head() == vk) {
  //     return false;
  //   }
  // }

  return true;
}

///////////////////////////////////////////////////
// Triangles //////////////////////////////////////
///////////////////////////////////////////////////
void Triangle::UpdateNeighborEdges() {
  edges_[0] = h_right()->e_parallel();
  edges_[1] = h_right()->h_next()->e_parallel();
  edges_[2] = h_right()->h_next()->h_next()->e_parallel();
}

///////////////////////////////////////////////////
// TriMesh ////////////////////////////////////////
///////////////////////////////////////////////////

// } // namespace half_edge