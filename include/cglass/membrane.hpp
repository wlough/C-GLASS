#pragma once

#include "half_edge_primitives.hpp"
#include "mesh.hpp"
#include "meshbrane/membrane.hpp"

class Membrane : public Mesh, public meshbrane::Membrane {
public:
  std::vector<Vertex> vrts_;
  membrane_parameters *sparams_;
  std::string ply_path_{"none"};
  Membrane(unsigned long seed);
  void SetParameters();
  void Init(membrane_parameters *sparams);
  void Draw(std::vector<graph_struct *> &graph_array);
  void LoadPly();
  void SyncWithMats();
  void UpdatePosition();
};

typedef std::vector<Membrane>::iterator membrane_iterator;
typedef std::vector<
    std::pair<std::vector<Membrane>::iterator, std::vector<Membrane>::iterator>>
    membrane_chunk_vector;