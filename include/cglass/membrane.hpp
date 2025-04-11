#pragma once

#include "mesh.hpp"
#include "meshbrane/membrane.hpp"

class Membrane : public Mesh, public meshbrane::Membrane {
public:
  membrane_parameters *sparams_;
  Membrane(unsigned long seed);
  void SetParameters();
  void Init(membrane_parameters *sparams);
};

// typedef std::vector<Membrane>::iterator membrane_iterator;
// typedef std::vector<
//     std::pair<std::vector<Membrane>::iterator, std::vector<Membrane>::iterator>>
//     membrane_chunk_vector;