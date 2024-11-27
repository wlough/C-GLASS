#pragma once

#include <cglass/object.hpp>
#include <meshbrane/matrix_mesh.hpp>
#include <meshbrane/meshbrane_data_types.hpp>

namespace mbrn = meshbrane;

/**
 * @brief A dynamically triangulated membrane.
 * 
 */
class Membrane : public Object, public mbrn::MatrixMesh {
public:
  Membrane() = default;
  ~Membrane() = default;
};
