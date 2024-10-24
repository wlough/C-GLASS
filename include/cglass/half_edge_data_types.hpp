#ifndef HALF_EDGE_DATA_TYPES_HPP
#define HALF_EDGE_DATA_TYPES_HPP

#include <Eigen/Dense> //
#include <tuple>       // std::tuple

using Coords3d = Eigen::Vector3d;
using Samplesi = Eigen::Matrix<int, Eigen::Dynamic, 1>;
using Samples2i = Eigen::Matrix<int, Eigen::Dynamic, 2>;
using Samples3i = Eigen::Matrix<int, Eigen::Dynamic, 3>;
using Samples3d = Eigen::Matrix<double, Eigen::Dynamic, 3>;
using HalfEdgeSamples = std::tuple<Samples3d, Samplesi, Samplesi, Samplesi,
                                   Samplesi, Samplesi, Samplesi, Samplesi>;
using VertexFaceSamples = std::tuple<Samples3d, Samples3i>;
using VertexEdgeFaceSamples = std::tuple<Samples3d, Samples2i, Samples3i>;
#endif /* HALF_EDGE_DATA_TYPES_HPP */