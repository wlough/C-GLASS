#ifndef _CGLASS_TRIANGLE_MESH_
#define _CGLASS_TRIANGLE_MESH_

// #include "common_libs.hpp"
// #include "definitions.hpp"
#include "rng.hpp"
#include "site.hpp"

// TODO add param storage and sync site seeds

struct Triangle;
struct Vertex : public Site {
  int seed{0};
  double pos_[3];

  int n_tris_ = 0;
  Triangle *tris_[10]; // triangles this vertex is a part of
  int n_neighbs_ = 0;
  std::vector<Vertex *> neighbs_;

  // bend energy jawn -- summed over all dem neighbs
  double sum_lsqT_;         // scalar
  double sum_del_lsqT_[3];  // vec; scalar in each dim
  double sum_rT_[3];        // vec; scalr in each dim
  double sum_del_rT_[3][3]; // tensor; vec. in each dim

  Vertex() : Site(seed) {
    pos_[0] = pos_[1] = pos_[2] = 0;
    Site::SetPositionXYZ(0, 0, 0);
  }
  Vertex(double x, double y, double z) : Site(seed) {
    pos_[0] = x;
    pos_[1] = y;
    pos_[2] = z;
    Site::SetPositionXYZ(x, y, z);
  }
  friend bool operator==(const Vertex &lhs, const Vertex &rhs) {
    return (lhs.pos_[0] == rhs.pos_[0] and lhs.pos_[1] == rhs.pos_[1] and
            lhs.pos_[2] == rhs.pos_[2]);
  }
  friend bool operator!=(const Vertex &lhs, const Vertex &rhs) {
    return !(lhs == rhs);
  }
  void SetPos(const double *const new_pos) {
    pos_[0] = position_[0] = new_pos[0];
    pos_[1] = position_[1] = new_pos[1];
    pos_[2] = position_[2] = new_pos[2];
  }
  void ZeroSums() {
    sum_lsqT_ = 0.0;
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      sum_rT_[i_dim] = 0.0;
      sum_del_lsqT_[i_dim] = 0.0;
      for (int j_dim{0}; j_dim < 3; j_dim++) {
        sum_del_rT_[i_dim][j_dim] = 0.0;
      }
    }
  }
};

struct Triangle {
  double area_;
  double nhat_[3];
  double color_[3];

  // euler angles
  double cosGamma_;
  double sinGamma_;
  double cosBeta_;
  double sinBeta_;
  double Zrot_;
  double XYrot_[2][3]; // [dim][i_vrt]

  Vertex *vrts_[3]; // verteces that compose this triangle
  Triangle *neighbs_[3];

  Triangle(Vertex *v1, Vertex *v2, Vertex *v3) {
    vrts_[0] = v1;
    vrts_[1] = v2;
    vrts_[2] = v3;
    color_[0] = rand() % 255;
    color_[1] = rand() % 255;
    color_[2] = rand() % 255;
  }
  // Returns the XOR between two booleans a, b
  inline bool b3dxor(bool a, bool b) { return ((a && !b) || (b && !a)); }
  // Inlined version of distance calc for ease
  // Function returns the parameter t for the point along the line
  // closest to x, y, z
  /**
 * @brief Inlined version of point distance to segment
 * @param (x) Point to find distance from
 * @param (a, b, c) Equation of line
 * @return parameter t for linear distance along line
 */
  inline double tri_point_to_line(double xStart, double yStart, double zStart,
                                  double a, double b, double c, double x,
                                  double y, double z) {
    return (a * (x - xStart) + b * (y - yStart) + c * (z - zStart)) /
           (a * a + b * b + c * c);
  }
  // inside triangle determines if xtest, ytest lies within the boundaries defined
  // by bx and by
  bool inside_triangle(double *bx, double *by, double xTest, double yTest) {
    bool inside_triangle = false;

    std::cout << "Testing (" << xTest << ", " << yTest << ") inside triangle\n";
    for (int iv = 0; iv < 3; ++iv) {
      std::cout << "Vertex[" << iv << "] (" << bx[iv] << ", " << by[iv]
                << ")\n";
    }
    // Compute signed area of each triangle between the point and an edge
    double area0 =
        (bx[0] - xTest) * (by[1] - yTest) - (bx[1] - xTest) * (by[0] - yTest);
    double area1 =
        (bx[1] - xTest) * (by[2] - yTest) - (bx[2] - xTest) * (by[1] - yTest);
    double area2 =
        (bx[2] - xTest) * (by[0] - yTest) - (bx[0] - xTest) * (by[2] - yTest);

    if (area0 != 0.0 && area1 != 0.0 && area2 != 0.0) {
      //std::cout << "Testing first case of all nonzero areas\n";
      inside_triangle = ((area0 > 0.0 && area1 > 0.0 && area2 > 0.0) ||
                         (area0 < 0.0 && area1 < 0.0 && area2 < 0.0));
      return inside_triangle;
    }
    if (area0 == 0.0 && area1 == 0.0 && area2 == 0.0) {
      std::cout << "ALL AREAS 0 NOT IMPLEMENTED YET!\n";
      exit(1);
    }
    //std::cout << "Testing third case of wacky triangles\n";
    inside_triangle = area0 == 0. and area1 > 0. and area2 > 0. or
                      area0 == 0. and area1 < 0. and area2 < 0. or
                      area1 == 0. and area0 > 0. and area2 > 0. or
                      area1 == 0. and area0 < 0. and area2 < 0. or
                      area2 == 0. and area0 > 0. and area1 > 0. or
                      area2 == 0. and area0 < 0. and area1 < 0. or
                      area0 == 0. and area1 == 0. or
                      area0 == 0. and area2 == 0. or
                      area0 == 1. and area2 == 0.;

    return inside_triangle;
  }
  // calcluate the distance between a pair of segments
  // INPUT:
  // (x1, y1, z1) describes the start/end of first segment
  // (x2, y2, z2) describes the start/end of second segment
  // OUTPUT:
  // (t1, t2) linear distance along segments of closest approach
  // dist the squared distance of closest approach
  void segment_dist(double xStart1, double yStart1, double zStart1,
                    double xEnd1, double yEnd1, double zEnd1, double xStart2,
                    double yStart2, double zStart2, double xEnd2, double yEnd2,
                    double zEnd2, double *t1, double *t2, double *dist) {
    //std::cout << "Segment distance\n";
    //std::cout << "  line1 ("
    //    << xStart1 << ", " << yStart1 << ", " << zStart1 << ")("
    //    << xEnd1 << ", " << yEnd1 << ", " << zEnd1 << ")\n";
    //std::cout << "  line2 ("
    //    << xStart2 << ", " << yStart2 << ", " << zStart2 << ")("
    //    << xEnd2 << ", " << yEnd2 << ", " << zEnd2 << ")\n";
    // First, find z1, z2 position of global minimum
    double a1 = xEnd1 - xStart1;
    double b1 = yEnd1 - yStart1;
    double c1 = zEnd1 - zStart1;
    double a2 = xEnd2 - xStart2;
    double b2 = yEnd2 - yStart2;
    double c2 = zEnd2 - zStart2;
    double sqr1 = a1 * a1 + b1 * b1 + c1 * c1;
    double sqr2 = a2 * a2 + b2 * b2 + c2 * c2;
    double crossTerm = -(a1 * a2 + b1 * b2 + c1 * c2);
    double const1 = a1 * (xStart2 - xStart1) + b1 * (yStart2 - yStart1) +
                    c1 * (zStart2 - zStart1);
    double const2 = -(a2 * (xStart2 - xStart1) + b2 * (yStart2 - yStart1) +
                      c2 * (zStart2 - zStart1));
    double den = sqr1 * sqr2 - crossTerm * crossTerm;
    double t1Numer = const1 * sqr2 - const2 * crossTerm;
    double t2Numer = sqr1 * const2 - const1 * crossTerm;

    if (fabs(den) < 1e-20 ||
        fabs(den) < 1e-6 * std::max(fabs(t1Numer), fabs(t2Numer))) {
      //std::cout << "  Found parallel lines\n";
      // Parallel lines: "just" check the 4 endpoints
      // start of line1 versus line2
      double t2Trunc, t1Trunc, tdist;
      t2Trunc = std::max(
          0.0,
          std::min(1.0, tri_point_to_line(xStart2, yStart2, zStart2, a2, b2, c2,
                                          xStart1, yStart1, zStart1)));
      *dist = (a2 * t2Trunc + xStart2 - xStart1) *
                  (a2 * t2Trunc + xStart2 - xStart1) +
              (b2 * t2Trunc + yStart2 - yStart1) *
                  (b2 * t2Trunc + yStart2 - yStart1) +
              (c2 * t2Trunc + zStart2 - zStart1) *
                  (c2 * t2Trunc + zStart2 - zStart1);
      *t1 = 0.0;
      *t2 = t2Trunc;
      // End of line1 versus line2
      t2Trunc = std::max(
          0.0, std::min(1.0, tri_point_to_line(xStart2, yStart2, zStart2, a2,
                                               b2, c2, xEnd1, yEnd1, zEnd1)));
      tdist =
          (a2 * t2Trunc + xStart2 - xEnd1) * (a2 * t2Trunc + xStart2 - xEnd1) +
          (b2 * t2Trunc + yStart2 - yEnd1) * (b2 * t2Trunc + yStart2 - yEnd1) +
          (c2 * t2Trunc + zStart2 - zEnd1) * (c2 * t2Trunc + zStart2 - zEnd1);
      if (tdist < *dist) {
        *dist = tdist;
        *t1 = 1.0;
        *t2 = t2Trunc;
      }
      // Start of line2 versus line1
      t1Trunc = std::max(
          0.0,
          std::min(1.0, tri_point_to_line(xStart1, yStart1, zStart1, a1, b1, c1,
                                          xStart2, yStart2, zStart2)));
      tdist = (a1 * t1Trunc + xStart1 - xStart2) *
                  (a1 * t1Trunc + xStart1 - xStart2) +
              (b1 * t1Trunc + yStart1 - yStart2) *
                  (b1 * t1Trunc + yStart1 - yStart2) +
              (c1 * t1Trunc + zStart1 - zStart2) *
                  (c1 * t1Trunc + zStart1 - zStart2);
      if (tdist < *dist) {
        *dist = tdist;
        *t1 = t1Trunc;
        *t2 = 0.0;
      }
      // End of line2 versus line1
      t1Trunc = std::max(
          0.0, std::min(1.0, tri_point_to_line(xStart1, yStart1, zStart1, a1,
                                               b1, c1, xEnd2, yEnd2, zEnd2)));
      tdist =
          (a1 * t1Trunc + xStart1 - xEnd2) * (a1 * t1Trunc + xStart1 - xEnd2) +
          (b1 * t1Trunc + yStart1 - yEnd2) * (b1 * t1Trunc + yStart1 - yEnd2) +
          (c1 * t1Trunc + zStart1 - zEnd2) * (c1 * t1Trunc + zStart1 - zEnd2);
      if (tdist < *dist) {
        *dist = tdist;
        *t1 = t1Trunc;
        *t2 = 1.0;
      }
    } else {
      // Non parallel lines
      //std::cout << "  Non parallel lines\n";
      *t1 = t1Numer / den;
      *t2 = t2Numer / den;
      bool out1 = (*t1 < 0.0 || *t1 > 1.0);
      bool out2 = (*t2 < 0.0 || *t2 > 1.0);
      if (out1 && out2) {
        // If both closest points are out of bounds, truncate each one to
        // its segment, then find closest point on other segment to that
        // truncated point. If this gives different answers, pick the pair
        // with the closest approach
        //std::cout << "  Both closest points out of bounds\n";
        *t1 = std::max(0.0, std::min(1.0, *t1));
        double t2Trunc = std::max(
            0.0, std::min(1.0, tri_point_to_line(xStart2, yStart2, zStart2, a2,
                                                 b2, c2, a1 * *t1 + xStart1,
                                                 b1 * *t1 + yStart1,
                                                 c1 * *t1 + zStart1)));
        *t2 = std::max(0.0, std::min(1.0, *t2));
        double t1Trunc = std::max(
            0.0, std::min(1.0, tri_point_to_line(xStart1, yStart1, zStart1, a1,
                                                 b1, c1, a2 * *t2 + xStart2,
                                                 b2 * *t2 + yStart2,
                                                 c2 * *t2 + zStart2)));
        if (*t1 != t1Trunc || *t2 != t2Trunc) {
          //std::cout << "  or not equals\n";
          *dist = (a1 * *t1 + xStart1 - a2 * t2Trunc - xStart2) *
                      (a1 * *t1 + xStart1 - a2 * t2Trunc - xStart2) +
                  (b1 * *t1 + yStart1 - b2 * t2Trunc - yStart2) *
                      (b1 * *t1 + yStart1 - b2 * t2Trunc - yStart2) +
                  (c1 * *t1 + zStart1 - c2 * t2Trunc - zStart2) *
                      (c1 * *t1 + zStart1 - c2 * t2Trunc - zStart2);
          double tdist = (a1 * t1Trunc + xStart1 - a2 * *t2 - xStart2) *
                             (a1 * t1Trunc + xStart1 - a2 * *t2 - xStart2) +
                         (b1 * t1Trunc + yStart1 - b2 * *t2 - yStart2) *
                             (b1 * t1Trunc + yStart1 - b2 * *t2 - yStart2) +
                         (c1 * t1Trunc + zStart1 - c2 * *t2 - zStart2) *
                             (c1 * t1Trunc + zStart1 - c2 * *t2 - zStart2);
          if (*dist < tdist) {
            *t2 = t2Trunc;
          } else {
            *dist = tdist;
            *t1 = t1Trunc;
          }
        } else {
          //std::cout << "   simple dist\n";
          *dist = (a1 * *t1 + xStart1 - a2 * *t2 - xStart2) *
                      (a1 * *t1 + xStart1 - a2 * *t2 - xStart2) +
                  (b1 * *t1 + yStart1 - b2 * *t2 - yStart2) *
                      (b1 * *t1 + yStart1 - b2 * *t2 - yStart2) +
                  (c1 * *t1 + zStart1 - c2 * *t2 - zStart2) *
                      (c1 * *t1 + zStart1 - c2 * *t2 - zStart2);
        }
      } else if (out1) {
        //std::cout << "  Segment point 1 outside bounds\n";
        *t1 = std::max(0.0, std::min(1.0, *t1));
        *t2 = std::max(
            0.0, std::min(1.0, tri_point_to_line(xStart2, yStart2, zStart2, a2,
                                                 b2, c2, a1 * *t1 + xStart1,
                                                 b1 * *t1 + yStart1,
                                                 c1 * *t1 + zStart1)));
        *dist = (a1 * *t1 + xStart1 - a2 * *t2 - xStart2) *
                    (a1 * *t1 + xStart1 - a2 * *t2 - xStart2) +
                (b1 * *t1 + yStart1 - b2 * *t2 - yStart2) *
                    (b1 * *t1 + yStart1 - b2 * *t2 - yStart2) +
                (c1 * *t1 + zStart1 - c2 * *t2 - zStart2) *
                    (c1 * *t1 + zStart1 - c2 * *t2 - zStart2);
      } else if (out2) {
        //std::cout << "  Segment point 2 outside bounds\n";
        *t2 = std::max(0.0, std::min(1.0, *t2));
        *t1 = std::max(
            0.0, std::min(1.0, tri_point_to_line(xStart1, yStart1, zStart1, a1,
                                                 b1, c1, a2 * *t2 + xStart2,
                                                 b2 * *t2 + yStart2,
                                                 c2 * *t2 + zStart2)));
        *dist = (a1 * *t1 + xStart1 - a2 * *t2 - xStart2) *
                    (a1 * *t1 + xStart1 - a2 * *t2 - xStart2) +
                (b1 * *t1 + yStart1 - b2 * *t2 - yStart2) *
                    (b1 * *t1 + yStart1 - b2 * *t2 - yStart2) +
                (c1 * *t1 + zStart1 - c2 * *t2 - zStart2) *
                    (c1 * *t1 + zStart1 - c2 * *t2 - zStart2);
      } else {
        //std::cout << "  Segment points both inside bounds\n";
        *dist = (a1 * *t1 + xStart1 - a2 * *t2 - xStart2) *
                    (a1 * *t1 + xStart1 - a2 * *t2 - xStart2) +
                (b1 * *t1 + yStart1 - b2 * *t2 - yStart2) *
                    (b1 * *t1 + yStart1 - b2 * *t2 - yStart2) +
                (c1 * *t1 + zStart1 - c2 * *t2 - zStart2) *
                    (c1 * *t1 + zStart1 - c2 * *t2 - zStart2);
      }
    }
  }
  // line segment to triangle calculation
  // INPUT
  // (x, y, z) segment line begin/end
  // tri, itriang triangle to calculate distance to
  // OUTPUT
  // (xRot, yRot) rotated point of closest approach to triangle
  // dist distance of closest approach
  // t fractional distance along line
  void MinDist_Segment(double xStart, double yStart, double zStart, double xEnd,
                       double yEnd, double zEnd, double *t, double *xRot,
                       double *yRot,
                       double *dist) { // Variables needed for running
    bool startIn, endIn;
    double temp;
    double xStartRot, yStartRot, zStartRot;
    double xEndRot, yEndRot, zEndRot;

    // Rotate the endpoints
    temp = xStart * cosGamma_ - yStart * sinGamma_;
    xStartRot = temp * cosBeta_ - zStart * sinBeta_;
    yStartRot = xStart * sinGamma_ + yStart * cosGamma_;
    zStartRot = temp * sinBeta_ + zStart * cosBeta_;
    temp = xEnd * cosGamma_ - yEnd * sinGamma_;
    xEndRot = temp * cosBeta_ - zEnd * sinBeta_;
    yEndRot = xEnd * sinGamma_ + yEnd * cosGamma_;
    zEndRot = temp * sinBeta_ + zEnd * cosBeta_;
    std::cout << "xStartRot (" << xStartRot << ", " << yStartRot << ", "
              << zStartRot << ")\n";
    std::cout << "xEndRot   (" << xEndRot << ", " << yEndRot << ", " << zEndRot
              << ")\n";
    startIn =
        inside_triangle(&(XYrot_[0][0]), &(XYrot_[1][0]), xStartRot, yStartRot);
    endIn = inside_triangle(&(XYrot_[0][0]), &(XYrot_[1][0]), xEndRot, yEndRot);
    //std::cout << "start/end in (" << startIn << ", " << endIn << ")\n";

    //// TEST ROTATION UNDO
    //double xprime = xStartRot * cosBeta[itriang] + zStartRot * sinBeta[itriang];
    //double x2 = xprime * cosGamma[itriang] + yStartRot * sinGamma[itriang];
    //double y2 = -xprime * sinGamma[itriang] + yStartRot * cosGamma[itriang];
    //double z2 = -xStartRot * sinBeta[itriang] + zStartRot * cosBeta[itriang];
    //std::cout << "start undo (" << x2 << ", " << y2 << ", " << z2 << ")\n";
    //// END TEST

    if (startIn && endIn) {
      // if both endpoints are over triangle, then one must be closest
      // unless line passes through triangle
      if (b3dxor(zStartRot > Zrot_, zEndRot > Zrot_)) {
        //std::cout << "b3dxor triggered\n";
        *dist = 0;
        *t = (Zrot_ - zStartRot) / (zEndRot - zStartRot);
        *xRot = (1.0 - *t) * xStartRot + *t * xEndRot;
        *yRot = (1.0 - *t) * yStartRot + *t * yEndRot;
      } else {
        //std::cout << "b3dxor not triggered\n";
        double distEnd = fabs(zStartRot - Zrot_);
        double distStart = fabs(zEndRot - Zrot_);
        if (distEnd < distStart) {
          *dist = distEnd;
          *xRot = xStartRot;
          *yRot = yStartRot;
          *t = 0.;
        } else {
          *dist = distStart;
          *xRot = xEndRot;
          *yRot = yEndRot;
          *t = 1.;
        }
      }
      return;
    }
    // Set dist to something absurd
    *dist = 1e10;

    // If one endpoint is over, it is a candidate
    if (startIn) {
      *xRot = xStartRot;
      *yRot = yStartRot;
      *dist = fabs(zStartRot - Zrot_);
      *t = 0.;
    }
    if (endIn) {
      *xRot = xEndRot;
      *yRot = yEndRot;
      *dist = fabs(zEndRot - Zrot_);
      *t = 1.;
    }

    // but still need to check each line segment
    double dmin = *dist * *dist;
    int ivMin = -1;
    int ivNextMin = -1;
    double t1AtMin, t2AtMin;
    for (int iv = 0; iv < 3; ++iv) {
      int ivNext = iv + 1;
      if (ivNext > 2)
        ivNext = 0;
      double t1, t2, dsqr;
      segment_dist(XYrot_[0][iv], XYrot_[1][iv], Zrot_, XYrot_[0][ivNext],
                   XYrot_[1][ivNext], Zrot_, xStartRot, yStartRot, zStartRot,
                   xEndRot, yEndRot, zEndRot, &t1, &t2, &dsqr);
      if (dsqr < dmin) {
        //std::cout << "  found new minimum distance " << dsqr << std::endl;
        ivMin = iv;
        ivNextMin = ivNext;
        dmin = dsqr;
        t1AtMin = t1;
        t2AtMin = t2;
        //std::cout << "  set t1 " << t1AtMin << ", t2 " << t2AtMin << std::endl;
      }
    }

    // If a segment was it, need squre root and rotated position
    if (ivMin >= 0) {
      *dist = sqrt(dmin);
      *xRot =
          (1.0 - t1AtMin) * XYrot_[0][ivMin] + t1AtMin * XYrot_[0][ivNextMin];
      *yRot =
          (1.0 - t1AtMin) * XYrot_[1][ivMin] + t1AtMin * XYrot_[1][ivNextMin];
      *t = t2AtMin;
    }
  }
};

class TriMesh {

private:
  size_t n_datapoints{0};
  bool file_open_{false};
  FILE *forces_{nullptr};
  system_parameters *params_{nullptr};
  double l_avg_{0.0};
  double gamma_{0.0};

  // params for radial force
  double kappa_B_{0.0};
  double l_max_{0.0};
  double l_min_{0.0};
  double l_c0_{0.0};
  double l_c1_{0.0};
  // params for bending force
  double kappa_{0.0};
  // params for area force
  double kappa_l_{0.0};
  double A_prime_{0.0};
  RNG *rng_; // SF TODO link with system RNG

public:
  std::vector<Object *> neighbs_;
  double r_sys_{0.0};
  std::vector<Vertex> vrts_;
  std::vector<Triangle> tris_;

private:
  void MakeIcosphere();
  void MakeIcosahedron();
  void DivideFaces();
  void ProjectToUnitSphere();
  void UpdateNeighbors();
  int SegmentToPolygon(double xStart, double yStart, double zStart, double xEnd,
                       double yEnd, double zEnd, double *t, double *xRot,
                       double *yRot, double *dist, double *rcontact);
  void UpdateMesh();
  void ApplyBoundaryForces();

public:
  TriMesh() {}
  void Init(system_parameters *params);
  int MinDist_Sphero(double *r_1, double *s_1, double *u_1, double length_1,
                     double *rmin, double *rminmag2, double *rcontact,
                     double *mu);
  void Draw(std::vector<graph_struct *> &graph_array);
  void UpdatePositions();
  void WriteOutputs();
};

#endif