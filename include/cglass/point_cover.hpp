#ifndef _CGLASS_POINT_COVER_H_
#define _CGLASS_POINT_COVER_H_

#include "sphere.hpp"

class PointCover : public Composite {
  std::vector<Sphere*> sphere_ptrs_;
};

#endif  // _CGLASS_POINT_COVER_H_