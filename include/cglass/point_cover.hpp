#ifndef _CGLASS_POINT_COVER_H_
#define _CGLASS_POINT_COVER_H_

#include "sphere.hpp"
#include "composite.hpp"

class PointCover : public Composite {
  protected:
    std::vector<Sphere*> sphere_ptrs_;
  public:
    PointCover(unsigned long seed);
    void AddSpherePtr(Sphere* s);
};

#endif  // _CGLASS_POINT_COVER_H_