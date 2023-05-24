#ifndef _CGLASS_POINT_COVER_H_
#define _CGLASS_POINT_COVER_H_

#include "object.hpp"
class Receptor;

class PointCover : public Object {
protected:
  std::vector<Receptor *> sphere_ptrs_;

public:
  PointCover(unsigned long seed);
  void AddSpherePtr(Receptor *s);
};

#endif // _CGLASS_POINT_COVER_H_