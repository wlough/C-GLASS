#ifndef _CGLASS_COMPOSITE_H_
#define _CGLASS_COMPOSITE_H_

#include "cglass/object.hpp"

class Composite : public Object {
 private:
  static int _next_comp_id_;
  static std::mutex _comp_mtx_;
  void InitCompID();
 public:
  Composite(unsigned long seed);
};

#endif  // _CGLASS_COMPOSITE_H_
