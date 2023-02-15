#ifndef _CGLASS_KINETOCHORE_H_
#define _CGLASS_KINETOCHORE_H_

#include "object.hpp"

class Kinetochore : public Object {

public:
  Kinetochore(unsigned long seed) : Object(seed) {
    SetSID(species_id::chromosome);
  }
  void Init() {}
};
#endif