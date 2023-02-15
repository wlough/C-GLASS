#ifndef _CGLASS_CHROMATID_H_
#define _CGLASS_CRROMATID_H_

#include "kinetochore.hpp"

class Chromatid : public Object {

private:
  Kinetochore kc;

public:
  Chromatid(unsigned long seed) : Object(seed), kc(seed) {
    SetSID(species_id::chromosome);
  }
  void Init() {}
};

#endif