#ifndef _CGLASS_CORTEX_H_
#define _CGLASS_CORTEX_H_

#include "point_cover.hpp"

// Class for assigning mesh of binding sites on the cell cortex

class Cortex : public PointCover {
private:
public:
  Cortex(unsigned long seed);
  void UpdateInteractors();
};

#endif
