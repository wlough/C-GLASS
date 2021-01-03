#ifndef _CGLASS_CORTEX_H_
#define _CGLASS_CORTEX_H_

#include "mesh.hpp"

// Class for assigning mesh of binding sites on the cell cortex

class Cortex : public Mesh {
private:
public:
  Cortex(unsigned long seed);
  void InitSiteAt(double *new_pos, double d);
  void UpdateInteractors();
  void Draw(std::vector<graph_struct *> &graph_array);
};

#endif
