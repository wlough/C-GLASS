#ifndef _CGLASS_CORTEX_H_
#define _CGLASS_CORTEX_H_

#include "mesh.hpp"

// Class for assigning mesh of binding sites on the cell cortex

class Cortex : public Mesh {
private:
  system_parameters *params_;
  double site_diameter_;
  double site_concentration_;
  int total_sites_;
public:
  Cortex(unsigned long seed);
  void Init(system_parameters *params);
  void SetParameters();
  void AddSites();
  
  void UpdateInteractors();
  void Draw(std::vector<graph_struct *> &graph_array);
};

#endif
