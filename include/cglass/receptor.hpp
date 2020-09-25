#ifndef _CGLASS_RECEPTOR_H_
#define _CGLASS_RECEPTOR_H_

#include "site.hpp"
#include "cortex.hpp"

// Receptors are modelled as sites on a mesh that can
// bind to one side of a crosslinking protein. 

class Receptor: public Site {
private:
  receptor_parameters *sparams_;
public:
  Receptor(unsigned long seed);
  void Init(receptor_parameters *sparams);
  void SetParameters();
};

#endif
