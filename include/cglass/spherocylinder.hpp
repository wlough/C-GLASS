#ifndef _CGLASS_SPHEROCYLINDER_H_
#define _CGLASS_SPHEROCYLINDER_H_

#include "br_rod.hpp"
class Spherocylinder : public BrRod {
 protected:
  spherocylinder_parameters *sparams_;
  void ApplyForcesTorques();
  void InsertSpherocylinder();
  void SetParameters();

 public:
  Spherocylinder(unsigned long seed);
  void Init(spherocylinder_parameters *sparams);
  void UpdatePosition();
};

#endif  // _CGLASS_SPHEROCYLINDER_H_
