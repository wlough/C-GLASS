/**
 * @author      : Adam Lamson (adam.lamson@colorado.edu)
 * @file        : optical_trap
 * @created     : Friday Jul 03, 2020 14:49:10 MDT
 */

#ifndef OPTICAL_TRAP_HPP

#define OPTICAL_TRAP_HPP

#include "anchor.hpp"

class OpticalTrap : public Object {
private:
  optical_trap_parameters *sparams_;
  double trap_spring_;

  double bead_diameter_;
  double bead_color_;

public:
  OpticalTrap(unsigned long seed);
  void Init(optical_trap_parameters *sparams);

  //virtual ~optical_trap();
};

#endif /* end of include guard OPTICAL_TRAP_HPP */

