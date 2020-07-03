/**
 * @author      : Adam Lamson (adam.lamson@colorado.edu)
 * @file        : optical_trap
 * @created     : Friday Jul 03, 2020 14:49:10 MDT
 */

#include "cglass/optical_trap.hpp"

OpticalTrap::OpticalTrap(unsigned long seed) : Object(seed) {
  SetSID(species_id::optical_trap);
}

/*! \brief Initialize optical trap object
 *
 *  Detailed description
 *
 * \param *sparams Parameter description
 * \return Return parameter description
 */
void OpticalTrap::Init(optical_trap_parameters *sparams) {
  sparams_ = sparams;
  trap_spring_ = sparams_->trap_spring;
  diameter_ = sparams_->trap_diameter;
  bead_diameter_ = sparams_->trap_diameter;
  color_ = sparams_->trap_color;
  bead_color_ = sparams_->trap_color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
}

