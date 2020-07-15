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

/*! \brief Set position of optical trap relative to obj minus end
 *
 *  Detailed description
 *
 * \param Object *obj Parameter description
 * \return Return parameter description
 */
void OpticalTrap::InsertAndAttach(Object *obj) {
  const double *obj_pos = obj->GetPosition();
  const double *obj_orient = obj->GetOrientation();
  const double obj_len = obj->GetLength() * .5;

  for (int i = 0; i < 3; ++i) {
    position_[i] = obj_pos[i] - obj_orient[i] * obj_len;
  }
  UpdatePeriodic();

  return;
}

