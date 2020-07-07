/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : optical_trap_manager
 * @created     : Tuesday Jul 07, 2020 11:19:30 MDT
 */

#include "cglass/optical_trap_manager.hpp"

void OpticalTrapManager::Init(system_parameters *params, space_struct *space,
                              std::vector<Object *> *objs) {
  objs_ = objs;
  update_ = false;
  params_ = params;
  space_ = space;
}

