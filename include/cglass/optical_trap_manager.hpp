/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : optical_trap_manager
 * @created     : Tuesday Jul 07, 2020 11:22:58 MDT
 */

#ifndef _CGLASS_OPTICAL_TRAP_MANAGER_HPP_

#define _CGLASS_OPTICAL_TRAP_MANAGER_HPP_

#include "output_manager.hpp"
#include "species_factory.hpp"

class OpticalTrapOutputManager : public OutputManagerBase<OpticalTrapSpecies> {
  /* Do not write thermo - base output manager will handle that */
  void WriteThermo() {}
  void ReadThermo() {}
  void InitThermo(std::string fname) {}
};

class OpticalTrapManager {
private:
  /* private data */
  system_parameters *params_;
  OpticalTrapOutputManager output_mgr_;
  std::vector<Object *> *objs_;
  space_struct *space_;
  bool update_;

public:
  void Init(system_parameters *params, space_struct *space,
            std::vector<Object *> *objs);
};

#endif

