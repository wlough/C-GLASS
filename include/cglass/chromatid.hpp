#ifndef _CGLASS_CHROMATID_H_
#define _CGLASS_CRROMATID_H_

#include "kinetochore.hpp"
class Chromosome;

class Chromatid : public Object {

private:
  Kinetochore kc;

  chromosome_parameters *sparams_;

  friend Chromosome;

public:
  Chromatid(unsigned long seed) : Object(seed), kc(seed) {
    SetSID(species_id::chromosome);
    printf("BONK\n");
  }
  void Init(chromosome_parameters *sparams) {
    sparams_ = sparams;
    color_ = sparams->color;
    draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
    length_ = sparams->length;
    diameter_ = sparams->diameter;
  }
};

#endif