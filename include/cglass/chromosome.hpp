#ifndef _CGLASS_CHROMOSOME_H_
#define _CGLASS_CHROMOSOME_H_

#include "object.hpp"
#include "tracker.hpp"

class Chromosome : public Object {
private:
  chromosome_parameters *sparams_;
  bool zero_temperature_ = false;
  double gamma_trans_ = 0.0;
  double gamma_rot_ = 0.0;
  double noise_tr_ = 0.0;
  double noise_rot_ = 0.0;
  double diffusion_ = 0.0;
  double diffusion_rot_ = 0.0;

public:
private:
public:
  Chromosome(unsigned long seed);
  void Init(chromosome_parameters *sparams);

  void UpdatePosition();
};

#endif