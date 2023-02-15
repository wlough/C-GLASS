#ifndef _CGLASS_CHROMOSOME_H_
#define _CGLASS_CHROMOSOME_H_

#include "chromatid.hpp"
#include "object.hpp"
#include "tracker.hpp"

class Chromosome : public Object {
private:
  bool zero_temperature_ = false;

  std::vector<Chromatid> sisters_;

  chromosome_parameters *sparams_;

  // VV MOVE TO CHROMATID
  double gamma_trans_ = 0.0;
  double gamma_perp_ = 0;
  double gamma_rot_ = 0.0;
  double noise_tr_ = 0.0;
  double noise_rot_ = 0.0;
  double diffusion_ = 0.0;
  double diffusion_rot_ = 0.0;
  double body_frame_[6];

public:
private:
  void GetBodyFrame();

public:
  Chromosome(unsigned long seed);
  void Init(chromosome_parameters *sparams);

  void UpdatePosition();
};

#endif