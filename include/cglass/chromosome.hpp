#ifndef _CGLASS_CHROMOSOME_H_
#define _CGLASS_CHROMOSOME_H_

#include "chromatid.hpp"
#include "object.hpp"
#include "tracker.hpp"

// By convention, the Object that 'Chromosome' represents
// is the centromere since it is a logical mid-way point
class Chromosome : public Object {
private:
  bool zero_temperature_ = false;

  std::vector<Chromatid> sisters_;

  // Add centromere as its own distinct object?

  chromosome_parameters *sparams_;

  // VV MOVE TO CHROMATID ??
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

  void SetPosition(double const *const new_pos) {
    Object::SetPosition(new_pos);
    for (int i_sis{0}; i_sis < sisters_.size(); i_sis++) {
      double pos_sis[3]{new_pos[0], new_pos[1], new_pos[2]};
      for (int i{0}; i < 3; i++) {
        double sign{i_sis == 0 ? -1.0 : 1.0};
        pos_sis[i] += 0.8 * sign * diameter_;
      }
      sisters_[i_sis].SetPosition(pos_sis);
    }
    // for (auto &&sis : sisters_) {
    //   sis.SetPosition(new_pos);
    // }
  }
  void Draw(std::vector<graph_struct *> &graph_array) {
    Object::Draw(graph_array);
    for (auto &&sis : sisters_) {
      sis.Draw(graph_array);
    }
  }
  void GetInteractors(std::vector<Object *> &ixors) {
    ixors.push_back(this);
    for (auto &&sis : sisters_) {
      ixors.push_back(&sis);
    }
  }
  void UpdatePosition();
};

#endif