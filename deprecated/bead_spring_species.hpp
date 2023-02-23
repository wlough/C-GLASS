#ifndef _CGLASS_BEAD_SPRING_SPECIES_H_
#define _CGLASS_BEAD_SPRING_SPECIES_H_

#include "species.hpp"
#include "bead_spring.hpp"

typedef std::vector<BeadSpring>::iterator bs_iterator;
typedef std::vector<std::pair<std::vector<BeadSpring>::iterator,
                              std::vector<BeadSpring>::iterator> >
    bead_spring_chunk_vector;

class BeadSpringSpecies : public Species<BeadSpring, species_id::bead_spring> {
 protected:
  // Analysis structures
  bool midstep_;
  double e_bend_, tot_angle_, mse2e_, mse2e2_, avg_clen_;
  int **theta_histogram_;
  int time_, n_bins_, n_samples_;
  std::fstream spiral_file_, theta_file_, mse2e_file_, diffusion_file_;

 public:
  BeadSpringSpecies() : Species() { SetSID(species_id::bead_spring); }
  void Init(system_parameters *params, species_base_parameters *sparams, SpaceBase *space) {
    Species::Init(params, sparams, space);
    midstep_ = true;
    if (sparams_.packing_fraction > 0) {
      if (sparams_.length <= 0) {
        Logger::Error(
            "Packing fraction with polydisperse lengths not implemented yet\n");
      }
      if (params_->n_dim == 2) {
        double fil_vol =
            sparams_.length * sparams_.diameter +
            0.25 * M_PI * SQR(sparams_.diameter);
        sparams_.num =
            sparams_.packing_fraction * space_->volume / fil_vol;
      } else {
        double fil_vol = 0.25 * M_PI * SQR(sparams_.diameter) *
                             sparams_.length +
                         M_PI * CUBE(sparams_.diameter) / 6.0;
        sparams_.num =
            sparams_.packing_fraction * space_->volume / fil_vol;
      }
    }
  }
  void InitAnalysis();
  void InitDiffusionAnalysis();
  void InitThetaAnalysis();
  void InitMse2eAnalysis();
  void RunAnalysis();
  void RunThetaAnalysis();
  void RunMse2eAnalysis();
  void FinalizeAnalysis();
  void FinalizeDiffusionAnalysis();
  void FinalizeMse2eAnalysis();
  void FinalizeThetaAnalysis();
  void UpdatePositions() {
#ifdef ENABLE_OPENMP
    int max_threads = omp_get_max_threads();
    bead_spring_chunk_vector chunks;
    chunks.reserve(max_threads);
    size_t chunk_size = members_.size() / max_threads;
    bs_iterator cur_iter = members_.begin();
    for (int i = 0; i < max_threads - 1; ++i) {
      bs_iterator last_iter = cur_iter;
      std::advance(cur_iter, chunk_size);
      chunks.push_back(std::make_pair(last_iter, cur_iter));
    }
    chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks)
    {
#pragma omp for
      for (int i = 0; i < max_threads; ++i)
        for (auto it = chunks[i].first; it != chunks[i].second; ++it)
          it->UpdatePosition(midstep_);
    }
#else
    for (bs_iterator it = members_.begin(); it != members_.end(); ++it)
      it->UpdatePosition(midstep_);
#endif

    midstep_ = !midstep_;
  }
};

#endif
