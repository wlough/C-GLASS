#include "cglass/membrane_species.hpp"

MembraneSpecies::MembraneSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::membrane);
}

void MembraneSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
}

// void MembraneSpecies::InitMembers() {
//   // what the fuck is this
//   return;
//   if ((sparams_.radius_of_curvature > 0 || sparams_.intrinsic_curvature != 0) &&
//       sparams_.randomize_intrinsic_curvature_handedness) {
//     int n_neg = n_members_ / 2;
//     if (n_members_ % 2 == 1 && rng_.RandomUniform() > 0.5) {
//       n_neg += 1;
//     }
//     int *neg = new int[n_members_];
//     for (int i = 0; i < n_members_; ++i) {
//       if (i >= n_neg)
//         neg[i] = 0;
//       else
//         neg[i] = 1;
//     }
//     rng_.Shuffle<int>(neg, n_members_);
//     int k = 0;

//     for (int i = 0; i < n_members_; ++i) {
//       if (neg[i] > 0) {
//         k++;
//         const double c = members_[i].GetCurvature();
//         members_[i].SetCurvature(-c);
//         members_[i].SetColor(-c > 0 ? sparams_.color : sparams_.color + M_PI,
//                              draw_type::fixed);
//       }
//     }
//     delete[] neg;
//   }
// }

// void MembraneSpecies::CleanUp() {
//   Species::CleanUp();
//   if (error_file_.is_open()) {
//     error_file_.close();
//   }
// }

// void MembraneSpecies::UpdatePositions() {
// #ifdef ENABLE_OPENMP
//   int max_threads = omp_get_max_threads();
//   membrane_chunk_vector chunks;
//   chunks.reserve(max_threads);
//   size_t chunk_size = members_.size() / max_threads;
//   membrane_iterator cur_iter = members_.begin();
//   for (int i = 0; i < max_threads - 1; ++i) {
//     membrane_iterator last_iter = cur_iter;
//     std::advance(cur_iter, chunk_size);
//     chunks.push_back(std::make_pair(last_iter, cur_iter));
//   }
//   chunks.push_back(std::make_pair(cur_iter, members_.end()));

// #pragma omp parallel shared(chunks)
//   {
// #pragma omp for
//     for (int i = 0; i < max_threads; ++i)
//       for (auto it = chunks[i].first; it != chunks[i].second; ++it)
//         it->UpdatePosition();
//   }
// #else
//   for (membrane_iterator it = members_.begin(); it != members_.end(); ++it)
//     it->UpdatePosition();
// #endif
//   if (sparams_.error_analysis) {
//     RunErrorAnalysis();
//   }
//   if (sparams_.spiral_init_flag && sparams_.spiral_number_fail_condition >= 0) {
//     int n_failed_spirals = 0;
//     for (auto it = members_.begin(); it != members_.end(); ++it) {
//       if (ABS(it->GetSpiralNumber()) < sparams_.spiral_number_fail_condition) {
//         // Failed spiral
//         n_failed_spirals++;
//       }
//     }
//     // If all spirals have failed, end simulation early
//     if (n_failed_spirals == n_members_) {
//       early_exit = true;
//     }
//   }
// }

// void MembraneSpecies::InitErrorAnalysis() {
//   std::string fname =
//       params_->run_name + "_membrane_" + sparams_.name + ".error.analysis";
//   error_file_.open(fname, std::ios::out);
//   if (!error_file_.is_open()) {
//     Logger::Error("Membrane error analysis file %s failed to open!",
//                   fname.c_str());
//   }
//   error_file_ << "n_dim delta length persistence_length min_bond_length "
//                  "driving_factor\n";
//   error_file_ << params_->n_dim << " " << params_->delta << " "
//               << sparams_.length << " " << sparams_.persistence_length << " "
//               << sparams_.min_bond_length << " " << sparams_.driving_factor
//               << "\n";
//   error_file_ << "steps_to_renorm\n";
// }

// void MembraneSpecies::RunErrorAnalysis() {
//   if (!error_file_.is_open()) {
//     Logger::Error("Error analysis file failed to open in MembraneSpecies");
//   }
//   std::vector<int> error_rates;
//   for (membrane_iterator it = members_.begin(); it != members_.end(); ++it)
//     it->GetErrorRates(error_rates);
//   for (auto it = error_rates.begin(); it != error_rates.end(); ++it) {
//     error_file_ << *it << "\n";
//   }
// }

// void MembraneSpecies::Reserve() {
//   int max_insert = GetNInsert();
//   if (packing_fraction_ > 0) {
//     double min_length = 2 * sparams_.min_bond_length;
//     double diameter = sparams_.diameter;
//     double min_vol = 0;
//     if (params_->n_dim == 2) {
//       min_vol = diameter * min_length + 0.25 * M_PI * diameter * diameter;
//     }
//     if (params_->n_dim == 3) {
//       min_vol = 0.25 * M_PI * diameter * diameter * min_length +
//                 1.0 / 6.0 * M_PI * diameter * diameter * diameter;
//     }
//     max_insert = (int)ceil(packing_fraction_ * space_->volume / min_vol);
//   }
//   members_.reserve(max_insert);
//   Logger::Debug("Reserving memory for %d members in MembraneSpecies",
//                 max_insert);
// }

// void MembraneSpecies::AddMember() {
//   Species::AddMember();
//   if (packing_fraction_ > 0) {
//     double vol = members_.back().GetVolume();
//     fill_volume_ += vol;
//     /* if we are still short on volume for the target packing fraction, then
//        request more members */
//     if (fill_volume_ < packing_fraction_ * space_->volume &&
//         members_.size() == sparams_.num) {
//       sparams_.num++;
//     }
//   }
// }

// void MembraneSpecies::PopMember() {
//   if (packing_fraction_ > 0) {
//     double vol = members_.back().GetVolume();
//     fill_volume_ -= vol;
//   }
//   Species::PopMember();
// }

// const double MembraneSpecies::GetSpecLength() const {
//   if (sparams_.dynamic_instability_flag) {
//     return 2 * sparams_.min_bond_length;
//   } else {
//     return 1.5 * sparams_.min_bond_length;
//   }
// }

// // Overloaded to catch error for membranes with receptor cover
// void MembraneSpecies::CalcPCPosition(int i, double s, double *pos) {
//   Logger::Error("Receptor PointCovers not set up for flexible membranes.");
// }

// void MembraneSpecies::LoadAnalysis() {
//   if (sparams_.msd_analysis) {
//     MembraneAnalysis *msd = new MSDAnalysis;
//     analysis_.push_back(msd);
//   }
//   if (sparams_.curvature_cluster_analysis) {
//     MembraneAnalysis *ccluster = new CurvatureClusterAnalysis;
//     analysis_.push_back(ccluster);
//   }
//   if (sparams_.spiral_analysis) {
//     MembraneAnalysis *spiral = new SpiralAnalysis;
//     analysis_.push_back(spiral);
//   }
//   if (sparams_.theta_analysis) {
//     MembraneAnalysis *theta = new AngleDistributionAnalysis;
//     analysis_.push_back(theta);
//   }
//   /* PolarOrderAnalysis before FlockingAnalysis, so we can calculate polar order
//      only once if both analyses are being done  */
//   if (sparams_.polar_order_analysis) {
//     MembraneAnalysis *polar_order = new PolarOrderAnalysis;
//     analysis_.push_back(polar_order);
//   }
//   if (sparams_.flocking_analysis) {
//     MembraneAnalysis *flock = new FlockingAnalysis;
//     analysis_.push_back(flock);
//   }
//   if (sparams_.orientation_corr_analysis) {
//     MembraneAnalysis *ocorr = new OrientationCorrelationAnalysis;
//     analysis_.push_back(ocorr);
//   }
//   if (sparams_.global_order_analysis) {
//     MembraneAnalysis *global_order = new GlobalOrderAnalysis;
//     analysis_.push_back(global_order);
//   }
//   if (sparams_.number_fluctuation_analysis) {
//     MembraneAnalysis *gnf = new GiantNumberFluctuationAnalysis;
//     analysis_.push_back(gnf);
//   }
//   if (sparams_.lp_analysis) {
//     MembraneAnalysis *mse2e = new EndToEndFluctuationAnalysis;
//     analysis_.push_back(mse2e);
//   }
//   if (sparams_.in_out_analysis) {
//     MembraneAnalysis *in_out = new IncomingOutgoingAngleAnalysis;
//     analysis_.push_back(in_out);
//   }
//   if (sparams_.crossing_analysis) {
//     MembraneAnalysis *crossing = new BarrierCrossingAnalysis;
//     analysis_.push_back(crossing);
//   }
// }
