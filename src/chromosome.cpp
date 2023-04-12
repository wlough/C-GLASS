#include "cglass/chromosome.hpp"
#include <iostream>

Chromosome::Chromosome(unsigned long seed) : Object(seed) {
  printf("NEW chromosome\n");
  SetSID(species_id::chromosome);
  sisters_.emplace_back(Chromatid(seed));
  sisters_.emplace_back(Chromatid(seed));
}

void Chromosome::Init(chromosome_parameters *sparams) {

  sparams_ = sparams;
  color_ = sparams->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = 0.0; //sparams->length;
  diameter_ = sparams->diameter;
  double pos[3] = {0, 0, 0};
  double u[3] = {0, 0, 0};
  rng_.RandomCoordinate(space_, pos, diameter_);
  rng_.RandomUnitVector(n_dim_, u);
  InsertAt(pos, u);
  // SF TODO change this to GetBodyFrame here and in Centromere()
  // Get vectors that define plane perpendicular to u
  double v[3] = {0, 0, 0}; // aux vector; defines plane of SPB together with w
  double w[3] = {0, 0, 0}; // aux vector; defines plane of SPB together with v
  // SF: not sure if this is general or b/c we initialize aligned to zhat
  double xhat[3] = {1.0, 0.0, 0.0};
  double yhat[3] = {0.0, 1.0, 0.0};
  if (1.0 - ABS(u[0]) > 1e-2) {
    cross_product(u, xhat, v, 3);
  } else {
    cross_product(u, yhat, v, 3);
  }
  double norm_factor{sqrt(1.0 / dot_product(3, v, v))};
  for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
    v[i_dim] *= norm_factor;
  }
  cross_product(u, v, w, 3);
  double pos_sis_one[3];
  double pos_sis_two[3];
  for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
    pos_sis_one[i_dim] = pos[i_dim] + v[i_dim] * diameter_ / 2.0;
    pos_sis_two[i_dim] = pos[i_dim] - v[i_dim] * diameter_ / 2.0;
  }
  for (auto &&sis : sisters_) {
    sis.Init(sparams);
  }
  sisters_[0].InsertAt(pos_sis_one, u);
  sisters_[1].InsertAt(pos_sis_two, u);
}

void Chromosome::UpdatePosition() {

  // ! SF TODO check interKC forces; currently unaccounted for

  // SF TODO check:
  // orientation -> v_
  // v_ -> orientation

  // First: calculate and add force/torques to chromatin sisters
  for (auto &&sis : sisters_) {
    sis.ZeroForce();
    sis.GetBodyFrame();
  }
  int n_dim{params_->n_dim};
  double k, r0, kth, kv;
  //SF TODO put these in parameter structure/yaml file
  k = 6.0;
  r0 = 4.0;
  kth = 450.0;
  kv = 450.0;
  double r[3], rhat[3];
  // Find the initial r vector r = r_A - r_B
  double rmag2 = 0.0;
  for (int i = 0; i < n_dim; ++i) {
    r[i] = sisters_[0].position_[i] - sisters_[1].position_[i];
    rmag2 += SQR(r[i]);
  }
  double rmag = sqrt(rmag2);
  for (int i = 0; i < n_dim; ++i) {
    rhat[i] = r[i] / rmag;
  }

  // Calculate angular stuff
  // double dotA = dot_product(n_dim, sisters_[0].orientation_, rhat);
  // double dotB = dot_product(n_dim, sisters_[1].orientation_, rhat);
  double dotA = dot_product(n_dim, sisters_[0].v_, rhat);
  double dotB = dot_product(n_dim, sisters_[1].v_, rhat);
  double thetaA = safe_acos(dotA);
  double thetaB = safe_acos(dotB);
  double sinthetaA = sin(thetaA);
  double sinthetaB = sin(thetaB);

  // Linear factor
  double linearfactor = -k * (rmag - r0);

  // theta A
  double thetaAoversinA = 0.0;
  double rhatcrossuA[3] = {0.0};
  double rhatcrossrhatcrossuA[3] = {0.0};
  if (sinthetaA != 0.0) {
    thetaAoversinA = thetaA / sinthetaA;
  }
  // cross_product(rhat, sisters_[0].orientation_, rhatcrossuA,
  cross_product(rhat, sisters_[0].v_, rhatcrossuA,
                n_dim); // cross product takes pointers, so can't just use again
  cross_product(rhat, rhatcrossuA, rhatcrossrhatcrossuA, n_dim);

  double thetaAfactor = 0.0;
  if (rmag > 0.0) {
    thetaAfactor = -kth / rmag * thetaAoversinA;
  }

  // theta B
  double thetaBoversinB = 0.0;
  double rhatcrossuB[3] = {0.0};
  double rhatcrossrhatcrossuB[3] = {0.0};
  if (sinthetaB != 0.0) {
    thetaBoversinB = thetaB / sinthetaB;
  }
  // cross_product(rhat, sisters_[1].orientation_, rhatcrossuB,
  cross_product(rhat, sisters_[1].v_, rhatcrossuB,
                n_dim); // cross product takes pointers, so can't just use again
  cross_product(rhat, rhatcrossuB, rhatcrossrhatcrossuB, n_dim);

  double thetaBfactor = 0.0;
  if (rmag > 0.0) {
    thetaBfactor = -kth / rmag * thetaBoversinB;
  }

  // update the forces
  double force[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < n_dim; ++i) {
    force_[i] =
        (linearfactor * rhat[i] + thetaAfactor * rhatcrossrhatcrossuA[i] +
         thetaBfactor * rhatcrossrhatcrossuB[i]);
    // f_chromosome[2 * ic + 1][i] -=
    //     (linearfactor * rhat[i] + thetaAfactor * rhatcrossrhatcrossuA[i] +
    //      thetaBfactor * rhatcrossrhatcrossuB[i]);
    // f_interkc[2 * ic][i] +=
    //     (linearfactor * rhat[i] + thetaAfactor * rhatcrossrhatcrossuA[i] +
    //      thetaBfactor * rhatcrossrhatcrossuB[i]);
    // f_interkc[2 * ic + 1][i] -=
    //     (linearfactor * rhat[i] + thetaAfactor * rhatcrossrhatcrossuA[i] +
    //      thetaBfactor * rhatcrossrhatcrossuB[i]);
  }
  sisters_[0].AddForce(force);
  sisters_[1].SubForce(force);

  // Other stuff required for torque
  double dotV =
      // dot_product(n_dim, sisters_[0].v_, sisters_[1].v_);
      dot_product(n_dim, sisters_[0].orientation_, sisters_[1].orientation_);
  double thetaV = safe_acos(dotV);
  double sinthetaV = sin(thetaV);
  double thetaVoversinV = 0.0;
  if (sinthetaV != 0.0) {
    thetaVoversinV = thetaV / sinthetaV;
  }
  double vhatAcrossvhatB[3] = {0.0};
  // cross_product(sisters_[0].v_, sisters_[1].v_,
  cross_product(sisters_[0].orientation_, sisters_[1].orientation_,
                vhatAcrossvhatB, n_dim);

  // torques
  double torque_one[3];
  double torque_two[3];
  for (int i = 0; i < n_dim; ++i) {
    torque_one[i] = (-kth * thetaAoversinA * rhatcrossuA[i] +
                     kv * thetaVoversinV * vhatAcrossvhatB[i]);
    torque_two[i] = (-kth * thetaBoversinB * rhatcrossuB[i] -
                     kv * thetaVoversinV * vhatAcrossvhatB[i]);
  }
  sisters_[0].AddTorque(torque_one);
  sisters_[1].AddTorque(torque_two);

  // Next: integrate to update their position (and add thermal noise)
  for (auto &&sis : sisters_) {
    sis.Integrate();
  }
  // Final: update position of centromere to be exactly in between two sisters
  double r_com[3];
  for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
    r_com[i_dim] =
        (sisters_[0].position_[i_dim] + sisters_[1].position_[i_dim]) / 2.0;
  }
  SetPosition(r_com);
}