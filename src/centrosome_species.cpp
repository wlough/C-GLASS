#include "cglass/centrosome_species.hpp"

CentrosomeSpecies::CentrosomeSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::centrosome);
  printf("NEW centrosome!\n");
}

void CentrosomeSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  if (GetNInsert() <= 0) {
    return;
  }
  printf("Initializing centrosome\n");
  int n_anchors{GetNInsert()};
  // TODO make sure you delete these in deconstructor
  // All of the new rotation information for this class
  q_ = new arma::vec4[n_anchors];
  A_ = new arma::mat33[n_anchors];
  mu_tb_ = new arma::mat33[n_anchors];
  mu_rb_ = new arma::mat33[n_anchors];
  sqrt_mu_tb_ = new arma::mat33[n_anchors];
  sqrt_mu_rb_ = new arma::mat33[n_anchors];
  bsub_ = new arma::mat::fixed<4, 3>[n_anchors];

  double attach_diameter[n_anchors];

  // initialize geometry stuff lol
  for (int i_spb{0}; i_spb < n_anchors; i_spb++) {
    double spb_diffusion = 0.256;
    attach_diameter[i_spb] = 6.5;
    // Now that we have the diffusion scale, we need to encode the proper SPB
    // translation and rotation matrices into the matrices (including sqrt)
    // for proper rigid body dynamics
    arma::mat33 translation;
    arma::mat33 rotation;
    translation.eye();
    rotation.eye();

    //   parameters->spb_diffusion * diffusion_scale_spbs;

    translation(0, 0) = 0.5 * spb_diffusion;
    translation(1, 1) = spb_diffusion;
    translation(2, 2) = spb_diffusion;
    rotation(0, 0) = spb_diffusion / SQR(0.5 * attach_diameter[i_spb]);
    rotation(1, 1) = 0.01 * spb_diffusion / SQR(0.5 * attach_diameter[i_spb]);
    rotation(2, 2) = 0.01 * spb_diffusion / SQR(0.5 * attach_diameter[i_spb]);

    // Set the mobility matrices (assuming they don't change) and the sqrt of same
    mu_tb_[i_spb] = translation;
    mu_rb_[i_spb] = rotation;
    sqrt_mu_tb_[i_spb] = arma::sqrtmat_sympd(translation);
    sqrt_mu_rb_[i_spb] = arma::sqrtmat_sympd(rotation);
  }
}

void CentrosomeSpecies::Reserve() {
  SpeciesBase::Reserve();
  members_.reserve(GetNInsert());
}

void CentrosomeSpecies::PopMember() { Species::PopMember(); }
void CentrosomeSpecies::AddMember() {
  Species::AddMember();
  double u[3] = {members_.back().GetOrientation()[0],
                 members_.back().GetOrientation()[1],
                 members_.back().GetOrientation()[2]};
  int i_spb{members_.size() - 1};
  // Create the rotatio matrix from body-space coordinate transformation
  CreateRotationMatrixBodySpace(A_[i_spb], u, members_.back().v_,
                                members_.back().w_);

  // Create the associated quaternion describing the initial orientation of the
  // rigid body
  QuaternionFromRotationMatrix(q_[i_spb], A_[i_spb]);
}

void CentrosomeSpecies::UpdatePositions() {

  const double thermal_diff_scale = 10;
  double delta_t = params_->delta;
  double kT = 1.0;
  const arma::vec3 nx = {1.0, 0.0, 0.0};
  const arma::vec3 ny = {0.0, 1.0, 0.0};
  const arma::vec3 nz = {0.0, 0.0, 1.0};
  for (int idx{0}; idx < members_.size(); idx++) {

    Centrosome *centro{&members_[idx]};

    // Set up armadillo versions of all of our information of interest
    arma::vec3 r = {centro->GetPosition()[0], centro->GetPosition()[1],
                    centro->GetPosition()[2]};
    arma::vec3 u = {centro->GetOrientation()[0], centro->GetOrientation()[1],
                    centro->GetOrientation()[2]};
    arma::vec3 v = {centro->v_[0], centro->v_[1], centro->v_[2]};
    arma::vec3 w = {centro->w_[0], centro->w_[1], centro->w_[2]};

    arma::vec3 force = {centro->GetForce()[0], centro->GetForce()[1],
                        centro->GetForce()[2]};
    arma::vec3 torque = {centro->GetTorque()[0], centro->GetTorque()[1],
                         centro->GetTorque()[2]};

    // Deterministic translation
    arma::vec3 dr = (((A_[idx] * mu_tb_[idx]) * A_[idx].t()) * force) * delta_t;
    // Random thermal translation motion
    arma::vec3 theta = {thermal_diff_scale * rng_.RandomNormal(1.0),
                        thermal_diff_scale * rng_.RandomNormal(1.0),
                        thermal_diff_scale * rng_.RandomNormal(1.0)};
    dr +=
        ((A_[idx] * sqrt_mu_tb_[idx]) * theta) * std::sqrt(2.0 * kT * delta_t);

    // for (int i{0}; i < 3; i++) {
    //   printf("%g\n", dr(i));
    // }
    // Now handle the rotation part
    bsub_[idx](0, 0) = -q_[idx][1];
    bsub_[idx](0, 1) = -q_[idx][2];
    bsub_[idx](0, 2) = -q_[idx][3];

    bsub_[idx](1, 0) = q_[idx][0];
    bsub_[idx](1, 1) = -q_[idx][3];
    bsub_[idx](1, 2) = q_[idx][2];

    bsub_[idx](2, 0) = q_[idx][3];
    bsub_[idx](2, 1) = q_[idx][0];
    bsub_[idx](2, 2) = -q_[idx][1];

    bsub_[idx](3, 0) = -q_[idx][2];
    bsub_[idx](3, 1) = q_[idx][1];
    bsub_[idx](3, 2) = q_[idx][0];

    // Deterministic rotation part
    arma::vec4 dq =
        (((bsub_[idx] * mu_rb_[idx]) * A_[idx].t()) * torque) * delta_t;

    // Random reorientation
    arma::vec3 theta_rotation = {thermal_diff_scale * rng_.RandomNormal(1.0),
                                 thermal_diff_scale * rng_.RandomNormal(1.0),
                                 thermal_diff_scale * rng_.RandomNormal(1.0)};
    dq += ((bsub_[idx] * sqrt_mu_rb_[idx]) * theta_rotation) *
          std::sqrt(2.0 * kT * delta_t);

    // Compute the lagrange correction
    arma::vec4 qtilde = q_[idx] + dq;
    // Compute the corrected quaternion

    // double lambdaq = compute_lagrange_correction(q[idx], qtilde);
    double lambdaq = ComputeLagrangeCorrection(q_[idx], qtilde);
    if (std::isnan(lambdaq)) {
      std::cerr << " quaternion correction error: " << lambdaq << std::endl;
    }

    // Actually update the information
    r += dr;
    q_[idx] = qtilde + lambdaq * q_[idx];

    // Update the rotation matrix and orientation vectors
    // rotation_matrix_from_quaternion(A[idx], q[idx]);
    RotationMatrixFromQuaternion(A_[idx], q_[idx]);
    u = A_[idx] * nx;
    v = A_[idx] * ny;
    w = A_[idx] * nz;

    double r_new[3];
    double u_new[3];
    double v_new[3];
    double w_new[3];

    // Write back into the anchor structure the proper double variables
    for (int i = 0; i < 3; ++i) {
      r_new[i] = r(i);
      printf("r[%i] = %g\n", i, r_new[i]);
      u_new[i] = u(i);
      v_new[i] = v(i);
      w_new[i] = w(i);
    }

    centro->SetPosition(r_new);
    // printf("BONK\n");
    centro->SetOrientation(u_new);
    std::copy(v_new, v_new + 3, centro->v_);
    std::copy(w_new, w_new + 3, centro->w_);

    // Now we have to update the positions of the MTs on the anchor
    /*
    for (al_list::iterator p = properties->anchors.anchor_list[idx].begin();
         p < properties->anchors.anchor_list[idx].end(); p++) {
      double factor_u = dot_product(3, p->pos_rel, u_anchor[idx]);
      double factor_v = dot_product(3, p->pos_rel, v_anchor[idx]);
      double factor_w = dot_product(3, p->pos_rel, w_anchor[idx]);

      // Calculate new lab frame coordinate
      for (int i = 0; i < 3; ++i) {
        p->pos[i] = r_anchor[idx][i] + factor_u * u_anchor[idx][i] +
                    factor_v * v_anchor[idx][i] + factor_w * w_anchor[idx][i];
      }
      // Make sure that we properly attack the orientation of the anchor list
      // Orientations!
      double original_u[3] = {0.0};
      for (int i = 0; i < 3; ++i) {
        original_u[i] = p->u[i];
        p->u[i] = p->u_rel[0] * u_anchor[idx][i] +
                  p->u_rel[1] * v_anchor[idx][i] +
                  p->u_rel[2] * w_anchor[idx][i];
      }
      // Create a check if the position isn't working correctly
      for (int i = 0; i < 3; ++i) {
        if (std::isnan(p->pos[i])) {
          // Print out the information of the centrosome and exit
          std::cerr << "Encountered an error in SPB update code, printing "
                       "information then exiting\n";
          std::cerr << "  step: " << properties->i_current_step << std::endl;
          print_centrosomes_information(idx, &(properties->anchors));
          std::cerr << "  Force = " << force.as_row();
          std::cerr << "  Torque = " << torque.as_row();
          exit(1);
        }
      }
    }
    */
  }
}

double CentrosomeSpecies::ComputeLagrangeCorrection(const arma::vec4 &q,
                                                    const arma::vec4 &qtilde) {
  double qqtilde = arma::dot(q, qtilde);
  double qtilde2 = arma::dot(qtilde, qtilde);

  // Solve the quadratic formulat for what is correct. Doesn't matter which
  // solution we take, as the rotations are equivalent
  double d = (2.0 * qqtilde) * (2.0 * qqtilde) - 4.0 * (qtilde2 - 1.0);
  double lambdaq = (-2.0 * qqtilde + std::sqrt(d)) / 2.0;
  return lambdaq;
}

void CentrosomeSpecies::CreateRotationMatrixBodySpace(arma::mat33 &A, double *u,
                                                      double *v, double *w) {
  arma::vec3 u_ = {u[0], u[1], u[2]};
  arma::vec3 v_ = {v[0], v[1], v[2]};
  arma::vec3 w_ = {w[0], w[1], w[2]};
  A.col(0) = u_;
  A.col(1) = v_;
  A.col(2) = w_;
}

void CentrosomeSpecies::QuaternionFromRotationMatrix(arma::vec4 &q,
                                                     const arma::mat33 R) {
  if ((R(1, 1) >= -R(2, 2)) && (R(0, 0) >= -R(1, 1)) && (R(0, 0) >= -R(2, 2))) {
    q[0] = std::sqrt(1 + R(0, 0) + R(1, 1) + R(2, 2));
    q[1] = (R(2, 1) - R(1, 2)) / std::sqrt(1 + R(0, 0) + R(1, 1) + R(2, 2));
    q[2] = (R(0, 2) - R(2, 0)) / std::sqrt(1 + R(0, 0) + R(1, 1) + R(2, 2));
    q[3] = (R(1, 0) - R(0, 1)) / std::sqrt(1 + R(0, 0) + R(1, 1) + R(2, 2));
  } else if ((R(1, 1) <= -R(2, 2)) && (R(0, 0) >= R(1, 1)) &&
             (R(0, 0) >= R(2, 2))) {
    q[0] = (R(2, 1) - R(1, 2)) / std::sqrt(1 + R(0, 0) - R(1, 1) - R(2, 2));
    q[1] = std::sqrt(1 + R(0, 0) - R(1, 1) - R(2, 2));
    q[2] = (R(1, 0) + R(0, 1)) / std::sqrt(1 + R(0, 0) - R(1, 1) - R(2, 2));
    q[3] = (R(2, 0) + R(0, 2)) / std::sqrt(1 + R(0, 0) - R(1, 1) - R(2, 2));
  } else if ((R(1, 1) >= R(2, 2)) && (R(0, 0) <= R(1, 1)) &&
             (R(0, 0) <= -R(2, 2))) {
    q[0] = (R(0, 2) - R(2, 0)) / std::sqrt(1 - R(0, 0) + R(1, 1) - R(2, 2));
    q[1] = (R(1, 0) + R(0, 1)) / std::sqrt(1 - R(0, 0) + R(1, 1) - R(2, 2));
    q[2] = std::sqrt(1 - R(0, 0) + R(1, 1) - R(2, 2));
    q[0] = (R(2, 1) + R(1, 2)) / std::sqrt(1 - R(0, 0) + R(1, 1) - R(2, 2));
  } else if ((R(1, 1) <= R(2, 2)) && (R(0, 0) <= -R(1, 1)) &&
             (R(0, 0) <= R(2, 2))) {
    q[0] = (R(1, 0) - R(0, 1)) / std::sqrt(1 - R(0, 0) - R(1, 1) + R(2, 2));
    q[1] = (R(2, 0) + R(0, 2)) / std::sqrt(1 - R(0, 0) - R(1, 1) + R(2, 2));
    q[2] = (R(2, 1) + R(1, 2)) / std::sqrt(1 - R(0, 0) - R(1, 1) + R(2, 2));
    q[3] = std::sqrt(1 - R(0, 0) - R(1, 1) + R(2, 2));
  } else {
    std::cerr
        << "ERROR: quaterion_from_rotation_matrix did not work correctly!\n";
    exit(1);
  }

  // Have to multiply by 1/2 on all cases
  for (auto i = 0; i < 4; ++i) {
    q[i] = 0.5 * q[i];
  }
}

void CentrosomeSpecies::RotationMatrixFromQuaternion(arma::mat33 &R,
                                                     const arma::vec4 &q) {
  R(0, 0) = 1. - 2. * (q[2] * q[2] + q[3] * q[3]);
  R(0, 1) = 2. * (q[1] * q[2] - q[3] * q[0]);
  R(0, 2) = 2. * (q[1] * q[3] + q[2] * q[0]);

  R(1, 0) = 2. * (q[1] * q[2] + q[3] * q[0]);
  R(1, 1) = 1 - 2. * (q[1] * q[1] + q[3] * q[3]);
  R(1, 2) = 2. * (q[2] * q[3] - q[1] * q[0]);

  R(2, 0) = 2. * (q[1] * q[3] - q[2] * q[0]);
  R(2, 1) = 2. * (q[2] * q[3] + q[1] * q[0]);
  R(2, 2) = 1 - 2. * (q[1] * q[1] + q[2] * q[2]);
}

void CentrosomeSpecies::GetInteractors(std::vector<Object *> &ixors) {
  for (auto &&centro : members_) {
    centro.GetInteractors(ixors);
  }
}