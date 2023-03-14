#include "cglass/centrosome_species.hpp"
#include "cglass/crosslink_species.hpp"
#include "cglass/filament_species.hpp"

// TODO: add overlap detection, wca potential between centros

CentrosomeSpecies::CentrosomeSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::centrosome);
}

void CentrosomeSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  if (GetNInsert() <= 0) {
    return;
  }
  // std::string filament_species_name = sparams_.filament_species_name;
  // printf("filament species name: %s\n", filament_species_name.c_str());
  // // First, store filament parameters to later pass to each individual centrosome
  // species_base_parameters *sparams{parser.GetNewSpeciesParameters(
  //     species_id::filament, filament_species_name)};
  // fparams_ = *dynamic_cast<filament_parameters *>(sparams);
  // delete sparams;

  // NAB uses 'anchor' to refer to SPBs. Anchors are unique objects in C-GLASS tho (FIXME)
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
    double spb_diffusion = 2.56; //  x10
    //   parameters->spb_diffusion * diffusion_scale_spbs;
    attach_diameter[i_spb] = 6.5;
    // Now that we have the diffusion scale, we need to encode the proper SPB
    // translation and rotation matrices into the matrices (including sqrt)
    // for proper rigid body dynamics
    arma::mat33 translation;
    arma::mat33 rotation;
    translation.eye();
    rotation.eye();
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
  // Centrosome::i_spb_ = 0;
  Species::AddMember();
  // Not worth it to handle filaments entirely independently
  // Use FilamentSpecies and link at initialization
  // members_.back().InitFilaments(&fparams_);
  int i_spb{(int)members_.size() - 1};
  double u[3] = {members_.back().GetU(0), members_.back().GetU(1),
                 members_.back().GetU(2)};
  // Create the rotatio matrix from body-space coordinate transformation
  CreateRotationMatrixBodySpace(A_[i_spb], u, members_.back().v_,
                                members_.back().w_);
  // Create the associated quaternion describing initial orientation of the rigid body
  QuaternionFromRotationMatrix(q_[i_spb], A_[i_spb]);
}

void CentrosomeSpecies::AnchorFilaments(FilamentSpecies *filas,
                                        CrosslinkSpecies *teths) {

  int n_filaments{filas->GetNMembers()};
  if (n_filaments != GetNMembers() * sparams_.num_filaments_ea) {
    printf("oh no\n");
    exit(1);
  }

  // need to:
  // Create tether (AddMember())
  // Insert it proper
  // Bind one head to SPB
  // Insert filament
  // Bind other to filament
  // printf("no 1\n");
  // printf("FOUND SPB tethers - '%s'\n", teths->GetSpeciesName().c_str());
  // printf("no 2\n");

  // SF TODO: generalize this for more than 1 SPB / anchor lol
  size_t i_spb{0};
  for (int i_fila{0}; i_fila < filas->GetNMembers(); i_fila++) {
    Filament *fil{dynamic_cast<Filament *>(filas->GetMember(i_fila))};
    size_t i_anchor{i_fila < sparams_.num_filaments_ea
                        ? i_fila
                        : i_fila - sparams_.num_filaments_ea};
    if (i_fila == sparams_.num_filaments_ea) {
      i_spb++;
    }
    // FIXME for discrete anchor site locations
    const double *const anchor_u = members_[i_spb].anchors_[i_anchor].u_;
    const double *const anchor_pos = members_[i_spb].anchors_[i_anchor].pos_;
    double new_pos[3] = {0, 0, 0};
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      new_pos[i_dim] =
          anchor_pos[i_dim] + 0.5 * fil->GetLength() * anchor_u[i_dim];
      printf("%g = %g + %g\n", new_pos[i_dim], anchor_pos[i_dim],
             0.5 * fil->GetLength() * anchor_u[i_dim]);
    }
    fil->InsertAt(new_pos, anchor_u);
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      printf("(%g)\n", fil->GetTailPosition()[i_dim]);
    }
    // teths->AddMember();
    // teths->GetMember(0)->InsertAt(new_pos, anchor_u);
    members_[i_spb].anchors_[i_anchor].filament_ = fil;
    // fil->SetOrientation(u);
  }

  // int i_fila{0};
  // for (auto &&fil : filas->GetMembers()) {
  //   double u[3] = {0, 0, 1}; // Initially align with z
  //   fil.SetOrientation(u);
  //   printf("FILA #%i\n", i_fila++);
  // }

  // UPDATE POS AND U OF ANCHORSITE
  // UPDATE POS AND U OF FILAMENT
  // UPDATE REL_POS AND REL_U OF ANCHORSITE
}

void CentrosomeSpecies::UpdatePositions() {

  double kT = 1.0;
  const double thermal_diff_scale = 1;
  double sys_radius = space_->radius;
  double delta_t = params_->delta;
  const arma::vec3 nx = {1.0, 0.0, 0.0};
  const arma::vec3 ny = {0.0, 1.0, 0.0};
  const arma::vec3 nz = {0.0, 0.0, 1.0};
  for (int idx{0}; idx < members_.size(); idx++) {
    Centrosome *centro{&members_[idx]};
    // Construct the amount (and direction) the SPB is outside the radius
    // If outside the nucleus, negative
    double rmag = 0.0;
    for (int i = 0; i < 3; ++i) {
      rmag += SQR(centro->GetR(i));
    }
    rmag = std::sqrt(rmag);
    double rhat[3] = {0.0};
    for (int i = 0; i < 3; ++i) {
      rhat[i] = centro->GetR(i) / rmag;
    }
    // Construct the force vector
    double forcevec[3] = {0.0};

    /*
    // If we are using a harmonic force
    if (properties->anchors.spb_confinement_type == 0) {
      double k = properties->anchors.centrosome_confinement_radial_k_;

      // Guard against pathological behavior by capping the force at some level
      double factor = -k * (rmag - sys_radius);
      double fcutoff = 0.1 / parameters->delta /
                       MAX(properties->anchors.mu_tb_[idx](0, 0),
                           properties->anchors.mu_tb_[idx](1, 1));
      if (factor < 0.0) {
        if (-1.0 * factor > fcutoff) {
          factor = -1.0 * fcutoff;
          std::cerr << " *** Force exceeded fcutoff "
                       "spb_confinement_potential_harmonic_bd radial ***\n";
        }
      } else {
        if (factor > fcutoff) {
          factor = fcutoff;
          std::cerr << " *** Force exceeded fcutoff "
                       "spb_confinement_potential_harmonic_bd radial ***\n";
        }
      }

      for (int i = 0; i < 3; ++i) {
        forcevec[i] = factor * rhat[i];
      }
    } else if (properties->anchors.spb_confinement_type == 1) {
        */

    // properties->anchors.centrosome_confinement_f0_;
    double f0 = 103.406326034025;
    double delta_r = ABS(sys_radius - rmag);
    double ne_ratio{24.61538};
    double factor = CalcNonMonotonicWallForce(ne_ratio, f0, delta_r);
    // Check the sign of the force, want to move outwards if we are in the nucleoplasm
    if (rmag > sys_radius) {
      for (int i = 0; i < 3; ++i) {
        forcevec[i] = -1.0 * factor * rhat[i];
      }
    } else {
      for (int i = 0; i < 3; ++i) {
        forcevec[i] = 1.0 * factor * rhat[i];
      }
    }
    // }
    //Based on the potential, calculate the torques on the system
    double uhat_dot_rhat = dot_product(3, rhat, centro->GetOrientation());
    double uhat_cross_rhat[3] = {0.0};
    cross_product(centro->GetOrientation(), rhat, uhat_cross_rhat, 3);
    double kr = 9999.0; //properties->anchors.centrosome_confinement_angular_k_;
    // The torque is thusly
    double torquevec[3] = {0.0};
    for (int i = 0; i < 3; ++i) {
      torquevec[i] = -kr * (uhat_dot_rhat + 1.0) * uhat_cross_rhat[i];
    }
    // Add the contribution to the total forces
    centro->AddForce(forcevec);
    centro->AddTorque(torquevec);

    /* BELOW CODE IS FOR 'FREE' SPB, I.E., NOT CONFINED TO MEMBRANE */
    // Set up armadillo versions of all of our information of interest
    arma::vec3 r = {centro->GetR(0), centro->GetR(1), centro->GetR(2)};
    arma::vec3 u = {centro->GetU(0), centro->GetU(1), centro->GetU(2)};
    arma::vec3 v = {centro->GetV(0), centro->GetV(1), centro->GetV(2)};
    arma::vec3 w = {centro->GetW(0), centro->GetW(1), centro->GetW(2)};
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
    double lambdaq = ComputeLagrangeCorrection(q_[idx], qtilde);
    if (std::isnan(lambdaq)) {
      std::cerr << " quaternion correction error: " << lambdaq << std::endl;
    }
    // Actually update the information
    r += dr;
    q_[idx] = qtilde + lambdaq * q_[idx];
    // Update the rotation matrix and orientation vectors
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
      u_new[i] = u(i);
      v_new[i] = v(i);
      w_new[i] = w(i);
    }

    centro->UpdatePosition(r_new, u_new, v_new, w_new);

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

double CentrosomeSpecies::CalcNonMonotonicWallForce(double ne_ratio, double f0,
                                                    double delta_r) {
  //Constants from paper
  double a1 = 0.746;
  double a2 = 0.726;
  double alpha1 = 0.347;
  double alpha2 = 3.691;
  double a = a1 * a2;
  double b = sqrt(ne_ratio);
  double c = alpha1 + alpha2;
  // Calculate location of maximum force
  double xm = b * (((7.0 / 4.0) * M_PI) - c); //value of dr that gives max force
  // Calculate wall force
  // Exponential prefactor added to make forces conform to boundary conditions.
  // Characteristic length chosen to be 1/30th the value of delta_r that gives max force
  double factor =
      (1 - exp(-delta_r * 30 / xm)) *
      (2.0 * f0 * a * exp(-1.0 * delta_r / b) * cos((delta_r / b) + c) + f0);
  return factor;
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
