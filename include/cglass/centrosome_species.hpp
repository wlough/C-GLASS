#ifndef _CGLASS_CENTROSOME_SPECIES_H_
#define _CGLASS_CENTROSOME_SPECIES_H_

#include "centrosome.hpp"
#include "species.hpp"

class FilamentSpecies;
class CrosslinkSpecies;
class CentrosomeSpecies : public Species<Centrosome, species_id::centrosome> {

private:
  arma::vec4 *q_;  // Quaternion describing current rotation of the SPB
  arma::mat33 *A_; // Rotation matrix describing SPB body-frame to fixed-frame
  arma::mat33 *mu_tb_;           // Translation mobility matrix
  arma::mat33 *mu_rb_;           // Rotation mobility matrix
  arma::mat33 *sqrt_mu_tb_;      // Sqrt translation mobility matrix
  arma::mat33 *sqrt_mu_rb_;      // Sqrt rotation mobility matrix
  arma::mat::fixed<4, 3> *bsub_; // Transform matrix for quaternion rotation

  // filament_parameters fparams_;

protected:
  bool midstep_;

public:
  CentrosomeSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void PopMember();
  void AddMember();
  void Reserve();
  std::string GetFilamentSpeciesName() {
    return sparams_.filament_species_name;
  }
  std::string GetCrosslinkSpeciesName() {
    return sparams_.crosslink_species_name;
  }
  void AnchorFilaments(FilamentSpecies *filas, CrosslinkSpecies *teths);

  void UpdatePositions();

  void ApplyInteractions() {

    /*
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
          // printf("f[%i] = %g\n", i, forcevec[i]);
        }
      }
      // }
      //Based on the potential, calculate the torques on the system
      double uhat_dot_rhat = dot_product(3, rhat, centro->GetOrientation());
      double uhat_cross_rhat[3] = {0.0};
      cross_product(centro->GetOrientation(), rhat, uhat_cross_rhat, 3);
      double kr =
          9999.0; //properties->anchors.centrosome_confinement_angular_k_;
      // The torque is thusly
      double torquevec[3] = {0.0};
      for (int i = 0; i < 3; ++i) {
        torquevec[i] = -kr * (uhat_dot_rhat + 1.0) * uhat_cross_rhat[i];
      }
      // Add the contribution to the total forces
      centro->AddForce(forcevec);
      centro->AddTorque(torquevec);
    }
    */
    for (auto &&centro : members_) {
      centro.ApplyInteractions();
    }
  }

  double CalcNonMonotonicWallForce(double ne_ratio, double f0, double delta_r);

  double ComputeLagrangeCorrection(const arma::vec4 &q,
                                   const arma::vec4 &qtilde);

  void CreateRotationMatrixBodySpace(arma::mat33 &A, double *u, double *v,
                                     double *w);

  void QuaternionFromRotationMatrix(arma::vec4 &q, const arma::mat33 R);
  void RotationMatrixFromQuaternion(arma::mat33 &R, const arma::vec4 &q);

  void GetInteractors(std::vector<Object *> &ixors);
  void GetLastInteractors(std::vector<Object *> &ixors) {}
};

#endif
