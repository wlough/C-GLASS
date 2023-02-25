#ifndef _CGLASS_CENTROSOME_SPECIES_H_
#define _CGLASS_CENTROSOME_SPECIES_H_

#include "centrosome.hpp"
#include "species.hpp"

class CentrosomeSpecies : public Species<Centrosome, species_id::centrosome> {

private:
  arma::vec4 *q_; ///< Quaternion describing current rotation of the SPB
  arma::mat33 *
      A_; ///< Rotation matrix describing the body-frame to fixed-frame of the SPB
  arma::mat33 *mu_tb_;      ///< Translation mobility matrix
  arma::mat33 *mu_rb_;      ///< Rotation mobility matrix
  arma::mat33 *sqrt_mu_tb_; ///< Sqrt translation mobility matrix
  arma::mat33 *sqrt_mu_rb_; ///< Sqrt rotation mobility matrix
  arma::mat::fixed<4, 3>
      *bsub_; ///< Transformation matrix needed for rotation of quaternions
protected:
  bool midstep_;

public:
  CentrosomeSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void PopMember();
  void AddMember();

  void Reserve();

  void UpdatePositions();

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
