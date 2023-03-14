#ifndef _CGLASS_CENTROSOME_H_
#define _CGLASS_CENTROSOME_H_

#include "anchor.hpp"
#include "filament.hpp"
// #include "receptor.hpp"
#include <armadillo>

struct AnchorSite {
  double k_{0.0};
  double kr_{0.0};
  double r0_{0.0};
  double pos_[3];
  double pos_rel_[3];
  double u_[3];
  double u_rel_[3];
  Filament *filament_;
  graph_struct g_;
};

class CentrosomeSpecies;

class Centrosome : public Object {
protected:
  bool zero_temperature_ = false;

  bool alignment_potential_, fixed_spacing_;
  int n_filaments_, n_filaments_min_, n_filaments_max_;
  double k_spring_, k_align_, spring_length_, anchor_distance_, gamma_trans_,
      gamma_rot_, diffusion_;

  double gamma_perp_ = 0;
  double noise_tr_ = 0.0;
  double noise_rot_ = 0.0;
  double diffusion_rot_ = 0.0;
  double body_frame_[6];

  static size_t i_spb_; // lol this is bad

  double r_[3];
  double u_[3];
  double v_[3];
  double w_[3];

  // std::vector<Filament> filaments_;
  std::vector<AnchorSite> anchors_;
  // Receptor anchor_;

  // filament_parameters *fparams_;
  centrosome_parameters *sparams_;

  // FIXME this is probably unnecessary and a sloppy middleman fix
  friend CentrosomeSpecies;

public:
protected:
  void ApplyForcesTorques();
  void ApplyBoundaryForces();
  void InsertCentrosome();
  void GenerateAnchorSites();
  void SetDiffusion();
  void Translate();
  void Rotate();
  void Integrate();
  void RandomizeAnchorPosition(int i_fil);

public:
  Centrosome(unsigned long seed);
  void Init(centrosome_parameters *sparams);
  // void InitFilaments(filament_parameters *fparams);
  void SetPosition(double const *const new_pos) {}
  // SF TODO only draws single tether
  void Draw(std::vector<graph_struct *> &graph_array) {
    for (auto &&anch : anchors_) {
      // std::copy(scaled_position_, scaled_position_ + 3, g_.r);
      double r_teth{0.0};
      for (int i = space_->n_periodic; i < n_dim_; ++i) {
        r_teth += SQR(anch.pos_[i] - anch.filament_->GetTailPosition()[i]);
        anch.g_.r[i] = anch.pos_[i];
      }
      r_teth = sqrt(r_teth);
      // std::copy(position_, position_+3, g_.r);
      std::copy(anch.u_, anch.u_ + 3, anch.g_.u);
      anch.g_.color = 1.8 * M_PI;
      // if (params_->graph_diameter > 0) {
      //   g_.diameter = params_->graph_diameter;
      // } else {
      //   g_.diameter = diameter_;
      // }
      anch.g_.diameter = 0.3;
      anch.g_.length = r_teth; //length_;
      anch.g_.draw = draw_;
      graph_array.push_back(&anch.g_);
    }
    // for (auto &&fila : filaments_) {
    //   fila.Draw(graph_array);
    // }
  }
  void GetInteractors(std::vector<Object *> &ixors) {
    // for (auto &&fila : filaments_) {
    //   fila.GetInteractors(ixors);
    // }
  }
  void UpdatePosition() {
    // for (auto &&fila : filaments_) {
    //   fila.UpdatePosition();
    // }
  }

  double GetR(int i_dim) { return r_[i_dim]; }
  double GetU(int i_dim) { return u_[i_dim]; }
  double GetV(int i_dim) { return v_[i_dim]; }
  double GetW(int i_dim) { return w_[i_dim]; }

  void ApplyInteractions() {

    for (auto &&anchor : anchors_) {
      double dr[3];
      // printf("boink\n");
      // Now that anchor positions are updated, apply force to filaments
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        dr[i_dim] =
            anchor.pos_[i_dim] - anchor.filament_->GetTailPosition()[i_dim];
        // printf("dr[%i] = %g\n", i_dim, dr[i_dim]);
        // (0.5 * anchor.filament_->GetLength() *
        //  anchor.filament_->GetOrientation()[i_dim]);
        // printf("dr[%i] = %g - %g = %g\n", i_dim, anchor.pos_[i_dim],
        //        anchor.filament_->GetTailPosition()[i_dim], dr[i_dim]);
      }
      double dr_mag2 = dot_product(params_->n_dim, dr, dr);
      // printf("dr_mag2 = %g\n\n", dr_mag2);
      double factor{
          dr_mag2 > 0.0 ? anchor.k_ * (1.0 - anchor.r0_ / sqrt(dr_mag2)) : 0.0};
      double f_spring[3];
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        f_spring[i_dim] = factor * dr[i_dim];
        // printf("f[%i] = %g\n", i_dim, f_spring[i_dim]);
      }
      // apply forces -- SF TODO: integrate this into Interacte() routine
      anchor.filament_->AddForceTail(f_spring);
      SubForce(f_spring);
      // anchor.filament_->UpdatePosition();
      // UpdatePosition();
      // UpdatePosition();
    }
  }

  void UpdatePosition(double *r_new, double *u_new, double *v_new,
                      double *w_new) {
    for (int i{0}; i < 3; i++) {
      r_[i] = position_[i] = r_new[i];
      u_[i] = orientation_[i] = u_new[i];
      v_[i] = v_new[i];
      w_[i] = w_new[i];
    }

    int i_anchor{0};
    for (auto &&anchor : anchors_) {
      /*
      double factor_u = dot_product(3, anchor.pos_rel_, u_);
      double factor_v = dot_product(3, anchor.pos_rel_, v_);
      double factor_w = dot_product(3, anchor.pos_rel_, w_);
      // Calculate new lab frame coordinate
      for (int i = 0; i < 3; ++i) {
        anchor.pos_[i] =
            r_[i] + factor_u * u_[i] + factor_v * v_[i] + factor_w * w_[i];
      }
      // Make sure that we properly attack the orientation of the anchor list
      // Orientations!
      double original_u[3] = {0.0};
      for (int i = 0; i < 3; ++i) {
        original_u[i] = anchor.u_[i];
        anchor.u_[i] = anchor.u_rel_[0] * u_[i] + anchor.u_rel_[1] * v_[i] +
                       anchor.u_rel_[2] * w_[i];
      }
      */
      // Centrosome surface is 2-D,so we can define position w/ R and phi
      // For now, distribute them uniformly along circumference of a circle
      double r_on_spb = 1.25;
      double delta_phi = 2.0 * M_PI / double(n_filaments_);
      double base_pos[3];
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        base_pos[i_dim] = GetPosition()[i_dim] +
                          r_on_spb * sin(i_anchor * delta_phi) * v_[i_dim] +
                          r_on_spb * cos(i_anchor * delta_phi) * w_[i_dim];
      }
      i_anchor++;
      double r_teth{0.0};
      double u_new[3] = {0, 0, 0};
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        // anchor.u_[i_dim] = GetOrientation()[i_dim];
        anchor.pos_[i_dim] = base_pos[i_dim] + 0.5 * anchor.r0_ * u_[i_dim];
        u_new[i_dim] =
            anchor.pos_[i_dim] - anchor.filament_->GetTailPosition()[i_dim];
        r_teth += SQR(u_new[i_dim]);
      }
      r_teth = sqrt(r_teth);
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        u_new[i_dim] = u_new[i_dim] / r_teth;
        anchor.u_[i_dim] = u_new[i_dim];
      }
      // Create a check if the position isn't working correctly
      for (int i = 0; i < 3; ++i) {
        if (std::isnan(anchor.pos_[i])) {
          printf("woah buddy\n");
          exit(1);
          /*
          // Print out the information of the centrosome and exit
          std::cerr << "Encountered an error in SPB update code, printing "
                       "information then exiting\n";
          std::cerr << "  step: " << properties->i_current_step << std::endl;
          print_centrosomes_information(idx, &(properties->anchors));
          std::cerr << "  Force = " << force.as_row();
          std::cerr << "  Torque = " << torque.as_row();
          exit(1);
          */
        }
      }
      /*
      // Now that anchor positions are updated, apply force to filaments
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        dr[i_dim] = anchor.filament_->GetPosition()[i_dim] -
                    anchor.pos_[i_dim] -
                    (0.5 * anchor.filament_->GetLength() *
                     anchor.filament_->GetOrientation()[i_dim]);
      }
      double dr_mag2 = dot_product(params_->n_dim, dr, dr);
      double factor = anchor.k_ * (1.0 - anchor.r0_ / sqrt(dr_mag2));
      double f_spring[3];
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        f_spring[i_dim] = factor * dr[i_dim];
        // printf("f[%i] = %g\n", i_dim, f_spring[i_dim]);
      }
      // apply forces -- SF TODO: integrate this into Interacte() routine
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      }
      anchor.filament_->SubForce(f_spring);
      AddForce(f_spring);
      // anchor.filament_->UpdatePosition();
      // UpdatePosition();
      */
    }
    // UpdatePosition();
  }
};

#endif // _CGLASS_CENTROSOME_H_
