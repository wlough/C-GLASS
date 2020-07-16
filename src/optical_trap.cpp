/**
 * @author      : Adam Lamson (adam.lamson@colorado.edu)
 * @file        : optical_trap
 * @created     : Friday Jul 03, 2020 14:49:10 MDT
 */

#include "cglass/optical_trap.hpp"

OpticalTrap::OpticalTrap(unsigned long seed) : Object(seed), bead_(seed) {
  SetSID(species_id::optical_trap);
}

/*! \brief Initialize optical trap object
 *
 *  Detailed description
 *
 * \param *sparams Parameter description
 * \return Return parameter description
 */
void OpticalTrap::Init(optical_trap_parameters *sparams) {
  sparams_ = sparams;
  trap_spring_ = sparams_->trap_spring;
  diameter_ = sparams_->trap_diameter;
  color_ = sparams_->trap_color;

  bead_.SetSpace(space_);
  bead_.SetNDim(n_dim_);
  bead_.SetDiameter(sparams_->bead_diameter);
  bead_.SetColor(sparams_->bead_color, draw_type::fixed);
  //bead_=Object(0)
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
}

/*! \brief Set position of optical trap relative to obj minus end
 *
 *  Detailed description
 *
 * \param Object *obj Parameter description
 * \return Return parameter description
 */
void OpticalTrap::InsertAndAttach(Object *obj) {
  attach_obj_ = obj;
  SetMeshID(-1);
  AttachObjRelLambda(0); //TODO attach to minus end for now
  const double *obj_pos = obj->GetPosition();
  const double *obj_orient = obj->GetOrientation();
  const double obj_len = obj->GetLength() * .5;

  for (int i = 0; i < 3; ++i) {
    position_[i] = obj_pos[i] - obj_orient[i] * obj_len;
  }

  UpdatePeriodic();
  bead_.SetPosition(position_);
  bead_.UpdatePeriodic();

  return;
}

/*! \brief Apply force from optical trap to bond
 *
 * \return void
 */
void OpticalTrap::ApplyOpticalTrapForce() {
  if (!bond_) {
    Logger::Error("Anchor attempted to apply forces to nullptr bond");
  }
}

void OpticalTrap::AttachObjRelLambda(double lambda) {
  //if (attach_obj_->GetType() != +obj_type::bond) {
  //  Logger::Error(
  //      "Optical traps bound to non-bond objects not yet implemented in "
  //      "AttachObjLambda.");
  //}

  if (attach_obj_->GetType() == +obj_type::bond) {
    bond_ = dynamic_cast<Bond *>(attach_obj_);
    bond_length_ = bond_->GetLength();
    if (bond_ == nullptr) {
      Logger::Error(
          "Object ptr passed to optical trap was not referencing a bond!");
    }
    mesh_ = dynamic_cast<Mesh *>(bond_->GetMeshPtr());
    if (mesh_ != nullptr) {
      mesh_n_bonds_ = mesh_->GetNBonds();
      mesh_length_ = mesh_->GetTrueLength();
      mesh_lambda_ = mesh_length_ * lambda;
      bond_lambda_ = mesh_lambda_ - bond_->GetMeshLambda();
      SetMeshID(bond_->GetMeshID());
      //  Logger::Error("Object ptr passed to anchor was not referencing a mesh!");
    } else {
      bond_lambda_ = lambda * bond_length_;
    }
  } else if (attach_obj_->IsMesh()) {
    mesh_ = dynamic_cast<Mesh *>(attach_obj_);
    mesh_n_bonds_ = mesh_->GetNBonds();
    mesh_length_ = mesh_->GetTrueLength();
    mesh_lambda_ = mesh_length_ * lambda;
    bond_ = mesh_->GetBondAtLambda(mesh_lambda_);
    bond_length_ = bond_->GetLength();
    bond_lambda_ = mesh_lambda_ - bond_->GetMeshLambda();
    SetMeshID(mesh_->GetMeshID());
  } else {
    Logger::Error(
        "Optical traps bound to non-bond objects not yet implemented in "
        "AttachObjLambda.");
  }

  if (bond_lambda_ < 0 || bond_lambda_ > bond_length_) {
    printf("bond_lambda: %2.2f\n", bond_lambda_);
    Logger::Error("Lambda passed to anchor does not match length of "
                  "corresponding bond! lambda: %2.2f, bond_length: %2.2f ",
                  bond_lambda_, bond_length_);
  }
}
void OpticalTrap::Draw(std::vector<graph_struct *> &graph_array) {
  Object::Draw(graph_array);
  double bead_pos[3] = {};
  double const *const bond_position = bond_->GetPosition();
  double const *const bond_orientation = bond_->GetOrientation();
  for (int i = 0; i < n_dim_; ++i) {
    bead_pos[i] = bond_position[i] +
                  (bond_lambda_ - 0.5 * bond_length_) * bond_orientation[i];
  }
  bead_.SetPosition(bead_pos);
  bead_.UpdatePeriodic();
  bead_.Draw(graph_array);
}
