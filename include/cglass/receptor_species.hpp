#ifndef _CGLASS_RECEPTOR_SPECIES_H_
#define _CGLASS_RECEPTOR_SPECIES_H_

#include "receptor.hpp"
#include "species.hpp"

/* receptor_species.hpp declares the ReceptorSpecies Class to represent 
 * ReceptorSpecies as a special type of species that can set the mesh that 
 * the receptors are on and control number of receptors placed on the mesh 
 * via the concentration option.
*/

class ReceptorSpecies: public Species<Receptor, species_id::receptor> {

private:
  double concentration_; // concentration (# receptors/(object surface area))
  std::string component_; // concentration (# receptors/(object surface area))

  PointCover* pc_ = nullptr; // ptr to the PointCover object the receptor is on
  SpeciesBase* pc_species_ = nullptr; // The species that the PointCover corresponds to
  double s_; // The position (length) along object of the receptor (for species-PointCover)
  int i_; // The index of the object that the receptor is on (for species-PointCover)

public:
  ReceptorSpecies(unsigned long seed); // Constructor with RNG seed input

  // Initialize with parameters from the input yaml file
  void Init(std::string spec_name, ParamsParser &parser);

  // Set the mesh that the receptor belongs to (currently only works with
  // Cortex pointcover)
  void SetPC(Cortex* cx, std::vector<SpeciesBase *> &species);

  // Reserve space for the vector of receptors
  void Reserve();

  // Create a receptor object and save it as a site on mesh_
  void AddMember();
};

#endif
