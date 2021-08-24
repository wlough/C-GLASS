#ifndef _CGLASS_RECEPTOR_H_
#define _CGLASS_RECEPTOR_H_

#include "site.hpp"
#include "cortex.hpp"
#include "species_base.hpp"

/* receptor.hpp declares the Receptor class. The Receptor class represents
 * biological receptors that can attach to crosslinking and motor proteins. 
 * In this code, Receptor's are modelled as sites on a mesh that can attach 
 * to the Anchor of a Crosslinker.
*/

class Receptor: public Sphere {
private:
  receptor_parameters *sparams_; // Input parameters unique to receptors
  SpeciesBase* pc_species_ = nullptr; // The species that the PointCover corresponds to
  Object* pc_object_ = nullptr; // The object that the receptor is on
  double s_; // The position (length) along object of the receptor (for species-PointCover)
  int i_; // The index of the object that the receptor is on (for species-PointCover)
  //Receptor* prevN_ = nullptr;
  //Receptor* nextN_ = nullptr;
public:
  Receptor(unsigned long seed); // Constructor with RNG seed input

  // Save yaml input file parameters locally
  void Init(receptor_parameters *sparams);

  // Initialize Receptor members with yaml input file parameters
  void SetParameters();

  // Initialize positions along PointCover objects
  void SetLocations(int i, double s);
  //void SetNeighbors(Receptor* prevN, Receptor* nextN);
  // Overloaded UpdatePositions to follow PointCover objects.
  void UpdatePosition();

  //Receptor* GetPlusNeighbor();
  int GetMinusNeighbor();

  // Setters
  void SetPCSpecies(SpeciesBase* pc_species);
  void SetPCObject(Object* pc_object);

  // Add/Sub forces/torques to objects receptors are on
  void CalcTorque();
  void AddForce(const double *const force);
  void AddTorque(const double *const torque);
  void SubForce(const double *const force);
  void SubTorque(const double *const torque);

  // Read/write binaries
  void WriteSpec(std::fstream &ospec);
  void ReadSpec(std::fstream &ispec);

  // Convert binary data to text. Static to avoid needing to istantiate
  // species members.
  static void ConvertSpec(std::fstream &ispec, std::fstream &otext);
  static void WriteSpecTextHeader(std::fstream &otext);
};

#endif
