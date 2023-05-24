#ifndef _CGLASS_RECEPTOR_H_
#define _CGLASS_RECEPTOR_H_

#include "cortex.hpp"
#include "site.hpp"
#include "species_base.hpp"

/* receptor.hpp declares the Receptor class. The Receptor class represents
 * biological receptors that can attach to crosslinking and motor proteins. 
 * In this code, Receptor's are modelled as sites on a mesh that can attach 
 * to the Anchor of a Crosslinker.
*/

class Receptor : public Object {
protected:
  bool induces_catastrophe_ = false; //  upon double binding
private:
  receptor_parameters *sparams_;      // Input parameters unique to receptors
  SpeciesBase *pc_species_ = nullptr; // The species of PointCover
  Receptor *prev_r_ = nullptr; // The previous receptor (in - direction) on fila
  Receptor *next_r_ = nullptr; // The next receptor (in + direction) on fila
  Object *pc_object_ = nullptr; // point cover object sphere is on
  double s_ = 0;                //how far the sphere is along a filament
  double step_size_ = 0;
  int i_; // The index of the object that the receptor is on (for species-PointCover)
public:
  Receptor(unsigned long seed); // Constructor with RNG seed input
  void Init(receptor_parameters *sparams);
  void SetParameters();

  bool InducesCatastrophe() { return induces_catastrophe_; }

  void SetSphereS(double s);
  void SetLocations(int i, double s);
  void SetStepSize(double step_size);
  void SetPCSpecies(SpeciesBase *pc_species);
  void SetPCObject(Object *pc_object);
  void SetPCObjectForSphere(Object *pc_object);
  void SetNeighbors(Receptor *prevN, Receptor *nextN);

  double GetSphereS();
  double GetStepSize();
  Object *GetPCObject();
  Object *GetPCObjectForSphere();
  Receptor *GetPlusNeighbor();
  Receptor *GetMinusNeighbor();

  // Overloaded UpdatePositions to follow PointCover objects.
  void UpdatePosition();

  // Add/Sub forces/torques to objects receptors are on
  void CalcTorque();
  void AddForce(const double *const force);
  void AddTorque(const double *const torque);
  void SubForce(const double *const force);
  void SubTorque(const double *const torque);

  // Read/write binaries
  void ReadSpec(std::fstream &ispec);
  void WriteSpec(std::fstream &ospec);

  // Convert binary data to text. Static to avoid needing to istantiate
  // species members.
  static void ConvertSpec(std::fstream &ispec, std::fstream &otext);
  static void WriteSpecTextHeader(std::fstream &otext);
};

#endif
