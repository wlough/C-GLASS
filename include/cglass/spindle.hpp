#ifndef _CGLASS_SPINDLE_H_
#define _CGLASS_SPINDLE_H_

#include "br_rod.hpp"
#include "filament.hpp"
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

class Spindle : public BrRod {
protected:
  spindle_parameters *sparams_;
  filament_parameters *fparams_;
  std::vector<Filament> filaments_;
  std::vector<Site> nuc_sites_;
  std::vector<Site *> fil_sites_;
  std::vector<double> theta_; // reference coordinates for nucleation site
  std::vector<double> phi_;   // positions in body frame
  bool midstep_;
  bool alignment_potential_;
  bool fixed_spacing_;
  int n_filaments_bud_;
  int n_filaments_mother_;
  int n_filaments_;
  double k_spring_;
  double k_align_;
  double spring_length_;
  double anchor_distance_;
  double spb_diameter_;
  void ApplyForcesTorques();
  void ApplyNucleationSiteForces();
  void ApplySpindleForces();
  void ApplyBoundaryForces();
  void InsertSpindle();
  void Integrate();
  void ResetSitePositions();
  void SetParameters();

public:
  Spindle(unsigned long seed);
  void Init(spindle_parameters *sparams);
  void GenerateNucleationSites();
  void InsertFilament(int i_fil);
  void UpdatePosition() {}
  void UpdatePosition(bool midstep);
  virtual void GetInteractors(std::vector<Object *> &ix);
  virtual int GetCount();
  virtual void Draw(std::vector<graph_struct *> &graph_array);
  virtual void ZeroForce();
  void InitFilamentParameters(filament_parameters *fparams);
  void UpdateDrTot();
  void ZeroDrTot();
  const double GetDrTot();
  const int GetNFilaments();
  const bool CheckInteractorUpdate();
  void WriteCheckpoint(std::fstream &ocheck);
  void WriteSpec(std::fstream &ospec);
  void ReadCheckpoint(std::fstream &icheck);
  void ReadSpec(std::fstream &ispec);

  // Convert binary data to text. Static to avoid needing to istantiate
  // species members.
  static void ConvertSpec(std::fstream &ispec, std::fstream &otext);
  static void WriteSpecTextHeader(std::fstream &otext);
};

#endif // _CGLASS_SPINDLE_H_
