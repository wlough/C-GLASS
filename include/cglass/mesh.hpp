#ifndef _CGLASS_MESH_H_
#define _CGLASS_MESH_H_

#include "bond.hpp"
#include "composite.hpp"

typedef std::vector<Bond>::iterator bond_iterator;
typedef std::vector<Site>::iterator site_iterator;

class Mesh : public Composite {
 private:
  static int _next_mesh_id_;
  static std::mutex _mesh_mtx_;
  void InitMeshID();

 protected:
  bool anchored_ = false;
  bool midstep_ = true;
  bool posits_only_ = false;
  bool dynamic_instability_flag_ = false;
  bool polydispersity_flag_ = false;
  int n_sites_ = 0;
  int n_bonds_ = 0;
  int n_bonds_max_ = 0;
  std::vector<Site> sites_;
  std::vector<Site*> site_ptrs_; // for meshes w/ sites that are species
  std::vector<Bond> bonds_;
  double bond_length_ = -1;
  double true_length_ = -1;
  Bond *GetRandomBond();
  virtual void UpdateInteractors();
  void UpdateSiteOrientations();
  void RelocateMesh(double const *const new_pos, double const *const u);
  
 public:
  Mesh(unsigned long seed);
  void AddRandomBondToSite(double l, int i_site);
  void AddBondToSite(double *u, double l, int i_site);
  void AddSite(Site s);
  void AddSitePtr(Site *s);
  void AddBond(Site *s1, Site *s2);
  void InitSiteAt(double *new_pos, double d);
  void InitBondAt(double *new_pos, double *u, double l, double d);
  void InitRandomSite(double d);
  void InitRandomBond(double d);
  void InitRandomBondOriented(double *u, double d);
  void AddRandomBondAnywhere(double l, double d = 1);
  void AddRandomBondToTip(double l);
  void AddBondToTip(double *u, double l);
  void SetBondLength(double l);
  void RemoveBondFromTip();
  void ReportSites();
  void ReportBonds();
  void Report();
  void SubReport();
  void UpdateBondPositions();
  void UpdatePrevPositions();
  virtual void ResetPreviousPosition();
  virtual void Draw(std::vector<graph_struct *> &graph_array);
  virtual void Reserve();
  void Clear();
  void DoubleGranularityLinear();
  void HalfGranularityLinear();
  int GetNBonds() { return n_bonds_; }
  Bond *GetBondAtLambda(double lambda);
  Site *GetSite(int i);
  Bond *GetBond(int i);
  virtual void ZeroForce();
  virtual void GetInteractors(std::vector<Object *> &ix);
  virtual int GetCount();
  virtual void ReadPosit(std::fstream &ip);
  virtual void WritePosit(std::fstream &op);
  virtual void ReadSpec(std::fstream &ip);
  virtual void WriteSpec(std::fstream &op);
  virtual void ReadCheckpoint(std::fstream &ip);
  virtual void WriteCheckpoint(std::fstream &op);
  virtual void ScalePosition();
  virtual void UpdateDrTot();
  virtual double const GetDrTot();
  virtual void ZeroDrTot();
  virtual void SetNBondsMax(int n) {n_bonds_max_ = n;}
  virtual void SetPosition(double const *const new_pos);
  //virtual std::vector<Interaction *> *GetInteractions();
  virtual void ClearInteractions();
  virtual void GetAvgPosition(double *ap);
  virtual void GetAvgScaledPosition(double *asp);
  virtual void GetAvgOrientation(double *au);
  virtual void SetAvgPosition();
  virtual void GetContactNumbers(std::vector<double> &cn);
  virtual void CalcPolarOrder();
  virtual void GetPolarOrders(std::vector<double> &po);
  std::pair<double, double> GetAvgOrientationCorrelation();
  virtual void ZeroOrientationCorrelations();
  virtual const double GetBondLength() const;
  virtual const bool CheckInteractorUpdate();
  virtual const double GetLambdaAtBond(int bond_oid);
  virtual const double GetTrueLength() const;

  // Convert binary data to text. Static to avoid needing to istantiate
  // species members.
  static void ConvertSpec(std::fstream &ispec, std::fstream &otext);
  static void WriteSpecTextHeader(std::fstream &otext);
};

#endif  // _CGLASS_MESH_H_
