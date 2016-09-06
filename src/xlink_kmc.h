#ifndef _SIMCORE_XLINK_KMC_H_
#define _SIMCORE_XLINK_KMC_H_

#include "auxiliary.h"
#include "kmc_base.h"
#include "lookup_table.h"
#include "xlink_harmonic.h"

class Xlink;
class XlinkHead;

class XlinkKMC : public KMCBase {
  protected:
    double eps_eff_0_1_[2];
    double eps_eff_1_2_[2];
    double on_rate_0_1_[2];
    double on_rate_1_2_[2];
    double f_stall_[2];
    double velocity_[2];
    double velocity_p_scale_[2];
    double velocity_ap_scale_[2];
    double alpha_;
    double barrier_weight_;
    double k_stretch_;
    double mrcut_;
    double mrcut2_;
    double max_length_;
    double rcutoff_0_1_;
    double rcutoff_1_2_;
    double r_equil_;
    double polar_affinity_;

    bool end_pause_[2];
    bool first_potential_use = false;

    int nfree_;
    int stall_type_;
    int nbound1_[2];
    int nbound2_[2];

    int nsimples_;
    bool write_event_ = false;

    std::ostringstream kmc_file_name_;
    std::ofstream kmc_file;

    std::vector<Simple*>* simples_;
    std::unordered_map<int, int>* oid_position_map_;
    nl_list *neighbors_;
    LookupTable n_exp_lookup_;

    void KMC_0_1();
    void KMC_1_0();
    void KMC_1_2();
    void KMC_2_1();
    void KMC_2_1_ForceIndep();
    void KMC_2_1_ForceDep();
    void Update_0_1(Xlink *xit);
    void Update_1_2(Xlink *xit);
    void UpdateStage1(Xlink *xit);
    void UpdateStage2(Xlink *xit);
    void ApplyStage2Force(Xlink *xit);

    // Helper functions to detach crosslinks
    void Detach_1_0(Xlink *xit, XlinkHead *freehead, XlinkHead *boundhead);
    void Detach_2_1(Xlink *xit, int headtype);
    void Detach_2_0(Xlink *xit);

    // Specifics to 2 stage crosslinks (building tables, numbers, etc)
    double XKMCErfinv(double x);
    void BuildTables();

    void CalcCutoff();

  public:
    virtual void Init(space_struct *pSpace, ParticleTracking *pTracking,
        PotentialManager *pPotentials,
        SpeciesBase *spec1, SpeciesBase *spec2, int ikmc, YAML::Node &node,
        long seed);
    virtual void Print();
    virtual void Dump();
    virtual double GetMaxRcut() {
      return rcutoff_1_2_;
    }

    virtual void PrepKMC();
    virtual void StepKMC();
    virtual void UpdateKMC();
    virtual void TransferForces();

    virtual void PrepOutputs();
    virtual void WriteOutputs(int istep);
    void WriteEvent(const std::string &pString);

};

#endif
