#include "cglass/analysis.hpp"
#include "cglass/space_base.hpp"
void AnalysisBase::SetParams(system_parameters *params) { params_ = params; }
void AnalysisBase::SetSpace(SpaceBase *space) { space_ = space; }
const SpaceBase *AnalysisBase::space_ = nullptr;
const system_parameters *AnalysisBase::params_ = nullptr;
