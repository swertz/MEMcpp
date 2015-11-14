#include "Math/Vector4D.h"

#include "MEEvent.h"

void MEEvent::SetVectors(const ROOT::Math::PtEtaPhiEVector &lep1, const ROOT::Math::PtEtaPhiEVector &lep2, const ROOT::Math::PtEtaPhiEVector &bjet1, const ROOT::Math::PtEtaPhiEVector &bjet2, const ROOT::Math::PtEtaPhiEVector &met){
  _p3 = lep1;
  _p5 = lep2;
  _p4 = bjet1;
  _p6 = bjet2;
  _Met = met;
}

void MEEvent::SetLepType(LepType lepType){
  _lepType = lepType;
}
