#include "Math/Vector4D.h"

#include "MEEvent.h"

void MEEvent::SetVectors(const ROOT::Math::PtEtaPhiEVector &ep, const ROOT::Math::PtEtaPhiEVector &mum, const ROOT::Math::PtEtaPhiEVector &b, const ROOT::Math::PtEtaPhiEVector &bbar, const ROOT::Math::PtEtaPhiEVector &met){
  _p3 = ep;
  _p5 = mum;
  _p4 = b;
  _p6 = bbar;
  _Met = met;
}
