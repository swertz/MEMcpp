#include "Math/Vector4D.h"

#include "MEEvent.h"

void MEEvent::SetVectors(const ROOT::Math::PtEtaPhiEVector &ep, const ROOT::Math::PtEtaPhiEVector &mum, const ROOT::Math::PtEtaPhiEVector &b, const ROOT::Math::PtEtaPhiEVector &bbar, const ROOT::Math::PtEtaPhiEVector &met){
  p3 = ep;
  p5 = mum;
  p4 = b;
  p6 = bbar;
  Met = met;
}
