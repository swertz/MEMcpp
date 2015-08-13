#ifndef _INC_MEEVENT
#define _INC_MEEVENT

#include "Math/Vector4D.h"

class MEEvent{
  public:

  void SetVectors(const ROOT::Math::PtEtaPhiEVector &ep, const ROOT::Math::PtEtaPhiEVector &mum, const ROOT::Math::PtEtaPhiEVector &b, const ROOT::Math::PtEtaPhiEVector &bbar, const ROOT::Math::PtEtaPhiEVector &met);

  inline const ROOT::Math::PtEtaPhiEVector& GetP3() const { return p3; }
  inline const ROOT::Math::PtEtaPhiEVector& GetP4() const { return p4; }
  inline const ROOT::Math::PtEtaPhiEVector& GetP5() const { return p5; }
  inline const ROOT::Math::PtEtaPhiEVector& GetP6() const { return p6; }
  inline const ROOT::Math::PtEtaPhiEVector& GetMet() const { return Met; }

  private:

  ROOT::Math::PtEtaPhiEVector p3, p4, p5, p6, Met;
};

#endif
