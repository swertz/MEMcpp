#ifndef _INC_MEEVENT
#define _INC_MEEVENT

#include "Math/Vector4D.h"

class MEEvent{
  public:

  void SetVectors(const ROOT::Math::PtEtaPhiEVector &ep, const ROOT::Math::PtEtaPhiEVector &mum, const ROOT::Math::PtEtaPhiEVector &b, const ROOT::Math::PtEtaPhiEVector &bbar, const ROOT::Math::PtEtaPhiEVector &met);

  inline const ROOT::Math::PtEtaPhiEVector& GetP3() const { return _p3; }
  inline const ROOT::Math::PtEtaPhiEVector& GetP4() const { return _p4; }
  inline const ROOT::Math::PtEtaPhiEVector& GetP5() const { return _p5; }
  inline const ROOT::Math::PtEtaPhiEVector& GetP6() const { return _p6; }
  inline const ROOT::Math::PtEtaPhiEVector& GetMet() const { return _Met; }

  private:

  ROOT::Math::PtEtaPhiEVector _p3, _p4, _p5, _p6, _Met;
};

#endif
