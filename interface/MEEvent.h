#ifndef _INC_MEEVENT
#define _INC_MEEVENT

#include "Math/Vector4D.h"
#include "TopEffTh/TopEffTh.h"

class MEEvent{
  public:

  void SetVectors(const ROOT::Math::PtEtaPhiEVector &lep1, const ROOT::Math::PtEtaPhiEVector &lep2, const ROOT::Math::PtEtaPhiEVector &bjet1, const ROOT::Math::PtEtaPhiEVector &bjet2, const ROOT::Math::PtEtaPhiEVector &met);
  void SetLepType(LepType lepType);

  inline const ROOT::Math::PtEtaPhiEVector& GetP3() const { return _p3; }
  inline const ROOT::Math::PtEtaPhiEVector& GetP4() const { return _p4; }
  inline const ROOT::Math::PtEtaPhiEVector& GetP5() const { return _p5; }
  inline const ROOT::Math::PtEtaPhiEVector& GetP6() const { return _p6; }
  inline const ROOT::Math::PtEtaPhiEVector& GetMet() const { return _Met; }

  inline const LepType GetLepType() const { return _lepType; }

  private:

  ROOT::Math::PtEtaPhiEVector _p3, _p4, _p5, _p6, _Met;
  LepType _lepType;
};

#endif
