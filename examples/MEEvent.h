#ifndef _INC_MEEVENT
#define _INC_MEEVENT

#include "TLorentzVector.h"
#include "TH1D.h"

class MEEvent{
  public:

  MEEvent();
  void SetVectors(const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met);
  ~MEEvent();

  inline const TLorentzVector& GetP3() const { return p3; }
  inline const TLorentzVector& GetP4() const { return p4; }
  inline const TLorentzVector& GetP5() const { return p5; }
  inline const TLorentzVector& GetP6() const { return p6; }
  inline const TLorentzVector& GetMet() const { return Met; }

  void writeHists(void);

  TH1D* GetTTbar();
  
  private:

  TLorentzVector p3, p4, p5, p6, Met;

  TH1D *hst_TTbar;
};

#endif
