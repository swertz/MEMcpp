#ifndef _INC_MEEVENT
#define _INC_MEEVENT

#include "TLorentzVector.h"

class MEEvent{
  public:

  MEEvent();
  void SetVectors(const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met);
  ~MEEvent();

  inline TLorentzVector GetP3(void) const { return p3; }
  inline TLorentzVector GetP4(void) const { return p4; }
  inline TLorentzVector GetP5(void) const { return p5; }
  inline TLorentzVector GetP6(void) const { return p6; }
  inline TLorentzVector GetMet(void) const { return Met; }

  //void writeHists(void);

  //TH1D* GetTTbar();
  
  private:

  TLorentzVector p3, p4, p5, p6, Met;

  //TH1D *hst_TTbar;
};

#endif
