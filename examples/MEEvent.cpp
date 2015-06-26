#include "TLorentzVector.h"

#include "MEEvent.h"

void MEEvent::SetVectors(const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met){
  p3 = ep;
  p5 = mum;
  p4 = b;
  p6 = bbar;
  Met = met;
}
