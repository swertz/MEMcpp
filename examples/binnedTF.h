#ifndef _INC_BINNEDTF
#define _INC_BINNEDTF

#include <string>

#include "TH2.h"
#include "TFile.h"

class BinnedTF{
  public:

  BinnedTF(const std::string particleName, const std::string histName, TFile* file);
  ~BinnedTF();
  double Evaluate(const double Erec, const double Egen) const;
  double GetDeltaRange() const;
  double GetDeltaMin() const;
  double GetDeltaMax() const;
  void SetDeltaRange(double min, double max);

  private:

  std::string _particleName;
  double _histDeltaMin, _histDeltaMax;
  double _deltaMin, _deltaMax, _deltaRange;
  double _EgenMax;
  const TH2D* _TF;
};

#endif
