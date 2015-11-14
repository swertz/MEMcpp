#ifndef _INC_BINNEDTF
#define _INC_BINNEDTF

#include <string>
#include <algorithm>
#include <memory>

#include "TH2.h"
#include "TFile.h"

class BinnedTF{
  public:

  BinnedTF(const std::string particleName, const std::string histName, std::unique_ptr<TFile>& file);
  ~BinnedTF() {};
  inline double Evaluate(const double &Erec, const double &Egen) const;
  inline double GetDeltaRange(const double &Erec) const;
  inline double GetDeltaMin(const double &Erec) const;
  inline double GetDeltaMax(const double &Erec) const;

  private:

  std::string _particleName;
  double _histDeltaMin, _histDeltaMax;
  double _deltaMin, _deltaMax, _deltaRange;
  double _EgenMax, _EgenMin;
  std::unique_ptr<TH2D> _TF;
};

inline double BinnedTF::Evaluate(const double &Erec, const double &Egen) const {
  //std::cout << "Evaluating TF for particle " << _particleName << ": Erec = " << Erec << ", Egen = " << Egen << std::endl;
  
  double delta = Erec - Egen;
  
  if(Egen < _EgenMin || Egen > _EgenMax || delta > _deltaMax || delta < _deltaMin){
    //std::cout << "... Out of range!" << std::endl;
    return 0.;
  }

  // Use ROOT's global bin number "feature" for 2-dimensional histograms
  const int bin = _TF->FindFixBin(Egen, delta);

  //std::cout << ", TF = " << _TF->GetBinContent(bin) << std::endl;
  return _TF->GetBinContent(bin);
}

inline double BinnedTF::GetDeltaRange(const double &Erec) const {
  return GetDeltaMax(Erec) - GetDeltaMin(Erec);
}

inline double BinnedTF::GetDeltaMin(const double &Erec) const {
  return std::max(_deltaMin, Erec - _EgenMax);
}

inline double BinnedTF::GetDeltaMax(const double &Erec) const {
  return std::min(_deltaMax, Erec - _EgenMin);
}

#endif
