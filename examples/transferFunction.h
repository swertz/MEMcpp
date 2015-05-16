#ifndef _INC_TRANSFERFUNCTION
#define _INC_TRANSFERFUNCTION

#include <map>
#include <string>

#include "TFile.h"

#include "binnedTF.h"

class TransferFunction{
  public:

  TransferFunction(const std::string file);
  ~TransferFunction();

  void DefineComponent(const std::string particleName, const std::string histName);

  double Evaluate(const std::string particleName, const double Erec, const double Egen);
  double GetDeltaRange(const std::string particleName);
  double GetDeltaMin(const std::string particleName);
  double GetDeltaMax(const std::string particleName);
  void SetDeltaRange(const std::string particleName, const double min, const double max);

  private:

  TFile* _file;

  std::map< std::string, BinnedTF* > _TF;
};

#endif
