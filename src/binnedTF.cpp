#include <string>
#include <iostream>
#include <memory>

#include "binnedTF.h"
#include "utils.h"

#include "TH2.h"
#include "TFile.h"

BinnedTF::BinnedTF(const std::string particleName, const std::string histName, std::unique_ptr<TFile>& file) : 
  _particleName(particleName),
  _TF( static_cast<TH2D*>( file->Get(histName.c_str()) ) )
  {
  if(!_TF){
    std::cerr << "Error when defining binned TF for particle " << particleName << ": unable to retrieve " << histName << " from file " << file->GetPath() << ".\n";
    exit(1);
  }
  _TF->SetDirectory(0);

  std::cout << "Creating TF component for " << particleName << " from histogram " << histName << ".\n";

  _deltaMin = _TF->GetYaxis()->GetXmin();
  _deltaMax = _TF->GetYaxis()->GetXmax();
  _deltaRange = _deltaMax - _deltaMin; 
  _EgenMax = _TF->GetXaxis()->GetXmax();
  _EgenMin = _TF->GetXaxis()->GetXmin();
 
  std::cout << "Delta range is " << _deltaRange << ", min. and max. values are " << _EgenMin << ", " <<_EgenMax << std::endl << std::endl;
}
