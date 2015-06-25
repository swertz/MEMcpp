#include <string>
#include <iostream>

#include "binnedTF.h"
#include "utils.h"

#include "TH2.h"
#include "TFile.h"

BinnedTF::BinnedTF(const std::string file, const std::string particleName, const std::string histName) : _particleName(particleName) {

  _file = new TFile(file.c_str());
  if(_file->IsZombie()){
    std::cerr << "Error opening transfer function file " << file << ".\n";
    exit(1);
  }

  _TF = dynamic_cast<TH2D*>( _file->Get(histName.c_str()) );
  if(!_TF){
    std::cerr << "Error when defining binned TF for paricle " << particleName << ": unable to retrieve " << histName << " from file " << _file->GetPath() << ".\n";
    exit(1);
  }

  std::cout << "Creating TF component for " << particleName << " from histogram " << histName << ".\n";

  _deltaMin = _TF->GetYaxis()->GetXmin();
  _deltaMax = _TF->GetYaxis()->GetXmax();
  _deltaRange = _deltaMax - _deltaMin; 
  _EgenMax = _TF->GetXaxis()->GetXmax();
  _EgenMin = _TF->GetXaxis()->GetXmin();
 
  std::cout << "Delta range is " << _deltaRange << ", min. and max. values are " << _EgenMin << ", " <<_EgenMax << std::endl << std::endl;
}

BinnedTF::~BinnedTF(){
  if(_file->IsOpen()){
    _file->Close();
    _file = nullptr;
  }
  delete _TF; _TF = nullptr;
}

