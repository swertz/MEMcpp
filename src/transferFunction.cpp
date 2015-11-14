#include <string>
#include <iostream>
#include <map>

#include "TFile.h"

#include "transferFunction.h"
#include "binnedTF.h"

TransferFunction::TransferFunction(const std::string file):
  _file( new TFile(file.c_str(), "READ") )  
  {
  if( _file->IsZombie() ){
    std::cerr << "Error opening TF file " << file << ".\n";
    exit(1);
  }
}

void TransferFunction::DefineComponent(const std::string particleName, const std::string histName){
  if( _TF.find(particleName) != _TF.end() ){
    std::cerr << "Error: TF component for " << particleName << " is already defined!" << std::endl;
    exit(1);
  }

  std::cout << "Adding TF component for " << particleName << " from histogram " << histName << ".\n";

  _TF[particleName] = std::unique_ptr<BinnedTF>( new BinnedTF(particleName, histName, _file) );
}

void TransferFunction::Close(){
  _file->Close();
}
