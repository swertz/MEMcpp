#include <string>
#include <iostream>
#include <map>

#include "TFile.h"

#include "transferFunction.h"
#include "binnedTF.h"

using namespace std;

TransferFunction::TransferFunction(const std::string file){
  _file = new TFile(file.c_str(), "READ");
  
  if( _file->IsZombie() ){
    cerr << "Error opening TF file " << file << ".\n";
    exit(1);
  }
}

TransferFunction::~TransferFunction(){
  for(auto &i: _TF){
    delete i.second; i.second = nullptr;
  }
  delete _file; _file = nullptr;
}
  
void TransferFunction::DefineComponent(const std::string particleName, const std::string histName){
  if( _TF.find(particleName) != _TF.end() ){
    cerr << "Error: TF component for " << particleName << " is already defined!" << endl;
    exit(1);
  }

  cout << "Adding TF component for " << particleName << " from histogram " << histName << ".\n";

  _TF[particleName] = new BinnedTF(particleName, histName, _file);
}

