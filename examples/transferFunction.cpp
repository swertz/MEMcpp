#include <string>
#include <iostream>
#include <map>

//#include "TFile.h"

#include "transferFunction.h"
#include "binnedTF.h"

/*template <typename tfType>
TransferFunction::TransferFunction(const std::string file){
  _file = new TFile(file.c_str(), "READ");
  
  if( _file->IsZombie() ){
    std::cerr << "Error opening TF file " << file << ".\n";
    exit(1);
  }
}*/

template <typename tfType>
TransferFunction<tfType>::~TransferFunction(){
  for(auto &i: _TF){
    delete i.second; i.second = nullptr;
  }
  //delete _file; _file = nullptr;
}
  
template <typename tfType>
void TransferFunction<tfType>::DefineComponent(const std::string fileName, const std::string particleName, const std::string histName){
  if( _TF.find(particleName) != _TF.end() ){
    std::cerr << "Error: TF component for " << particleName << " is already defined!" << std::endl;
    exit(1);
  }

  std::cout << "Adding TF component for " << particleName << " from histogram " << histName << ".\n";

  _TF[particleName] = new tfType(fileName, particleName, histName);
}

