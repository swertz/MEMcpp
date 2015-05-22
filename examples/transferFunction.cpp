#include <string>
#include <iostream>
#include <map>
#include <stdlib.h>

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
    delete i.second; i.second = NULL;
  }
  delete _file; _file = NULL;
}
  
void TransferFunction::DefineComponent(const std::string particleName, const std::string histName){
  if( _TF.find(particleName) != _TF.end() ){
    cerr << "Error: TF component for " << particleName << " is already defined!" << endl;
    exit(1);
  }

  cout << "Adding TF component for " << particleName << " from histogram " << histName << ".\n";

  _TF[particleName] = new BinnedTF(particleName, histName, _file);
}

double TransferFunction::Evaluate(const std::string particleName, const double Erec, const double Egen){
  // We might want to avoid to check this every time, since it might slow things down too much
  if( _TF.find(particleName) == _TF.end() ){
    cerr << "Error: TF component for " << particleName << " is not defined!" << endl;
    exit(1);
  }

  return _TF[particleName]->Evaluate(Erec, Egen);
}

double TransferFunction::GetDeltaRange(const std::string particleName, const double Erec){
  // We might want to avoid to check this every time, since it might slow things down too much
  if( _TF.find(particleName) == _TF.end() ){
    cerr << "Error: TF component for " << particleName << " is not defined!" << endl;
    exit(1);
  }

  return _TF[particleName]->GetDeltaRange(Erec);
}

double TransferFunction::GetDeltaMin(const std::string particleName){
  // We might want to avoid to check this every time, since it might slow things down too much
  if( _TF.find(particleName) == _TF.end() ){
    cerr << "Error: TF component for " << particleName << " is not defined!" << endl;
    exit(1);
  }

  return _TF[particleName]->GetDeltaMin();
}

double TransferFunction::GetDeltaMax(const std::string particleName, const double Erec){
  // We might want to avoid to check this every time, since it might slow things down too much
  if( _TF.find(particleName) == _TF.end() ){
    cerr << "Error: TF component for " << particleName << " is not defined!" << endl;
    exit(1);
  }

  return _TF[particleName]->GetDeltaMax(Erec);
}

/*void TransferFunction::SetDeltaRange(const std::string particleName, double min, double max){
  // We might want to avoid to check this every time, since it might slow things down too much
  if( _TF.find(particleName) == _TF.end() ){
    cerr << "Error: TF component for " << particleName << " is not defined!" << endl;
    exit(1);
  }

  return _TF[particleName]->SetDeltaRange(min, max);
}*/
