#ifndef _INC_TRANSFERFUNCTION
#define _INC_TRANSFERFUNCTION

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>

#include "TFile.h"

#include "binnedTF.h"

class TransferFunction{
  public:

  TransferFunction(const std::string file);
  ~TransferFunction() {};

  void DefineComponent(const std::string particleName, const std::string histName);

  inline double Evaluate(const std::string &particleName, const double &Erec, const double &Egen);
  inline double GetDeltaRange(const std::string &particleName, const double &Erec);
  inline double GetDeltaMin(const std::string &particleName, const double &Erec);
  inline double GetDeltaMax(const std::string &particleName, const double &Erec);
  void Close();

  private:

  std::unique_ptr<TFile> _file;

  std::map< std::string, std::unique_ptr<BinnedTF> > _TF;
};

inline double TransferFunction::Evaluate(const std::string &particleName, const double &Erec, const double &Egen){
  // We might want to avoid to check this every time, since it might slow things down too much
  /*if( _TF.find(particleName) == _TF.end() ){
    std::cerr << "Error: TF component for " << particleName << " is not defined!" << std::endl;
    exit(1);
  }*/

  return _TF.at(particleName)->Evaluate(Erec, Egen);
}

inline double TransferFunction::GetDeltaRange(const std::string &particleName, const double &Erec){
  // We might want to avoid to check this every time, since it might slow things down too much
  /*if( _TF.find(particleName) == _TF.end() ){
    std::cerr << "Error: TF component for " << particleName << " is not defined!" << std::endl;
    exit(1);
  }*/

  return _TF.at(particleName)->GetDeltaRange(Erec);
}

inline double TransferFunction::GetDeltaMin(const std::string &particleName, const double &Erec){
  // We might want to avoid to check this every time, since it might slow things down too much
  /*if( _TF.find(particleName) == _TF.end() ){
    std::cerr << "Error: TF component for " << particleName << " is not defined!" << std::endl;
    exit(1);
  }*/

  return _TF.at(particleName)->GetDeltaMin(Erec);
}

inline double TransferFunction::GetDeltaMax(const std::string &particleName, const double &Erec){
  // We might want to avoid to check this every time, since it might slow things down too much
  /*if( _TF.find(particleName) == _TF.end() ){
    std::cerr << "Error: TF component for " << particleName << " is not defined!" << std::endl;
    exit(1);
  }*/

  return _TF.at(particleName)->GetDeltaMax(Erec);
}

#endif
