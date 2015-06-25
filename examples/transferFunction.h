#ifndef _INC_TRANSFERFUNCTION
#define _INC_TRANSFERFUNCTION

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "binnedTF.h"

template <typename tfType>
class TransferFunction{
  public:

  //TransferFunction(const std::string file);
  ~TransferFunction();

  void DefineComponent(const std::string fileName, const std::string particleName, const std::string histName);

  inline double Evaluate(const std::string particleName, const double Erec, const double Egen);
  inline double GetDeltaRange(const std::string particleName, const double Erec);
  inline double GetDeltaMin(const std::string particleName, const double Erec);
  inline double GetDeltaMax(const std::string particleName, const double Erec);

  private:

  std::map< std::string, tfType > _TF;
};

template <typename tfType>
inline double TransferFunction<tfType>::Evaluate(const std::string particleName, const double Erec, const double Egen){
  // We might want to avoid to check this every time, since it might slow things down too much
  /*if( _TF.find(particleName) == _TF.end() ){
    std::cerr << "Error: TF component for " << particleName << " is not defined!" << std::endl;
    exit(1);
  }*/

  return _TF[particleName]->Evaluate(Erec, Egen);
}

template <typename tfType>
inline double TransferFunction<tfType>::GetDeltaRange(const std::string particleName, const double Erec){
  // We might want to avoid to check this every time, since it might slow things down too much
  /*if( _TF.find(particleName) == _TF.end() ){
    std::cerr << "Error: TF component for " << particleName << " is not defined!" << std::endl;
    exit(1);
  }*/

  return _TF[particleName]->GetDeltaRange(Erec);
}

template <typename tfType>
inline double TransferFunction<tfType>::GetDeltaMin(const std::string particleName, const double Erec){
  // We might want to avoid to check this every time, since it might slow things down too much
  /*if( _TF.find(particleName) == _TF.end() ){
    std::cerr << "Error: TF component for " << particleName << " is not defined!" << std::endl;
    exit(1);
  }*/

  return _TF[particleName]->GetDeltaMin(Erec);
}

template <typename tfType>
inline double TransferFunction<tfType>::GetDeltaMax(const std::string particleName, const double Erec){
  // We might want to avoid to check this every time, since it might slow things down too much
  /*if( _TF.find(particleName) == _TF.end() ){
    std::cerr << "Error: TF component for " << particleName << " is not defined!" << std::endl;
    exit(1);
  }*/

  return _TF[particleName]->GetDeltaMax(Erec);
}

#endif
