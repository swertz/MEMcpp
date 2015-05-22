#include <string>
#include <iostream>
#include <stdlib.h>

#include "binnedTF.h"
#include "utils.h"

#include "TH2.h"
#include "TFile.h"

using namespace std;

BinnedTF::BinnedTF(const std::string particleName, const std::string histName, TFile* file){
  _TF = dynamic_cast<TH2D*>( file->Get(histName.c_str()) );
  if(!_TF){
    cerr << "Error when defining binned TF for paricle " << particleName << ": unable to retrieve " << histName << " from file " << file->GetPath() << ".\n";
    exit(1);
  }

  cout << "Creating TF component for " << particleName << " from histogram " << histName << ".\n";

  _particleName = particleName;
  _histDeltaMin = _TF->GetYaxis()->GetXmin();
  _histDeltaMax = _TF->GetYaxis()->GetXmax();
  //_histDeltaMin = 0.; 
  //_histDeltaMax = 0.;
  _deltaMin = _histDeltaMin;
  _deltaMax = _histDeltaMax;
  _deltaRange = _deltaMax - _deltaMin; 
  _EgenMax = _TF->GetXaxis()->GetXmax();
 
  cout << "Delta range is " << _deltaRange << ", max. energy is " << _EgenMax << endl << endl;
}

BinnedTF::~BinnedTF(){
  delete _TF; _TF = NULL;
}

double BinnedTF::Evaluate(const double Erec, const double Egen) const {
  //cout << "Evaluating TF for particle " << _particleName << ": Erec = " << Erec << ", Egen = " << Egen;
  
  double delta = Erec - Egen;
  
  if(Egen < 0. || Erec < 0. || Egen > _EgenMax || delta > _deltaMax || delta < _deltaMin){
    //cout << "... Out of range!" << endl;
    return 0.;
  }

  const int xBin = _TF->GetXaxis()->FindBin(Egen);
  const int yBin = _TF->GetYaxis()->FindBin(delta);

  //cout << ", TF = " << _TF->GetBinContent(xBin, yBin) << endl;
  
  return _TF->GetBinContent(xBin, yBin);
}

double BinnedTF::GetDeltaRange(const double Erec) const {
  return min(_deltaMax, Erec) - _deltaMin;
}

double BinnedTF::GetDeltaMin() const {
  return _deltaMin;
}

double BinnedTF::GetDeltaMax(const double Erec) const {
  return min(_deltaMax, Erec);
}

/*void BinnedTF::SetDeltaRange(double min, double max){
  if(min < _histDeltaMin)
    min = _histDeltaMin;
  
  if(max > _histDeltaMax)
    max = _histDeltaMax;

  _deltaMin = min;
  _deltaMax = max;
  _deltaRange = _deltaMax - _deltaMin;
}*/
