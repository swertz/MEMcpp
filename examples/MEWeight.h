#ifndef _INC_MEWEIGHT
#define _INC_MEWEIGHT

#include <string>
#include <vector>

#include "SubProcesses/P0_Sigma_sm_gg_epvebmumvmxbx/CPPProcess.h"

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"

#include "TLorentzVector.h"

#include "transferFunction.h"
#include "MEEvent.h"

int CUBAIntegrand(const int *nDim, const double* Xarg, const int *nComp, double *Value, void *inputs, const int *nVec, const int *core, const double *weight);

class MEWeight{
  public:

  double Integrand(const double* Xarg, const double *weight);
  double ComputePdf(const int pid, const double x, const double q2);
  inline void setProcessMomenta(vector<double*> &p){ process.setMomenta(p); }
  inline void computeMatrixElements(){ process.sigmaKin(); }
  inline const double* const getMatrixElements() const { return process.getMatrixElements(); }
  double ComputeWeight(double &error);
  MEEvent* GetEvent();
  void SetEvent(const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met);
  void AddTF(const std::string particleName, const std::string histName);

  //void WriteHist();

  MEWeight(const std::string paramCardPath, const std::string pdfName, const std::string fileTF);
  ~MEWeight();

  private:

  //TH1D* hst_TTbar;

  CPPProcess process;
  LHAPDF::PDF* pdf;
  MEEvent* myEvent;
  TransferFunction* myTF;
};

#endif
