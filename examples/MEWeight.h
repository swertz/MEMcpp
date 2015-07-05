#ifndef _INC_MEWEIGHT
#define _INC_MEWEIGHT

#include <string>
#include <vector>

#include "SubProcesses/P0_Sigma_sm_gg_epvebmumvmxbx/CPPProcess.h"

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"

#include "Math/Vector4D.h"
#include "TH1D.h"

#include "transferFunction.h"
#include "MEEvent.h"

int CUBAIntegrand(const int *nDim, const double* Xarg, const int *nComp, double *Value, void *inputs, const int *nVec, const int *core, const double *weight);

class MEWeight{
  public:

  double Integrand(const double* Xarg, const double *weight);
  double ComputePdf(const int &pid, const double &x, const double &q2);
  inline void setProcessMomenta(vector<double*> &p){ process.setMomenta(p); }
  inline void computeMatrixElements(){ process.sigmaKin(); }
  inline const double* const getMatrixElements() const { return process.getMatrixElements(); }
  double ComputeWeight(double &error);
  MEEvent* GetEvent();
  void SetEvent(const ROOT::Math::PtEtaPhiEVector &ep, const ROOT::Math::PtEtaPhiEVector &mum, const ROOT::Math::PtEtaPhiEVector &b, const ROOT::Math::PtEtaPhiEVector &bbar, const ROOT::Math::PtEtaPhiEVector &met);
  void AddTF(const std::string particleName, const std::string histName);

  MEWeight(const std::string paramCardPath, const std::string pdfName, const std::string fileTF);
  ~MEWeight();

  private:

  CPPProcess process;
  LHAPDF::PDF* pdf;
  MEEvent* myEvent;
  TransferFunction* myTF;
};

#endif
