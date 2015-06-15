#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdlib.h>

#include "classes/DelphesClasses.h"

#include "src/HelAmps_sm.h"
#include "SubProcesses/P0_Sigma_sm_gg_epvebmumvmxbx/CPPProcess.h"
#include "src/rambo.h"
#include "src/Parameters_sm.h"

#include "TStopwatch.h"
#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"

#include "cuba.h"

#include "utils.h"
#include "jacobianD.h"

#define M_T 173.
#define G_T 1.4915

#define M_W 80.419
#define G_W 2.0476

#define SQRT_S 13000

#define BINNING 1750
#define START 250
#define STOP 2000

#define SSTR( x ) dynamic_cast< std::ostringstream & > \
        ( std::ostringstream() << std::dec << x ).str()

using namespace LHAPDF;
using namespace std;

int mycount = 0, count_wgt = 0, count_perm=1;

int MEFunct(const int *nDim, const double* Xarg, const int *nComp, double *Value, void *inputs, const int *nVec, const int *core, const double *weight);

unsigned int setFlags(char verbosity = 0, bool subregion = false, bool retainStateFile = false, unsigned int level = 0, bool smoothing = false, bool takeOnlyGridFromFile = true){
  // Set option flags for CUBA integrator
  // 
  // smoothing only used by Suave
  // smoothing and takeOnlyGridFromFile only used by Vegas
  //
  // verbosity = 0-3
  // subregion true = only last set of samples is used for final evaluation of integral
  // smoothing true = smoothe importance function
  // retainStateFile true => retain state file when integration ends
  // takeOnlyGridFromFile false => full state taken from file (if present), true => only grid is taken (e.g. to use it for another integrand)
  // level determines random-number generator:
  //  seed = 0 => Sobol quasi-random
  //  seed != 0 => level is used:
  //    level = 0 => Mersenne Twister pseudo-random
  //    level != 0 => Ranlux pseudo-random, level = determines period of generator
  //

  unsigned int flags = 0;

  unsigned int opt_subregion = 0x04; // bit 2
  unsigned int opt_smoothing = 0x08; // bit 3
  unsigned int opt_retainStateFile = 0x10; // bit 4
  unsigned int opt_takeOnlyGridFromFile = 0x20; // bit 5

  level <<= 8; // bits 8-31
  flags |= level | verbosity; // verbosity: bis 0-1
  if(subregion) flags |= opt_subregion;
  if(!smoothing) flags |= opt_smoothing; // careful true-false inverted
  if(retainStateFile) flags |= opt_retainStateFile;
  if(takeOnlyGridFromFile) flags |= opt_takeOnlyGridFromFile;

  cout << "Integrator flags = ";
  for(int i=31; i>=0; --i){
    bool bit = (flags >> i) & 1;
    cout << bit;
  }
  cout << endl;

  return flags;
}


class MEEvent{
  public:

  MEEvent();
  void SetVectors(const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met);
  ~MEEvent();

  inline TLorentzVector GetP3(void) const { return p3; }
  inline TLorentzVector GetP4(void) const { return p4; }
  inline TLorentzVector GetP5(void) const { return p5; }
  inline TLorentzVector GetP6(void) const { return p6; }
  inline TLorentzVector GetMet(void) const { return Met; }

  void writeHists(void);

  TH1D* GetTTbar();
  
  private:

  TLorentzVector p3, p4, p5, p6, Met;

  TH1D *hst_TTbar;
};

MEEvent::MEEvent(){
  hst_TTbar = new TH1D("test_TTbar", "test_TTbar", BINNING, START, STOP); 
}

void MEEvent::SetVectors(const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met){
  p3 = ep;
  p5 = mum;
  p4 = b;
  p6 = bbar;
  Met = met;

  if(hst_TTbar)
    hst_TTbar->Reset();
}

MEEvent::~MEEvent(){
  delete hst_TTbar; hst_TTbar = NULL;
}

TH1D* MEEvent::GetTTbar(){
  return hst_TTbar;
}

void MEEvent::writeHists(void){
  TCanvas *c = new TCanvas(TString("Mttbar")+"_"+SSTR(count_wgt)+"_"+SSTR(count_perm),"Canvas for plotting",600,600);
  c->cd();
  hst_TTbar->Draw();
  c->Write();
  c->Print(TString("plots/Mttbar_")+SSTR(count_wgt)+"_"+SSTR(count_perm)+".png");
  delete c; c = 0;
}

class MEWeight{
  public:

  inline double ComputePdf(const int pid, const double x, const double q2);
  inline void setProcessMomenta(vector<double*> &p){ process.setMomenta(p); }
  inline void computeMatrixElements(){ process.sigmaKin(); }
  inline const double* const getMatrixElements() const { return process.getMatrixElements(); }
  double ComputeWeight(double &error);
  MEEvent* GetEvent();
  void SetEvent(const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met);
  double GetTempAverage();
  double GetTempMaxLikelihood();
  void AddAndResetTempHist();

  void WriteHist();

  MEWeight(const string paramCardPath, const string pdfName);
  ~MEWeight();

  private:

  TH1D* hst_TTbar;
  TH1D* hst_TTbar_temp;

  CPPProcess process;
  PDF* pdf;
  MEEvent* myEvent;
}; 
  
MEWeight::MEWeight(const string paramCardPath, const string pdfName){

  cout << "Initizialing Matrix Element computation with:" << endl;
  cout << "Parameter card " << paramCardPath << endl;
  cout << "PDF " << pdfName << endl;

  process.initProc(paramCardPath);
  pdf = mkPDF(pdfName, 0);
  myEvent = new MEEvent();

  hst_TTbar = new TH1D("DMEM_TTbar", "DMEM  M_{tt}", BINNING, START, STOP);
  hst_TTbar_temp = new TH1D("DMEM_TTbar_temp", "DMEM  M_{tt}", BINNING, START, STOP);
}

double MEWeight::ComputePdf(const int pid, const double x, const double q2){
  // return f(pid,x,q2)
  if(x <= 0 || x >= 1 || q2 <= 0){
    cout << "WARNING: PDF x or Q^2 value out of bounds!" << endl;
    return 0.;
  }else{
    return pdf->xfxQ2(pid, x, q2)/x;
  }
}

MEEvent* MEWeight::GetEvent(){
  return myEvent;
}

void MEWeight::SetEvent(const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met){
  myEvent->SetVectors(ep, mum, b, bbar, met);
}

void MEWeight::WriteHist(){

  hst_TTbar->Scale(1./hst_TTbar->Integral());

  /*TCanvas *c = new TCanvas("DMEM_TTbar", "Canvas for plotting", 600, 600);
  c->cd();
  hst_TTbar->Draw();*/
  //c->Write();
  hst_TTbar->Write();
  //c->Print("plots/DMEM_ttbar.png", "png");
  //delete c; c = 0;

}

double MEWeight::GetTempAverage(){
  return hst_TTbar_temp->GetMean();
}

double MEWeight::GetTempMaxLikelihood(){
  return hst_TTbar_temp->GetBinCenter(hst_TTbar_temp->GetMaximumBin());
}

void MEWeight::AddAndResetTempHist(){
  if(hst_TTbar_temp->Integral() != 0.){
    hst_TTbar_temp->Scale(1./hst_TTbar_temp->Integral());
    hst_TTbar_temp->SetEntries(1);
    hst_TTbar->Add(hst_TTbar_temp);
  }
  hst_TTbar_temp->Reset();
}

double MEWeight::ComputeWeight(double &error){
  
  cout << "Initializing integration..." << endl;

  int neval, nfail;
  double mcResult=0, prob=0;
  
  char verbosity = 0; // 0-3
  bool subregion = false; // true = only last set of samples is used for final evaluation of integral
  bool smoothing = false;
  bool retainStateFile = false; // false => delete state file when integration ends
  bool takeOnlyGridFromFile = true; // false => full state taken from file (if present), true => only grid is taken (e.g. to use it for another integrand)
  unsigned int level = 0; 

  unsigned int flags = setFlags(verbosity, subregion, retainStateFile, level, smoothing, takeOnlyGridFromFile);

  cout << "Starting integration..." << endl;

  cubacores(0,0);           // This is mandatory if the integrand wants to *modify* something in the MEWeight object passed as argument
  Vegas(
    4,                      // (int) dimensions of the integrated volume
    1,                      // (int) dimensions of the integrand
    (integrand_t) MEFunct,  // (integrand_t) integrand (cast to integrand_t)
    (void*) this,           // (void*) pointer to additional arguments passed to integrand
    1,                      // (int) maximum number of points given the integrand in each invocation (=> SIMD) ==> PS points = vector of sets of points (x[ndim][nvec]), integrand returns vector of vector values (f[ncomp][nvec])
    0.005,                  // (double) requested relative accuracy |-> error < max(rel*value,abs)
    0.,                     // (double) requested absolute accuracy |
    flags,                  // (int) various control flags in binary format, see setFlags function
    0,                      // (int) seed (seed==0 => SOBOL; seed!=0 && control flag "level"==0 => Mersenne Twister)
    0,                  // (int) minimum number of integrand evaluations
    /*1,                  // (int) maximum number of integrand evaluations (approx.!)
    1,                  // (int) number of integrand evaluations per interations (to start)
    0,                      // (int) increase in number of integrand evaluations per interations
    1,                   // (int) batch size for sampling*/
    750000,                  // (int) maximum number of integrand evaluations (approx.!)
    50000,                  // (int) number of integrand evaluations per interations (to start)
    0,                      // (int) increase in number of integrand evaluations per interations
    10000,                   // (int) batch size for sampling
    0,                      // (int) grid number, 1-10 => up to 10 grids can be stored, and re-used for other integrands (provided they are not too different)
    "",                     // (char*) name of state file => state can be stored and retrieved for further refinement
    NULL,                   // (int*) "spinning cores": -1 || NULL <=> integrator takes care of starting & stopping child processes (other value => keep or retrieve child processes, probably not useful here)
    &neval,                 // (int*) actual number of evaluations done
    &nfail,                 // 0=desired accuracy was reached; -1=dimensions out of range; >0=accuracy was not reached
    &mcResult,              // (double*) integration result ([ncomp])
    &error,                 // (double*) integration error ([ncomp])
    &prob                   // (double*) Chi-square p-value that error is not reliable (ie should be <0.95) ([ncomp])
  );
  
  cout << "Integration done." << endl;

  cout << "Status: " << nfail << ", nr. fails: " << mycount << endl;
  mycount = 0;

  cout << " mcResult= " << mcResult << " +- " << error << " in " << neval << " evaluations. Chi-square prob. = " << prob << endl;

  if(myEvent->GetTTbar()->Integral()){
    myEvent->GetTTbar()->Scale(mcResult/myEvent->GetTTbar()->Integral());
    myEvent->GetTTbar()->SetEntries(1);
    hst_TTbar_temp->Add(myEvent->GetTTbar());
  }
  
  myEvent->writeHists();

  if(std::isnan(error))
  error = 0.;
  if(std::isnan(mcResult))
  mcResult = 0.;
  return mcResult;
}

MEWeight::~MEWeight(){
  cout << "Deleting PDF" << endl;
  delete pdf; pdf = NULL;
  cout << "Deleting myEvent" << endl;
  delete myEvent; myEvent = NULL;
  cout << "Deleting hst_TTbar" << endl;
  delete hst_TTbar; hst_TTbar = NULL;
}

int MEFunct(const int *nDim, const double* Xarg, const int *nComp, double *Value, void *inputs, const int *nVec, const int *core, const double *weight){
  //cout << endl << endl << endl << "########## Starting phase-space point ############" << endl << endl;

  //cout << "Inputs = [" << Xarg[0] << "," << Xarg[1] << "," << Xarg[2] << "," << Xarg[3] << endl;

  MEWeight* myWeight = (MEWeight*) inputs;
  MEEvent* myEvent = myWeight->GetEvent();

  TLorentzVector p3 = myEvent->GetP3();
  TLorentzVector p4 = myEvent->GetP4();
  TLorentzVector p5 = myEvent->GetP5();
  TLorentzVector p6 = myEvent->GetP6();
  TLorentzVector Met = myEvent->GetMet();

  *Value = 0.;

  for(int i=0; i<*nDim; ++i){
    if(Xarg[i] == 1.){
      mycount++;
      return 0;
    }
  }

  // We flatten the Breit-Wigners by doing a change of variable for each resonance separately
  // s = M G tan(y) + M^2
  // jac = M G / cos^2(y)
  // ==> BW(s(y))*jac(y) is flat in the variable y, as BW(s) = 1/((s-M^2)^2 - (GM)^2)
  // Where y = -arctan(M/G) + (pi/2+arctan(M/G))*x_foam (x_foam between 0 and 1 => s between 0 and infinity)

  const double range1 = TMath::Pi()/2. + TMath::ATan(M_W/G_W);
  const double y1 = - TMath::ATan(M_W/G_W) + range1 * Xarg[0];
  const double s13 = M_W * G_W * TMath::Tan(y1) + pow(M_W,2.);

  //cout << "y1=" << y1 << ", m13=" << TMath::Sqrt(s13) << endl;

  const double range2 = TMath::Pi()/2. + TMath::ATan(M_T/G_T);
  const double y2 = - TMath::ATan(M_T/G_T) + range2 * Xarg[1];
  const double s134 = M_T * G_T * TMath::Tan(y2) + pow(M_T,2.);

  //cout << "y2=" << y2 << ", m134=" << TMath::Sqrt(s134) << endl;

  const double range3 = TMath::Pi()/2. + TMath::ATan(M_W/G_W);
  const double y3 = - TMath::ATan(M_W/G_W) + range3 * Xarg[2];
  const double s25 = M_W * G_W * TMath::Tan(y3) + pow(M_W,2.);

  //cout << "y3=" << y3 << ", m25=" << TMath::Sqrt(s25) << endl;

  const double range4 = TMath::Pi()/2. + TMath::ATan(M_T/G_T);
  const double y4 = - TMath::ATan(M_T/G_T) + range4 * Xarg[3];
  const double s256 = M_T * G_T * TMath::Tan(y4) + pow(M_T,2.);

  //cout << "y4=" << y4 << ", m256=" << TMath::Sqrt(s256) << endl;
  
  double flatterJac = range1 * range2 * range3 * range4;
  flatterJac *= M_W*G_W * M_T*G_T * M_W*G_W * M_T*G_T;
  flatterJac /= pow(TMath::Cos(y1) * TMath::Cos(y2) * TMath::Cos(y3) * TMath::Cos(y4), 2.);

  // #### FOR NWA
  /*const double s13 = pow(M_W,2.);
  const double s134 = pow(M_T,2.);
  const double s25 = pow(M_W,2.);
  const double s256 = pow(M_T,2.);*/

  if(s13 > s134 || s25 > s256 || s13 < p3.M() || s25 < p5.M() || s134 < p4.M() || s256 < p6.M()){
    //cout << "Masses too small!" << endl;
    mycount++;
    return 0;
  }

  //cout << "weight = " << *weight << endl;

  // pT = transverse total momentum of the visible particles
  //TLorentzVector pT = p3 + p4 + p5 + p6;
  const TLorentzVector pT = -Met;

  const double p34 = p3*p4;
  const double p56 = p5*p6;
  const double p33 = p3.M2();
  const double p44 = p4.M2();
  const double p55 = p5.M2();
  const double p66 = p6.M2();

  //cout << "p34=" << p34 << ",p56=" << p56 << ",p33=" << p33 << ",p44=" << p44 << ",p55" << p55 << ",p66" << p66;
  //cout << "Met=" << Met.Pt() << endl;

  // A1 p1x + B1 p1y + C1 = 0, with C1(E1,E2)
  // A2 p1y + B2 p2y + C2 = 0, with C2(E1,E2)
  // ==> express p1x and p1y as functions of E1, E2
  
  const double A1 = 2.*( -p3.Px() + p3.Pz()*p4.Px()/p4.Pz() );
  const double A2 = 2.*( p5.Px() - p5.Pz()*p6.Px()/p6.Pz() );
  
  const double B1 = 2.*( -p3.Py() + p3.Pz()*p4.Py()/p4.Pz() );
  const double B2 = 2.*( p5.Py() - p5.Pz()*p6.Py()/p6.Pz() );

  const double Dx = B2*A1 - B1*A2;
  const double Dy = A2*B1 - A1*B2;

  const double X = 2*( pT.Px()*p5.Px() + pT.Py()*p5.Py() - p5.Pz()/p6.Pz()*( 0.5*(s25 - s256 + p66) + p56 + pT.Px()*p6.Px() + pT.Py()*p6.Py() ) ) + p55 - s25;
  const double Y = p3.Pz()/p4.Pz()*( s13 - s134 + 2*p34 + p44 ) - p33 + s13;

  // p1x = alpha1 E1 + beta1 E2 + gamma1
  // p1y = ...(2)
  // p1z = ...(3)
  // p2z = ...(4)
  // p2x = ...(5)
  // p2y = ...(6)
  
  const double alpha1 = -2*B2*(p3.E() - p4.E()*p3.Pz()/p4.Pz())/Dx;
  const double beta1 = 2*B1*(p5.E() - p6.E()*p5.Pz()/p6.Pz())/Dx;
  const double gamma1 = B1*X/Dx + B2*Y/Dx;

  const double alpha2 = -2*A2*(p3.E() - p4.E()*p3.Pz()/p4.Pz())/Dy;
  const double beta2 = 2*A1*(p5.E() - p6.E()*p5.Pz()/p6.Pz())/Dy;
  const double gamma2 = A1*X/Dy + A2*Y/Dy;

  const double alpha3 = (p4.E() - alpha1*p4.Px() - alpha2*p4.Py())/p4.Pz();
  const double beta3 = -(beta1*p4.Px() + beta2*p4.Py())/p4.Pz();
  const double gamma3 = ( 0.5*(s13 - s134 + p44) + p34 - gamma1*p4.Px() - gamma2*p4.Py() )/p4.Pz();

  const double alpha4 = (alpha1*p6.Px() + alpha2*p6.Py())/p6.Pz();
  const double beta4 = (p6.E() + beta1*p6.Px() + beta2*p6.Py())/p6.Pz();
  const double gamma4 = ( 0.5*(s25 - s256 + p66) + p56 + (gamma1 + pT.Px())*p6.Px() + (gamma2 + pT.Py())*p6.Py() )/p6.Pz();

  const double alpha5 = -alpha1;
  const double beta5 = -beta1;
  const double gamma5 = -pT.Px() - gamma1;

  const double alpha6 = -alpha2;
  const double beta6 = -beta2;
  const double gamma6 = -pT.Py() - gamma2;

  // a11 E1^2 + a22 E2^2 + a12 E1E2 + a10 E1 + a01 E2 + a00 = 0
  // id. with bij

  const double a11 = -1 + ( pow(alpha1,2.) + pow(alpha2,2.) + pow(alpha3,2.) );
  const double a22 = pow(beta1,2.) + pow(beta2,2.) + pow(beta3,2.);
  const double a12 = 2.*( alpha1*beta1 + alpha2*beta2 + alpha3*beta3 );
  const double a10 = 2.*( alpha1*gamma1 + alpha2*gamma2 + alpha3*gamma3 );
  const double a01 = 2.*( beta1*gamma1 + beta2*gamma2 + beta3*gamma3 );
  const double a00 = pow(gamma1,2.) + pow(gamma2,2.) + pow(gamma3,2.);

  const double b11 = pow(alpha5,2.) + pow(alpha6,2.) + pow(alpha4,2.);
  const double b22 = -1 + ( pow(beta5,2.) + pow(beta6,2.) + pow(beta4,2.) );
  const double b12 = 2.*( alpha5*beta5 + alpha6*beta6 + alpha4*beta4 );
  const double b10 = 2.*( alpha5*gamma5 + alpha6*gamma6 + alpha4*gamma4 );
  const double b01 = 2.*( beta5*gamma5 + beta6*gamma6 + beta4*gamma4 );
  const double b00 = pow(gamma5,2.) + pow(gamma6,2.) + pow(gamma4,2.);

  // Find the intersection of the 2 conics (at most 4 real solutions for (E1,E2))
  vector<double> E1, E2;
  //cout << "coefs=" << a11 << "," << a22 << "," << a12 << "," << a10 << "," << a01 << "," << a00 << endl;
  //cout << "coefs=" << b11 << "," << b22 << "," << b12 << "," << b10 << "," << b01 << "," << b00 << endl;
  solve2Quads(a11, a22, a12, a10, a01, a00, b11, b22, b12, b10, b01, b00, E1, E2, false);

  // For each solution (E1,E2), find the neutrino 4-momenta p1,p2, find the initial quark momenta,
  // evaluate the matrix element and the jacobian

  int countSol = 0;
  
  if(E1.size() == 0){
    //cout << "No solutions!" << endl;
    mycount++;
    return 0;
  }

  /*myWeight->fillWe(TMath::Sqrt(s13));
  myWeight->fillWm(TMath::Sqrt(s25));
  myWeight->fillT(TMath::Sqrt(s134));
  myWeight->fillTbar(TMath::Sqrt(s256));*/

  for(unsigned int i=0; i<E1.size(); i++){
    /*if(E1.size() >= 0){
      if(E1.at(0) > E1.at(1))
        i = 0;
      else
        i = 1;
    }*/

    const double e1 = E1.at(i);
    const double e2 = E2.at(i);

    //cout << endl << "## Evaluating Matrix Element based on solutions e1 = " << e1 << ", e2 = " << e2 << endl << endl;

    if(e1 < 0. || e2 < 0.){
      //mycount++;
      //cout << "Neg energies." << endl;
      continue;
      //break;
    }

    TLorentzVector p1,p2;

    p1.SetPx( alpha1*e1 + beta1*e2 + gamma1 );
    p1.SetPy( alpha2*e1 + beta2*e2 + gamma2 );
    p1.SetPz( alpha3*e1 + beta3*e2 + gamma3 );
    p1.SetE(e1);

    p2.SetPx( alpha5*e1 + beta5*e2 + gamma5 );
    p2.SetPy( alpha6*e1 + beta6*e2 + gamma6 );
    p2.SetPz( alpha4*e1 + beta4*e2 + gamma4 ); 
    p2.SetE(e2);

    const TLorentzVector p13 = p1 + p3;
    const TLorentzVector p134 = p1 + p3 + p4;
    const TLorentzVector p25 = p2 + p5;
    const TLorentzVector p256 = p2 + p5 + p6;

    /*myWeight->hst_We->Fill(p13.M());
    myWeight->hst_Wm->Fill(p25.M());
    myWeight->hst_t->Fill(p134.M());
    myWeight->hst_tbar->Fill(p256.M());*/
    /*cout << "Solution " << i << ":" << endl;
    cout << "Input: W+ mass=" << TMath::Sqrt(s13) << ", Top mass=" << TMath::Sqrt(s134) << ", W- mass=" << TMath::Sqrt(s25) << ", Anti-top mass=" << TMath::Sqrt(s256) << endl;
    cout << "Output: W+ mass=" << p13.M() << ", Top mass=" << p134.M() << ", W- mass=" << p25.M() << ", Anti-top mass=" << p256.M() << endl << endl;
    //cout << "Differences: W+ mass=" << TMath::Sqrt(s13)-p13.M() << ", Top mass=" << TMath::Sqrt(s134)-p134.M() << ", W- mass=" << TMath::Sqrt(s25)-p25.M() << ", Anti-top mass=" << TMath::Sqrt(s256)-p256.M() << endl << endl;
    cout << "Accuracies: " << endl;
    cout << "W+ mass=" << (TMath::Sqrt(s13)-p13.M())/p13.M() << endl;
    cout << "Top mass=" << (TMath::Sqrt(s134)-p134.M())/p134.M() << endl;
    cout << "W- mass=" << (TMath::Sqrt(s25)-p25.M())/p25.M() << endl;
    cout << "Anti-top mass=" << (TMath::Sqrt(s256)-p256.M())/p256.M() << endl << endl;*/
    
    /*cout << "Electron (E,Px,Py,Pz) = ";
    cout << p3.E() << "," << p3.Px() << "," << p3.Py() << "," << p3.Pz() << endl;
    cout << "Electron neutrino (E,Px,Py,Pz) = ";
    cout << p1.E() << "," << p1.Px() << "," << p1.Py() << "," << p1.Pz() << endl;
    cout << "b quark (E,Px,Py,Pz) = ";
    cout << p4.E() << "," << p4.Px() << "," << p4.Py() << "," << p4.Pz() << endl;
    cout << "Muon (E,Px,Py,Pz) = ";
    cout << p5.E() << "," << p5.Px() << "," << p5.Py() << "," << p5.Pz() << endl;
    cout << "Muon neutrino (E,Px,Py,Pz) = ";
    cout << p2.E() << "," << p2.Px() << "," << p2.Py() << "," << p2.Pz() << endl;
    cout << "Anti b quark (E,Px,Py,Pz) = ";
    cout << p6.E() << "," << p6.Px() << "," << p6.Py() << "," << p6.Pz() << endl << endl;*/
  
    const TLorentzVector tot = p1 + p2 + p3 + p4 + p5 + p6;
    
    const double ETot = tot.E();
    const double PzTot = tot.Pz();

    const double q1Pz = (PzTot + ETot)/2.;
    const double q2Pz = (PzTot - ETot)/2.;

    //cout << "===> Eext=" << ETot << ", Pzext=" << PzTot << ", q1Pz=" << q1Pz << ", q2Pz=" << q2Pz << endl << endl;
  
    if(q1Pz > SQRT_S/2. || q2Pz < -SQRT_S/2. || q1Pz < 0. || q2Pz > 0.){
      //cout << "Fail!" << endl;
      mycount++;
      continue;
      //break;
    }
    
    // Compute jacobian from change of variable:
    vector<TLorentzVector> momenta;
    momenta.push_back(p1);
    momenta.push_back(p2);
    momenta.push_back(p3);
    momenta.push_back(p4);
    momenta.push_back(p5);
    momenta.push_back(p6);
    const double jac = computeJacobianD(momenta, SQRT_S);
    if(jac <= 0.){
      cout << "Jac infinite!" << endl;
      mycount++;
      continue;
      //break;
    }
  
    // Compute the Pdfs
    const double x1 = TMath::Abs(q1Pz/(SQRT_S/2.));
    const double x2 = TMath::Abs(q2Pz/(SQRT_S/2.));

    const double pdf1_1 = myWeight->ComputePdf(21, x1, pow(M_T,2));
    const double pdf1_2 = myWeight->ComputePdf(21, x2, pow(M_T,2));
  
    // Compute flux factor 1/(2*x1*x2*s)
    const double PhaseSpaceIn = 1.0 / ( 2. * x1 * x2 * pow(SQRT_S,2)); 

    // Compute finale Phase Space for observed particles (not concerned by the change of variable)
    // dPhi = |P|^2 sin(theta)/(2*E*(2pi)^3)
    const double dPhip3 = pow(p3.P(),2.)*TMath::Sin(p3.Theta())/(2.0*p3.E()*pow(2.*TMath::Pi(),3));
    const double dPhip4 = pow(p4.P(),2.)*TMath::Sin(p4.Theta())/(2.0*p4.E()*pow(2.*TMath::Pi(),3));
    const double dPhip5 = pow(p5.P(),2.)*TMath::Sin(p5.Theta())/(2.0*p5.E()*pow(2.*TMath::Pi(),3));
    const double dPhip6 = pow(p6.P(),2.)*TMath::Sin(p6.Theta())/(2.0*p6.E()*pow(2.*TMath::Pi(),3));
    const double PhaseSpaceOut = dPhip5 * dPhip6 * dPhip3 * dPhip4;

    // momentum vector definition
    vector<double*> p;
    p.push_back(new double[4]);
    p[0][0] = q1Pz; p[0][1] = 0.0; p[0][2] = 0.0; p[0][3] = q1Pz;
    p.push_back(new double[4]);
    p[1][0] = TMath::Abs(q2Pz); p[1][1] = 0.0; p[1][2] = 0.0; p[1][3] = q2Pz;
    p.push_back(new double[4]);
    p[2][0] = p3.E(); p[2][1] = p3.Px(); p[2][2] = p3.Py(); p[2][3] = p3.Pz();
    p.push_back(new double[4]);
    p[3][0] = p1.E(); p[3][1] = p1.Px(); p[3][2] = p1.Py(); p[3][3] = p1.Pz();
    p.push_back(new double[4]);
    p[4][0] = p4.E(); p[4][1] = p4.Px(); p[4][2] = p4.Py(); p[4][3] = p4.Pz();
    p.push_back(new double[4]);
    p[5][0] = p5.E(); p[5][1] = p5.Px(); p[5][2] = p5.Py(); p[5][3] = p5.Pz();
    p.push_back(new double[4]);
    p[6][0] = p2.E(); p[6][1] = p2.Px(); p[6][2] = p2.Py(); p[6][3] = p2.Pz();
    p.push_back(new double[4]);
    p[7][0] = p6.E(); p[7][1] = p6.Px(); p[7][2] = p6.Py(); p[7][3] = p6.Pz();

    // Set momenta for this event
    myWeight->setProcessMomenta(p);

    // Evaluate matrix element
    myWeight->computeMatrixElements();
    const double* const matrix_elements1 = myWeight->getMatrixElements();
    
    // free up memory
    for(unsigned int j = 0; j < p.size(); ++j){
      delete [] p.at(j); p.at(j) = NULL;
    }

    const double thisSolResult = PhaseSpaceIn * matrix_elements1[0] * pdf1_1 * pdf1_2 * PhaseSpaceOut * jac * flatterJac;
    *Value += thisSolResult; 
    
    // Check whether the next solutions for (E1,E2) are the same => don't redo all this!
    int countEqualSol = 1;
    for(unsigned int j = i+1; j<E1.size(); j++){
      if(e1 == E1.at(j) && e2 == E2.at(j)){
        *Value += thisSolResult;
        countEqualSol++;
      }
    }

    //cout << "Found PDF1 = " << pdf1_1 << ", PDF2 = " << pdf1_2 << ", PS in = " << PhaseSpaceIn << ", PS out = " << PhaseSpaceOut << ", jac = " << jac << endl;
    //cout << "===> Matrix element = " << matrix_elements1[0] << ", prod = " << thisSolResult << ", multiplicity = " << countEqualSol << endl << endl; 

    // If we have included the next solutions already, skip them!
    i += countEqualSol - 1;
    countSol += countEqualSol;
    
    myEvent->GetTTbar()->Fill(tot.M(), *weight * (double)countEqualSol * thisSolResult);
    //cout << "Mttbar: " << tot.M() << endl;
  }

  if(*Value == 0.){
    mycount++;
    //cout << "Zero!" << endl;
    return 0;
  }

  return 0;
}



int main(int argc, char *argv[])
{
  TString inputFile(argv[1]);
  TString outputFile(argv[2]);
  int start_evt = atoi(argv[3]);
  int end_evt = atoi(argv[4]);

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Get pointers to branches used in this analysis
  TClonesArray *branchGen = NULL;
  chain.SetBranchAddress("Particle", &branchGen);

  cout << "Entries:" << chain.GetEntries() << endl;

  TFile* outFile = new TFile(outputFile, "RECREATE");
  TTree* outTree = chain.CloneTree(0);

  double Weight_TT_cpp, Weight_TT_Error_cpp;
  bool Weighted_TT_cpp;
  double time;
  outTree->Branch("Weight_TT_cpp", &Weight_TT_cpp);
  outTree->Branch("Weight_TT_Error_cpp", &Weight_TT_Error_cpp);
  outTree->Branch("Weighted_TT_cpp", &Weighted_TT_cpp);
  outTree->Branch("Weight_TT_cpp_time", &time);
  
  double mTruth, mAverage, mMaxL;
  outTree->Branch("MTTbar_MCTruth", &mTruth);
  outTree->Branch("MTTbar_DMEMAverage", &mAverage);
  outTree->Branch("MTTbar_DMEMMaxL", &mMaxL);

  if(end_evt >= chain.GetEntries())
    end_evt = chain.GetEntries()-1;
  
  MEWeight* myWeight = new MEWeight("/home/fynu/swertz/scratch/Madgraph/madgraph5/cpp_ttbar_epmum/Cards/param_card.dat", "cteq6l1");

  TH1D* truth_TTbar = new TH1D("MTruth_TTbar", "M_{tt}  Truth", BINNING, START, STOP);

  for(int entry = start_evt; entry <= end_evt ; ++entry){
    // Load selected branches with data from specified event
    chain.GetEntry(entry);

    TLorentzVector gen_ep, gen_mum, gen_b, gen_bbar, gen_Met;

    GenParticle *gen;

    for (int i = 0; i < branchGen->GetEntries(); i++){
      gen = (GenParticle*) branchGen->At(i);
      //cout << "Status=" << gen->Status << ", PID=" << gen->PID << ", E=" << gen->P4().E() << endl;
      if (gen->Status == 1){
        if (gen->PID == -11){
          gen_ep = gen->P4();
        }else if (gen->PID == 13){
          gen_mum = gen->P4();
        }
        else if (gen->PID == 12) gen_Met += gen->P4();
        else if (gen->PID == -14) gen_Met += gen->P4();
        else if (gen->PID == 5) gen_b = gen->P4();
        else if (gen->PID == -5) gen_bbar = gen->P4();
      }
    }

    cout << "From MadGraph:" << endl;
    cout << "Electron" << endl;
    cout << gen_ep.E() << "," << gen_ep.Px() << "," << gen_ep.Py() << "," << gen_ep.Pz() << endl;
    cout << "b quark" << endl;
    cout << gen_b.E() << "," << gen_b.Px() << "," << gen_b.Py() << "," << gen_b.Pz() << endl;
    cout << "Muon" << endl;
    cout << gen_mum.E() << "," << gen_mum.Px() << "," << gen_mum.Py() << "," << gen_mum.Pz() << endl;
    cout << "Anti b quark" << endl;
    cout << gen_bbar.E() << "," << gen_bbar.Px() << "," << gen_bbar.Py() << "," << gen_bbar.Pz() << endl;
    cout << "MET" << endl;
    cout << gen_Met.E() << "," << gen_Met.Px() << "," << gen_Met.Py() << "," << gen_Met.Pz() << endl;

    Weight_TT_cpp = 0.;
    Weight_TT_Error_cpp = 0.;
    time = 0.;
  
    TStopwatch chrono;
    chrono.Start();

    mTruth = (gen_ep + gen_b + gen_mum + gen_bbar + gen_Met).M();
    truth_TTbar->Fill(mTruth);
    cout << "M_ttbar truth: " << mTruth << endl;

    for(int permutation = 1; permutation <= 2; permutation ++){

      //permutation = 1;

      count_perm = permutation;
    
      double weight = 0;
      double error = 0;

      if(permutation == 1)
        myWeight->SetEvent(gen_ep, gen_mum, gen_b, gen_bbar, gen_Met);
      if(permutation == 2)
        myWeight->SetEvent(gen_ep, gen_mum, gen_bbar, gen_b, gen_Met);

      weight = myWeight->ComputeWeight(error)/2.; 
      /*weight = 0.;
      error = 0.;*/

      Weight_TT_cpp += weight;
      Weight_TT_Error_cpp += pow(error/2,2.);
    }
    time = chrono.CpuTime();

    Weight_TT_Error_cpp = TMath::Sqrt(Weight_TT_Error_cpp);
    Weighted_TT_cpp = true;

    cout << "====> Event " << entry << ": weight = " << Weight_TT_cpp << " +- " << Weight_TT_Error_cpp << endl;
    cout << "      CPU time : " << chrono.CpuTime() << "  Real-time : " << chrono.RealTime() << endl << endl;

    mAverage = myWeight->GetTempAverage();
    mMaxL = myWeight->GetTempMaxLikelihood();
    cout << "      Most likely M_ttbar : " << mMaxL << endl;
    cout << "      Average M_ttbar     : " << mAverage << endl << endl;
    myWeight->AddAndResetTempHist();
    outTree->Fill();

    count_wgt++;
  }

  truth_TTbar->Scale(1./truth_TTbar->Integral());
  truth_TTbar->Write();
  myWeight->WriteHist();

  outTree->Write();
  delete myWeight;
  delete truth_TTbar;
  delete outFile; outFile = NULL;
}

