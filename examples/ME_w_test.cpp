

#include "TLorentzVector.h"
#include "TMath.h"

#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "Riostream.h"

#include "/home/fynu/amertens/scratch/MadGraph/madgraph5/w_cpp/src/HelAmps_sm.h"
#include "/home/fynu/amertens/scratch/MadGraph/madgraph5/w_cpp/SubProcesses/P0_Sigma_sm_udx_epve/CPPProcess.h"
#include "/home/fynu/amertens/scratch/MadGraph/madgraph5/w_cpp/src/rambo.h"
#include "/home/fynu/amertens/scratch/MadGraph/madgraph5/w_cpp/src/Parameters_sm.h"

#include "TStopwatch.h"

#include "TString.h"
#include "TApplication.h"
#include "TChain.h"
#include "TFile.h"
#include "TObject.h"

#include "TClonesArray.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TRandom3.h"

class TFDISTR: public TFoamIntegrand {
private:
  CPPProcess process;

public:
  TFDISTR(std::string paramCardPath){
  process.initProc(paramCardPath);
  }


  Double_t Density(int nDim, Double_t *Xarg){
  // Integrand for mFOAM

  double weight;

  Double_t ENu=150*Xarg[0];
  Double_t EEl=150*Xarg[1];


  TLorentzVector p1,p2;
  p1.SetPtEtaPhiM(EEl,0,0,0);
  p2.SetPtEtaPhiM(ENu,0,TMath::Pi(),0);

  Double_t gEcm=(EEl+ENu)/2.0;

  // momentum vector definition
  vector<double*> p(1, new double[4]);
  p[0][0] = gEcm; p[0][1] = 0.0; p[0][2] = 0.0; p[0][3] = gEcm;
  p.push_back(new double[4]);
  p[1][0] = gEcm; p[1][1] = 0.0; p[1][2] = 0.0; p[1][3] = -gEcm;
  p.push_back(new double[4]);
  p[2][0] = p1.E(); p[2][1] = p1.Px(); p[2][2] = p1.Py(); p[2][3] = p1.Pz();
  p.push_back(new double[4]);
  p[3][0] = p2.E(); p[3][1] = p2.Px(); p[3][2] = p2.Py(); p[3][3] = p2.Pz();


  // Set momenta for this event
  process.setMomenta(p);

  // Evaluate matrix element
  process.sigmaKin();

  const Double_t* matrix_elements = process.getMatrixElements();
/*
  std::cout << "Momenta:" << std::endl;

  for(int i=0;i < process.nexternal; i++)
    std::cout << setw(4) << i+1
         << setiosflags(ios::scientific) << setw(14) <<" " << p[i][0]
         << setiosflags(ios::scientific) << setw(14) <<" " << p[i][1]
         << setiosflags(ios::scientific) << setw(14) <<" " << p[i][2]
         << setiosflags(ios::scientific) << setw(14) <<" " << p[i][3] << std::endl;
  std::cout << " -----------------------------------------------------------------------------" << std::endl;

  // Display matrix elements
  for(int i=0; i<process.nprocesses;i++)
    std::cout << " Matrix element = "
         << setiosflags(ios::fixed) << setprecision(17)
         << matrix_elements[i]
         << " GeV^" << -(2*process.nexternal-8) << std::endl;

  std::cout << " -----------------------------------------------------------------------------" << std::endl;
  */
  return matrix_elements[0];
  }


  //ClassDef(TFDISTR,1) //Class of testing functions for FOAM
};
//ClassImp(TFDISTR)

int main(int argc, char** argv){


  //process.initProc("/home/fynu/amertens/scratch/MadGraph/madgraph5/w_cpp/Cards/param_card.dat");

  //TH2D *hst = new TH2D("test","test",150,0,150,100,-50,0);
  TH2D *hst_Enu = new TH2D("test_2D", "test_2D", 150,0,150,150,0,150);

  TStopwatch chrono;
  std::cout << "initialising : " ;
  //chrono->Start();
  // TFoam Implementation
  Double_t *MCvect =new Double_t[2];
  TRandom3  *PseRan   = new TRandom3();
  PseRan->SetSeed(4357);  
  TFoam   *FoamX    = new TFoam("FoamX");
  FoamX->SetkDim(2);
  FoamX->SetnCells(500);      // No. of cells, can be omitted, default=2000
//  FoamX->SetRhoInt(Camel2);   // Set 2-dim distribution, included below
  TFoamIntegrand *rho= new TFDISTR("/home/fynu/amertens/scratch/MadGraph/madgraph5/w_cpp/Cards/param_card.dat");
  //rho->SetParamCard("/home/fynu/amertens/scratch/MadGraph/madgraph5/w_cpp/Cards/param_card.dat");
  FoamX->SetRho(rho);
  FoamX->SetPseRan(PseRan);   // Set random number generator
  FoamX->Initialize(); 
  
  std::cout << "CPU time : " << chrono.CpuTime() << "  Real-time : " << chrono.RealTime() << std::endl;
  chrono.Reset();
  chrono.Start();
//  std::cout << std::endl;
//  std::cout << "Looping over phase-space points" << std::endl;

//  chrono->Reset();
//  chrono->Start();

  for(Long_t loop=0; loop<100000; loop++){
    FoamX->MakeEvent();          // generate MC event
    FoamX->GetMCvect( MCvect);   // get generated vector (x,y)
    Double_t x=MCvect[0];
    Double_t y=MCvect[1];
    //if(loop<10) cout<<"(x,y) =  ( "<< x <<", "<< y <<" )"<<endl;
    hst_Enu->Fill(150*x, 150*y);           // fill scattergram
   }// loop

  Double_t mcResult, mcError;
  FoamX->GetIntegMC( mcResult, mcError);  // get MC integral, should be one
  
  std::cout << "CPU time : " << chrono.CpuTime() << "  Real-time : " << chrono.RealTime() << std::endl;
//  chrono->Print();
//  std::cout << std::endl;

  cout << " mcResult= " << mcResult << " +- " << mcError <<endl;
  // now hst_xy will be plotted visualizing generated distribution
  TCanvas *c = new TCanvas("c","Canvas for plotting",600,600);
   c->cd();
   hst_Enu->Draw();
   c->Print("Enu.png");
}

/*
Double_t MatrixElement(Int_t nDim, Double_t *Xarg ){


  CPPProcess process;

  // Read param_card and set parameters
  process.initProc("/home/fynu/amertens/scratch/MadGraph/madgraph5/w_cpp/Cards/param_card.dat");

  double weight;

  Double_t EtaNu=Xarg[0];
  Double_t ENu=Xarg[1];


  TLorentzVector p1,p2;
  p1.SetPtEtaPhiM(40,0,0,0);
  p2.SetPtEtaPhiM(ENu,EtaNu,TMath::Pi(),0);

  Double_t gEcm=80;

  // momentum vector definition
  vector<double*> p(1, new double[4]);
  p[0][0] = gEcm; p[0][1] = 0.0; p[0][2] = 0.0; p[0][3] = gEcm;
  p.push_back(new double[4]);
  p[1][0] = gEcm; p[1][1] = 0.0; p[1][2] = 0.0; p[1][3] = -gEcm;
  p.push_back(new double[4]);
  p[2][0] = p1.E(); p[2][1] = p1.Px(); p[2][2] = p1.Py(); p[2][3] = p1.Pz();
  p.push_back(new double[4]);
  p[3][0] = p2.E(); p[3][1] = p2.Px(); p[3][2] = p2.Py(); p[3][3] = p2.Pz();


  // Set momenta for this event
  process.setMomenta(p);

  // Evaluate matrix element
  process.sigmaKin();

  const Double_t* matrix_elements = process.getMatrixElements();

  std::cout << "Momenta:" << std::endl;
  for(int i=0;i < process.nexternal; i++)
    std::cout << setw(4) << i+1
         << setiosflags(ios::scientific) << setw(14) << p[i][0]
         << setiosflags(ios::scientific) << setw(14) << p[i][1]
         << setiosflags(ios::scientific) << setw(14) << p[i][2]
         << setiosflags(ios::scientific) << setw(14) << p[i][3] << std::endl;
  std::cout << " -----------------------------------------------------------------------------" << std::endl;

  // Display matrix elements
  for(int i=0; i<process.nprocesses;i++)
    std::cout << " Matrix element = "
         << setiosflags(ios::fixed) << setprecision(17)
         << matrix_elements[i]
         << " GeV^" << -(2*process.nexternal-8) << std::endl;

  std::cout << " -----------------------------------------------------------------------------" << std::endl;
  
  return matrix_elements;
}
*/

