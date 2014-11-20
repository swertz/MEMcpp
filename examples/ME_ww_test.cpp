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

#include "ExRootTreeReader.h"
#include "ExRootTreeWriter.h"
#include "ExRootTreeBranch.h"
#include "ExRootResult.h"
#include "DelphesClasses.h"

#include "/home/fynu/amertens/scratch/MadGraph/madgraph5/ww_cpp/src/HelAmps_sm.h"
#include "/home/fynu/amertens/scratch/MadGraph/madgraph5/ww_cpp/SubProcesses/P0_Sigma_sm_uux_epvemumvmx/CPPProcess.h"
#include "/home/fynu/amertens/scratch/MadGraph/madgraph5/ww_cpp/src/rambo.h"
#include "/home/fynu/amertens/scratch/MadGraph/madgraph5/ww_cpp/src/Parameters_sm.h"

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

#include "/home/fynu/amertens/scratch/Reconstruction/LHAPDF-6.1.4/include/LHAPDF/LHAPDF.h"
#include "/home/fynu/amertens/scratch/Reconstruction/LHAPDF-6.1.4/include/LHAPDF/PDFSet.h"


using namespace LHAPDF;


class TFDISTR: public TFoamIntegrand {
private:
  CPPProcess process;
  TLorentzVector p1, p2;
  TLorentzVector Met;

  const PDF *pdf;

public:
  TFDISTR(std::string paramCardPath, TLorentzVector ep, TLorentzVector mum, TLorentzVector met){
  process.initProc(paramCardPath);
  p1 = ep;
  p2 = mum;
  Met = met;
  pdf=LHAPDF::mkPDF("cteq6l1", 0);
  }
/*
  static void SetPDF(TString name, Int_t imem){
    *pdf = mkPDF("cteq6l1", 0);
    }
*/
  Double_t Density(int nDim, Double_t *Xarg){
  // Integrand for mFOAM

  double weight;

  Double_t Px1=-150+300*Xarg[0];
  Double_t Py1=-150+300*Xarg[1];
  Double_t Pz1=-150+300*Xarg[2];
  Double_t Px2=Met.Px()-Px1;
  Double_t Py2=Met.Py()-Py1;
  Double_t Pz2=-150+300*Xarg[3];


  Double_t E1= TMath::Sqrt(pow(Px1,2)+pow(Py1,2)+pow(Pz1,2));
  Double_t E2= TMath::Sqrt(pow(Px2,2)+pow(Py2,2)+pow(Pz2,2));

  TLorentzVector nu1,nu2;
  nu1.SetPxPyPzE(Px1,Py1,Pz1,E1);
  nu2.SetPxPyPzE(Px2,Py2,Pz2,E2);

 // Double_t gEcm=(EEl+ENu)/2.0;

  Double_t Eext=(p1+p2+nu1+nu2).E();
  Double_t Pzext=(p1+p2+nu1+nu2).Pz();

  Double_t q1Pz=(Pzext+Eext)/2;
  Double_t q2Pz=(Pzext-Eext)/2;


  

  // momentum vector definition
  vector<double*> p(1, new double[4]);
  p[0][0] = q1Pz; p[0][1] = 0.0; p[0][2] = 0.0; p[0][3] = q1Pz;
  p.push_back(new double[4]);
  p[1][0] = -q2Pz; p[1][1] = 0.0; p[1][2] = 0.0; p[1][3] = q2Pz;
  p.push_back(new double[4]);
  p[2][0] = p1.E(); p[2][1] = p1.Px(); p[2][2] = p1.Py(); p[2][3] = p1.Pz();
  p.push_back(new double[4]);
  p[3][0] = nu1.E(); p[3][1] = nu1.Px(); p[3][2] = nu1.Py(); p[3][3] = nu1.Pz();
  p.push_back(new double[4]);
  p[4][0] = p2.E(); p[4][1] = p2.Px(); p[4][2] = p2.Py(); p[4][3] = p2.Pz();
  p.push_back(new double[4]);
  p[5][0] = nu2.E(); p[5][1] = nu2.Px(); p[5][2] = nu2.Py(); p[5][3] = nu2.Pz();


  // Set momenta for this event
  process.setMomenta(p);

  // Evaluate matrix element
  process.sigmaKin();

  double pdf1 = ComputePdf(2,TMath::Abs(q1Pz/4000.0), pow(Eext,2)) / TMath::Abs(q1Pz/4000.0);
  double pdf2 = ComputePdf(-2,TMath::Abs(q2Pz/4000.0), pow(Eext,2)) / TMath::Abs(q2Pz/4000.0);

  //cout << "pdf1 : " << pdf1 << " pdf2 : " << pdf2 << endl;

  //double PhaseSpaceIn = 1 / (4*TMath::Abs(q1Pz*q2Pz-q1Pz*(-q2Pz)));

  double PhaseSpaceIn = 1 / (TMath::Abs(q1Pz/4000.0) *  TMath::Abs(q2Pz/4000.0) * pow(8000.0,2));
  double PhaseSpaceOut = 1 / (pow(2*TMath::Pi(),8) * (2*p1.E()*2*nu1.E()*2*p2.E()*2*nu2.E()));

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

  double matrix_element = matrix_elements[0] * pdf1 * pdf2 * PhaseSpaceIn * PhaseSpaceOut;
  return matrix_element;
  }


  //ClassDef(TFDISTR,1) //Class of testing functions for FOAM

//ClassImp(TFDISTR)

//int main(int argc, char** argv){
/*
Double_t ME(TLorentzVector ep, TLorentzVector mum, TLorentzVector Met){
};
//ClassImp(TFDISTR)
*/
//int main(int argc, char** argv){


//Double_t ComputePdf(char* argv[]){

  Double_t ComputePdf(int pid, double x, double q2){
  /*
    if (argc < 3) {
      cerr << "You must specify a PDF set and member number" << endl;
      return 1;
      }
  */
    //const string setname = "cteq6l1";//argv[1];
    //const string smem = 10042; //argv[2];
    //const int imem = boost::lexical_cast<int>(smem);
    //const PDF* pdf = mkPDF("cteq6l1", 0);

    const double xf = pdf->xfxQ2(pid, x, q2);


    //std::cout << "pdf : " << xf << std::endl;
    return xf;
  }
};


Double_t ME(TLorentzVector ep, TLorentzVector mum, TLorentzVector Met){

  //process.initProc("/home/fynu/amertens/scratch/MadGraph/madgraph5/w_cpp/Cards/param_card.dat");

  //TH2D *hst = new TH2D("test","test",150,0,150,100,-50,0);
  //TH2D *hst_Enu = new TH2D("test_2D", "test_2D", 150,0,150,150,0,150);
  TH1D *hst_Wm = new TH1D("test_1D", "test_1D", 150,0,150);


  TStopwatch chrono;
  std::cout << "initialising : " ;
  //chrono->Start();
  // TFoam Implementation
  Double_t *MCvect =new Double_t[4];
  TRandom3  *PseRan   = new TRandom3();
  PseRan->SetSeed(4357);  
  TFoam   *FoamX    = new TFoam("FoamX");
  FoamX->SetkDim(4);
  FoamX->SetnCells(500);      // No. of cells, can be omitted, default=2000
//  FoamX->SetRhoInt(Camel2);   // Set 2-dim distribution, included below

  TString pdfname = "cteq6l1";
  Int_t imem = 0;

  TFoamIntegrand *rho= new TFDISTR("/home/fynu/amertens/scratch/MadGraph/madgraph5/ww_cpp/Cards/param_card.dat", ep, mum, Met );


  //rho->SetParamCard("/home/fynu/amertens/scratch/MadGraph/madgraph5/w_cpp/Cards/param_card.dat");
  FoamX->SetRho(rho);
  FoamX->SetPseRan(PseRan);   // Set random number generator

  std::cout << "Starting initialisation..." << std::endl;
  FoamX->Initialize(); 
  
  std::cout << "CPU time : " << chrono.CpuTime() << "  Real-time : " << chrono.RealTime() << std::endl;
  chrono.Reset();
  chrono.Start();
//  std::cout << std::endl;
  std::cout << "Looping over phase-space points" << std::endl;

//  chrono->Reset();
//  chrono->Start();

  for(Long_t loop=0; loop<10000; loop++){
    FoamX->MakeEvent();          // generate MC event
    FoamX->GetMCvect( MCvect);   // get generated vector (x,y)
    Double_t px1=-150+300*MCvect[0];
    Double_t py1=-150+300*MCvect[1];
    Double_t pz1=-150+300*MCvect[2];
    Double_t pz2=-150+300*MCvect[3];
    //if(loop<10) cout<<"(x,y) =  ( "<< x <<", "<< y <<" )"<<endl;
    TLorentzVector nu1,nu2;
    nu1.SetPxPyPzE(px1,py1,pz1,TMath::Sqrt(pow(px1,2)+pow(py1,2)+pow(pz1,2)));
    
    hst_Wm->Fill((nu1+ep).M());           // fill scattergram
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
   hst_Wm->Draw();
   c->Print("Enu.png");

  return mcResult;
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



int main(int argc, char *argv[])
//const char *inputFile,const char *outputFile
{

  TString inputFile(argv[1]);
  TString outputFile(argv[2]);

  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchGen = treeReader->UseBranch("Particle");
  // Loop over all events

/*  TString pdfname = "cteq6l1";
  Int_t imem = 0;

  TFDISTR::SetPDF(pdfname, imem);
*/

  ofstream fout("ww_weights.out");
  for(Int_t entry = 0; entry < 10; ++entry)//numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    TLorentzVector gen_ep, gen_mum, gen_Met;

    GenParticle *gen;

    for (Int_t i = 0; i < branchGen->GetEntries(); i++)
    {
      gen = (GenParticle *) branchGen->At(i);
      if (gen->Status == 3)
      {
        if (gen->PID == -11) gen_ep = gen->P4();
        else if (gen->PID == 13) gen_mum = gen->P4();
	else if (gen->PID == 12) gen_Met += gen->P4();
	else if (gen->PID == -14) gen_Met += gen->P4();
      }
    }


  fout << entry << " " << ME(gen_ep, gen_mum, gen_Met) << endl;

  

  }
}
