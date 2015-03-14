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

#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "classes/DelphesClasses.h"

#include "src/HelAmps_sm.h"
#include "SubProcesses/P0_Sigma_sm_gg_epvebmumvmxbx/CPPProcess.h"
#include "src/rambo.h"
#include "src/Parameters_sm.h"

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

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"


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

  Double_t Px1=-500+1000*Xarg[0];
  Double_t Py1=-500+1000*Xarg[1];
  Double_t Pz1=-500+1000*Xarg[2];
  Double_t Px2=Met.Px()-Px1;
  Double_t Py2=Met.Py()-Py1;
  Double_t Pz2=-500+1000*Xarg[3];


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

  //cout << "q1Pz=" << q1Pz << ", q2Pz=" << q2Pz << endl;

  // momentum vector definition
  vector<double*> p(1, new double[4]);
  p[0][0] = q1Pz; p[0][1] = 0.0; p[0][2] = 0.0; p[0][3] = q1Pz;
  p.push_back(new double[4]);
  p[1][0] = TMath::Abs(q2Pz); p[1][1] = 0.0; p[1][2] = 0.0; p[1][3] = q2Pz;
  p.push_back(new double[4]);
  p[2][0] = p1.E(); p[2][1] = p1.Px(); p[2][2] = p1.Py(); p[2][3] = p1.Pz();
  p.push_back(new double[4]);
  p[3][0] = nu1.E(); p[3][1] = nu1.Px(); p[3][2] = nu1.Py(); p[3][3] = nu1.Pz();
  p.push_back(new double[4]);
  p[4][0] = p2.E(); p[4][1] = p2.Px(); p[4][2] = p2.Py(); p[4][3] = p2.Pz();
  p.push_back(new double[4]);
  p[5][0] = nu2.E(); p[5][1] = nu2.Px(); p[5][2] = nu2.Py(); p[5][3] = nu2.Pz();

  // Compute the Pdfs
  double pdf1_1 = ComputePdf( 2,TMath::Abs(q1Pz/4000.0), pow(Eext,2)) / TMath::Abs(q1Pz/4000.0);
  double pdf1_2 = ComputePdf(-2,TMath::Abs(q2Pz/4000.0), pow(Eext,2)) / TMath::Abs(q2Pz/4000.0);

  //double pdf2_1 = ComputePdf(-2,TMath::Abs(q1Pz/4000.0), pow(Eext,2)) / TMath::Abs(q1Pz/4000.0);
  //double pdf2_2 = ComputePdf( 2,TMath::Abs(q2Pz/4000.0), pow(Eext,2)) / TMath::Abs(q2Pz/4000.0);

  // Compute de Phase Space
  double PhaseSpaceIn = 1.0 / ( TMath::Abs(q1Pz/4000.0) *  TMath::Abs(q2Pz/4000.0) * pow(8000.0,2)); // ? factor 2?

  double dphip1 = p1.Pt()/(2.0*pow(2.0*TMath::Pi(),3));
  double dphip2 = p2.Pt()/(2.0*pow(2.0*TMath::Pi(),3));
  double dphinu1 = 1/(pow(2*TMath::Pi(),3)*2*nu1.E()); // nu1.Pt()/(2.0*pow(2.0*TMath::Pi(),3));
  double dphinu2 = 1/(pow(2*TMath::Pi(),3)*2*nu2.E()); // nu2.Pt()/(2.0*pow(2.0*TMath::Pi(),3));

  double PhaseSpaceOut = pow(2.0*TMath::Pi(),4) * 4./pow(8000.0,2) * dphip1 * dphip2 * dphinu1 * dphinu2;

  // Additional factor due to the integration range:
  double jac = pow(1000,4);

  //cout << "phase space=" << jac * PhaseSpaceIn * PhaseSpaceOut << ", pdfprod=" << pdf1_1*pdf1_2 << "\n\n";

  // Set momenta for this event
  process.setMomenta(p);

  // Evaluate matrix element
  process.sigmaKin();
  const Double_t* matrix_elements1 = process.getMatrixElements();
  
  // final integrand
  double matrix_element = jac * (matrix_elements1[0] * pdf1_1 * pdf1_2 /*+ matrix_elements1[0] * pdf2_1 * pdf2_2*/) * PhaseSpaceIn * PhaseSpaceOut;

  //cout << "|M|^2 = " << matrix_element << "     W :" << (p1+nu1).M() << " " << (p2+nu2).M() << endl;
  return matrix_element;
  }


  Double_t ComputePdf(int pid, double x, double q2){
    // return the xf(pid,x,q2), be careful: this value must be divided by x to obtain f
    const double xf = pdf->xfxQ2(pid, x, q2);
    return xf;
  }
};


Double_t ME(TLorentzVector ep, TLorentzVector mum, TLorentzVector Met){

  //TH1D *hst_Wm = new TH1D("test_1D", "test_1D", 150,0,150);


  TStopwatch chrono;
  std::cout << "initialising : " ;
  // TFoam Implementation
  Double_t *MCvect =new Double_t[4];
  TRandom3  *PseRan   = new TRandom3();
  PseRan->SetSeed(2245);  
  TFoam   *FoamX    = new TFoam("FoamX");
  FoamX->SetkDim(4);
  FoamX->SetnCells(4000);      // No. of cells, can be omitted, default=2000

  TString pdfname = "cteq6l1";
  //Int_t imem = 0;

  TFoamIntegrand *rho= new TFDISTR("/home/fynu/amertens/scratch/MadGraph/madgraph5/ww_cpp/Cards/param_card.dat", ep, mum, Met );

  FoamX->SetRho(rho);
  FoamX->SetPseRan(PseRan);   // Set random number generator

  std::cout << "Starting initialisation..." << std::endl;
  FoamX->Initialize(); 
  
  std::cout << "CPU time : " << chrono.CpuTime() << "  Real-time : " << chrono.RealTime() << std::endl;
  chrono.Reset();
  chrono.Start();
  std::cout << "Looping over phase-space points" << std::endl;

//  chrono->Reset();
//  chrono->Start();

  for(Long_t loop=0; loop<10000; loop++){
    FoamX->MakeEvent();          // generate MC event
    FoamX->GetMCvect( MCvect);   // get generated vector (x,y)
    /*
    Double_t px1=-500+1000*MCvect[0];
    Double_t py1=-500+1000*MCvect[1];
    Double_t pz1=-500+1000*MCvect[2];
    Double_t pz2=-500+1000*MCvect[3];
    //if(loop<10) cout<<"(x,y) =  ( "<< x <<", "<< y <<" )"<<endl;
    TLorentzVector nu1,nu2;
    nu1.SetPxPyPzE(px1,py1,pz1,TMath::Sqrt(pow(px1,2)+pow(py1,2)+pow(pz1,2)));
    //cout << "W mass : " << (nu1+ep).M() << endl;  
    hst_Wm->Fill((nu1+ep).M());           // fill scattergram
    */
   }// loop



  Double_t mcResult, mcError;
  FoamX->GetIntegMC( mcResult, mcError);  // get MC integral, should be one
  
  std::cout << "CPU time : " << chrono.CpuTime() << "  Real-time : " << chrono.RealTime() << std::endl;

  cout << " mcResult= " << mcResult << " +- " << mcError <<endl;
  // now hst_xy will be plotted visualizing generated distribution
  /*TCanvas *c = new TCanvas("c","Canvas for plotting",600,600);
   c->cd();
   hst_Wm->Draw();
   c->Print("Enu.png");*/

  return mcResult;
}


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
  //Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchGen = treeReader->UseBranch("Particle");
  // Loop over all events

/*  TString pdfname = "cteq6l1";
  Int_t imem = 0;

  TFDISTR::SetPDF(pdfname, imem);
*/

  double weight = 0;

  ofstream fout(outputFile);

  double madweight[10] = {5.70122210379e-16, 1.53206147539e-18, 9.33916724451e-17, 6.72210503551e-18, 2.27136432783e-16, 5.05485330753e-17, 4.05720785858e-16, 1.51874653102e-17, 6.68618012851e-17, 1.67333412761e-17};

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

  weight = ME(gen_ep, gen_mum, gen_Met);
  fout << entry << " " << weight << "   " << weight / (double) madweight[entry] << endl;

  

  }
}
