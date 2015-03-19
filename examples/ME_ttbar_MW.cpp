#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <complex>
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
#include "TLorentzVector.h"
#include "TMath.h"

#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TRandom3.h"

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"

#define SSTR( x ) dynamic_cast< std::ostringstream & > \
        ( std::ostringstream() << std::dec << x ).str()

using namespace LHAPDF;
using namespace std::literals;
using namespace std;

typedef complex<double> cd;

bool solveQuartic(double a, double b, double c, double d, double e, vector<double>& roots){
	// see http://en.wikipedia.org/wiki/Quartic_function#Solving_a_quartic_equation
	
	if(a == 0.)
		return false;
	
	/*double delta = 256.*pow(a,3.)*pow(e,3.) - 192.*pow(a,2.)*b*d*pow(e,2.) - 128.*pow(a*c*e,2.)
		+ 144.*c*e*pow(a*d,2.) - 27.*pow(a,2.)*pow(d,4.) + 144.*a*c*pow(b*e,2.) - 6.*a*e*pow(b*d,2.)
		- 80.*a*b*d*e*pow(c,2.) + 18.*a*b*c*pow(d,3.) + 16*a*e*pow(c,4.) - 4.a*pow(c,3.)*pow(d,2.)
		- 27.*pow(b,4.)*pow(e,2.) + 18*c*d*e*pow(b,3.) - 4.*pow(d*b,3.) - 4.*e*pow(b,2.)*pow(c,3.)
		+ pow(b*c*d,2.);*/
	
	cd delta0 = pow(c,2.) - 3.*b*d + 12.*a*e;
	cd delta1 = 2.*pow(c,3.) - 9.*b*c*d + 27.*pow(b,2)*e + 27.*a*pow(d,2.) - 72.*a*c*e;
	cd p = (8.*a*c - 3.*pow(b,2.))/(8.*pow(a,2.));
	cd q = (pow(b,3.) - 4.*a*b*c + 8.*pow(a,2.)*d)/(8*pow(a,3.));

	cd Q = pow( ( delta1 + sqrt( pow(delta1,2.) - 4.*pow(delta0,3.) ) )/2., 1./3.);
	cd S = 0.5 * sqrt( -2./3. * p + 1./(3.*a) * ( Q + delta0/Q ) );

	// four solutions
	
	cd x1 = -b/(4.*a) - S + 0.5 * sqrt( -4.*pow(S,2.) - 2.*p + q/S );
	cd x2 = -b/(4.*a) - S - 0.5 * sqrt( -4.*pow(S,2.) - 2.*p + q/S );

	cd x3 = -b/(4.*a) + S + 0.5 * sqrt( -4.*pow(S,2.) - 2.*p - q/S );
	cd x4 = -b/(4.*a) + S - 0.5 * sqrt( -4.*pow(S,2.) - 2.*p - q/S );

	// checking which ones are real
	
	if( abs(imag(x1)) < pow(10,-15) )
		roots.push_back(real(x1));
	
	if( abs(imag(x2)) < pow(10,-15) )
		roots.push_back(real(x2));
	
	if( abs(imag(x3)) < pow(10,-15) )
		roots.push_back(real(x3));
	
	if( abs(imag(x4)) < pow(10,-15) )
		roots.push_back(real(x4));

	return true;
}

bool solve2Quads(double a11, double a22, double a12, double a10, double a01, double a00,
				double b11, double b22, double b12, double b10, double b01, double b00,
				vector<double>& E1, vector<double>& E2){
	double alpha = b11*a22-a11*b22;
	double beta = b11*a12-a11*b12;
	double gamma = b11*a10-a11*b10;
	double delta = b11*a01-a11*b01;
	double omega = b11*a00-a11*b00;

	double a = a11*pow(alpha,2.) + a22*pow(delta,2.) - a12*alpha*beta;
	double b = 2.*a11*alpha*delta - a12*alpha*gamma - a12*delta*beta - a10*alpha*beta + 2.*a22*beta*gamma;
	double c = a11*pow(delta,2.) + 2.*a11*alpha*omega - a12*delta*gamma - a12*omega*beta - a10*alpha*gamma - a10*delta*beta
		+ a22*pow(gamma,2.) + 2.*a01*beta*gamma + a00*pow(beta,2.);
	double d = 2.*a11*delta*omega - a12*omega*gamma - a10*delta*gamma - a10*omega*beta + a01*pow(gamma,2.) + 2.*a00*beta*gamma;
	double e = a11*pow(omega,2.) - a10*omega*gamma + a00*pow(gamma,2.);

	solveQuartic(a, b, c, d, e, E2);

	for(unsigned int i =0; i < E2.size(); ++i){
		double e2 = E2.at(i);
		double e1 = -(alpha * pow(e2,2.) + delta*e2 + omega)/(beta*e2 + gamma);
		E1.push_back(e1);
	}

	return true;
}

int mycount = 0, count_wgt = 0, count_perm=1;

class TFDISTR: public TFoamIntegrand {
private:
  CPPProcess process;
  TLorentzVector p3, p4, p5, p6;
  TLorentzVector Met;

  const PDF *pdf;

public:
  TFDISTR(std::string paramCardPath, TLorentzVector ep, TLorentzVector mum, TLorentzVector b, TLorentzVector bbar, TLorentzVector met){
  process.initProc(paramCardPath);
  p3 = ep;
  p5 = mum;
  p4 = b;
  p6 = bbar;
  Met = met;
  pdf=LHAPDF::mkPDF("cteq6l1", 0);
  
  Double_t Density(int nDim, Double_t *Xarg){
  // Integrand for mFOAM
  
  double range1 = TMath::Pi()/2. + TMath::Arctan(M_W/G_W);
  double y1 = - TMath::Arctan(M_W/G_W) + range1 * Xarg[0];
  double s13 = M_W * G_W * TMath::Tan(y1) + pow(M_W,2.);


  double range2 = TMath::Pi()/2. + TMath::Arctan(M_T/G_T);
  double y2 = - TMath::Arctan(M_T/G_T) + range2 * Xarg[1];
  double s134 = M_T * G_T * TMath::Tan(y2) + pow(M_T,2.);


  double range3 = TMath::Pi()/2. + TMath::Arctan(M_W/G_W);
  double y3 = - TMath::Arctan(M_W/G_W) + range3 * Xarg[2];
  double s25 = M_W * G_W * TMath::Tan(y3) + pow(M_W,2.);


  double range4 = TMath::Pi()/2. + TMath::Arctan(M_T/G_T);
  double y4 = - TMath::Arctan(M_T/G_T) + range4 * Xarg[3];
  double s256 = M_T * G_T * TMath::Tan(y4) + pow(M_T,2.);


  TLorentzVector pT = p3 + p4 + p5 + p6;
  pT.SetPz(0.);

  double p34 = p3*p4;
  double p56 = p5*p6;
  double p33 = p5.M2();
  double p44 = p4.M2();
  double p55 = p5.M2();
  double p66 = p6.M2();

  double A1 = 2.*( -p3.Px() + p3.Pz()*p4.Px()/p4.Pz() );
  double A2 = 2.*( p5.Px() - p5.Pz()*p6.Px()/p6.Pz() );
  
  double B1 = 2.*( -p3.Py() + p3.Pz()*p4.Py()/p4.Pz() );
  double B2 = 2.*( p5.Py() - p5.Pz()*p6.Py()/p6.Pz() );

  double Dx = B2*A1 - B1*A2;
  double Dy = A2*B1 - A1*B2;

  double X = 2*( pT.Px()*p5.Px() + pT.Py()*p5.Py() - p5.Pz()/p6.Pz()*( 0.5*(s25 - s256 + p56 + p66) + pT.Px()*p6.Px() + pT.Py()*p6.Py() ) ) + p55 - s25;
  double Y = p3.Pz()/p4.Pz()*( s13 - s134 + p34 + p44 ) - p33 + s13;
  
  double alpha1 = -2*B2*(p3.E() - p4.E()*p3.Pz()/p4.Pz())/Dx;
  double beta1 = 2*B1*(p5.E() - p6.E()*p5.Pz()/p6.Pz())/Dx;
  double gamma1 = B1*X/Dx + B2*Y/Dx;

  double alpha2 = -2*A2*(p3.E() - p4.E()*p3.Pz()/p4.Pz())/Dy;
  double beta2 = 2*A1*(p5.E() - p6.E()*p5.Pz()/p6.Pz())/Dy;
  double gamma2 = A1*X/Dy + A2*Y/Dy;

  double alpha3 = (p4.E() - alpha1*p4.Px() - alpha2*p4.Py())/p4.Pz();
  double beta3 = -(beta1*p4.Px() + beta2*p4.Py())/p4.Pz();
  double gamma3 = ( 0.5*(s13 - s134 + p34 + p44) - gamma1*p4.Px() - gamma2*p4.Py() )/p4.Pz();

  double alpha4 = (alpha1*p6.Px() + alpha2*p6.Py())/p6.Pz();
  double beta4 = (p6.E() + beta1*p6.Px() + beta2*p6.Py())/p6.Pz();
  double gamma4 = ( 0.5*(s25 - s256 + p56 + p66) - (gamma1 + pT.Px())*p6.Px() - (gamma2 + pT.Py())*p6.Py() )/p6.Pz();

  double alpha5 = -alpha1;
  double beta5 = -beta1;
  double gamma5 = -pT.Px() - gamma1;

  double alpha6 = -alpha2;
  double beta6 = -beta2;
  double gamma6 = -pT.Py() - gamma2;

  double a11 = 1 - ( pow(alpha1,2.) + pow(alpha2,2.) + pow(alpha3,2.) );
  double a22 = - ( pow(beta1,2.) + pow(beta2,2.) + pow(beta3,2.) );
  double a12 = - 2.*( alpha1*beta1 + alpha2*beta2 + alpha3*beta3 );
  double a10 = - 2.*( alpha1*gamma1 + alpha2*gamma2 + alpha3*gamma3 );
  double a01 = - 2.*( beta1*gamma1 + beta2*gamma2 + beta3*gamma3 );
  double a00 = - ( pow(gamma1,2.) + pow(gamma2,2.) + pow(gamma3,2.) );

  double b11 = 1 - ( pow(alpha5,2.) + pow(alpha6,2.) + pow(alpha4,2.) );
  double b22 = - ( pow(beta5,2.) + pow(beta6,2.) + pow(beta4,2.) );
  double b12 = - 2.*( alpha5*beta5 + alpha6*beta6 + alpha4*beta4 );
  double b10 = - 2.*( alpha5*gamma5 + alpha6*gamma6 + alpha4*gamma4 );
  double b01 = - 2.*( beta5*gamma5 + beta6*gamma6 + beta4*gamma4 );
  double b00 = - ( pow(gamma5,2.) + pow(gamma6,2.) + pow(gamma4,2.) );

  vector<double> E1, E2;
  solve2Quads(a11, a22, a12, a10, a01, a00, b11, b22, b12, b10, b01, b00, E1, E2);

  TLorentzVector nu1,nu2;
  nu1.SetPxPyPzE(Px1,Py1,Pz1,E1);
  nu2.SetPxPyPzE(Px2,Py2,Pz2,E2);

  /*cout << "Electron" << endl;
  cout << pep.E() << "," << pep.Px() << "," << pep.Py() << "," << pep.Pz() << endl;
  cout << "Electron neutrino" << endl;
  cout << nu1.E() << "," << nu1.Px() << "," << nu1.Py() << "," << nu1.Pz() << endl;
  cout << "b quark" << endl;
  cout << pb.E() << "," << pb.Px() << "," << pb.Py() << "," << pb.Pz() << endl;
  cout << "Muon" << endl;
  cout << pmum.E() << "," << pmum.Px() << "," << pmum.Py() << "," << pmum.Pz() << endl;
  cout << "Muon neutrino" << endl;
  cout << nu2.E() << "," << nu2.Px() << "," << nu2.Py() << "," << nu2.Pz() << endl;
  cout << "Anti b quark" << endl;
  cout << pbbar.E() << "," << pbbar.Px() << "," << pbbar.Py() << "," << pbbar.Pz() << endl;*/

  TLorentzVector tot = nu1 + nu2 + pep + pmum + pb + pbbar;

 // Double_t gEcm=(EEl+ENu)/2.0;

  Double_t Eext=tot.E();
  Double_t Pzext=tot.Pz();

  Double_t q1Pz=(Pzext+Eext)/2;
  Double_t q2Pz=(Pzext-Eext)/2;

  //cout << " ===> Eext=" << Eext << ", Pzext=" << Pzext << ", q1Pz=" << q1Pz << ", q2Pz=" << q2Pz << endl;
  //cout << "q1Pz=" << q1Pz << ", q2Pz=" << q2Pz << endl;

  if(q1Pz > 6500. || q2Pz < -6500. || q1Pz < 0. || q2Pz > 0.){
	//cout << "Fail!" << endl;
	mycount++;
	return 0.;
  }


  // momentum vector definition
  vector<double*> p(1, new double[4]);
  p[0][0] = q1Pz; p[0][1] = 0.0; p[0][2] = 0.0; p[0][3] = q1Pz;
  p.push_back(new double[4]);
  p[1][0] = TMath::Abs(q2Pz); p[1][1] = 0.0; p[1][2] = 0.0; p[1][3] = q2Pz;
  p.push_back(new double[4]);
  p[2][0] = pep.E(); p[2][1] = pep.Px(); p[2][2] = pep.Py(); p[2][3] = pep.Pz();
  p.push_back(new double[4]);
  p[3][0] = nu1.E(); p[3][1] = nu1.Px(); p[3][2] = nu1.Py(); p[3][3] = nu1.Pz();
  p.push_back(new double[4]);
  p[4][0] = pb.E(); p[4][1] = pb.Px(); p[4][2] = pb.Py(); p[4][3] = pb.Pz();
  p.push_back(new double[4]);
  p[5][0] = pmum.E(); p[5][1] = pmum.Px(); p[5][2] = pmum.Py(); p[5][3] = pmum.Pz();
  p.push_back(new double[4]);
  p[6][0] = nu2.E(); p[6][1] = nu2.Px(); p[6][2] = nu2.Py(); p[6][3] = nu2.Pz();
  p.push_back(new double[4]);
  p[7][0] = pbbar.E(); p[7][1] = pbbar.Px(); p[7][2] = pbbar.Py(); p[7][3] = pbbar.Pz();

  // Compute the Pdfs
  double pdf1_1 = ComputePdf(21,TMath::Abs(q1Pz/6500.0), pow(tot.M(),2)) / TMath::Abs(q1Pz/6500.0);
  double pdf1_2 = ComputePdf(21,TMath::Abs(q2Pz/6500.0), pow(tot.M(),2)) / TMath::Abs(q2Pz/6500.0);

  // Compute de Phase Space
  double PhaseSpaceIn = 1.0 / ( TMath::Abs(q1Pz/6500.0) *  TMath::Abs(q2Pz/6500.0) * pow(13000.0,2)); // ? factor 2?

  double dphipep = pep.Pt()/(2.0*pow(2.0*TMath::Pi(),3));
  double dphipmum = pmum.Pt()/(2.0*pow(2.0*TMath::Pi(),3));
  double dphipb = pow(pb.P(),2.)*TMath::Sin(pb.Theta())/(2.0*pb.E()*pow(2.0*TMath::Pi(),3)); // massive b's!
  double dphipbbar = pow(pbbar.P(),2.)*TMath::Sin(pbbar.Theta())/(2.0*pbbar.E()*pow(2.0*TMath::Pi(),3));
  double dphinu1 = 1./(pow(2*TMath::Pi(),3)*2*nu1.E()); // nu1.Pt()/(2.0*pow(2.0*TMath::Pi(),3));
  double dphinu2 = 1./(pow(2*TMath::Pi(),3)*2*nu2.E()); // nu2.Pt()/(2.0*pow(2.0*TMath::Pi(),3));

  double PhaseSpaceOut = pow(2.0*TMath::Pi(),4) * 4./pow(13000.0,2) * dphipep * dphipmum * dphipb * dphipbbar * dphinu1 * dphinu2;

  // Additional factor due to the integration range:
  double jac = pow(TMath::Pi(),4.)/( pow(TMath::Cos(-TMath::Pi()/2. + TMath::Pi()*Xarg[0]), 2.) * pow(TMath::Cos(-TMath::Pi()/2. + TMath::Pi()*Xarg[1]), 2.)
	* pow(TMath::Cos(-TMath::Pi()/2. + TMath::Pi()*Xarg[2]), 2.) * pow(TMath::Cos(-TMath::Pi()/2. + TMath::Pi()*Xarg[3]), 2.) );

  //cout << "phase space=" << jac * PhaseSpaceIn * PhaseSpaceOut << ", pdfprod=" << pdf1_1*pdf1_2 << "\n\n";

  // Set momenta for this event
  process.setMomenta(p);

  // Evaluate matrix element
  process.sigmaKin();
  const Double_t* matrix_elements1 = process.getMatrixElements();
  
  // final integrand
  double matrix_element = jac * (matrix_elements1[0] * pdf1_1 * pdf1_2) * PhaseSpaceIn * PhaseSpaceOut;

  //cout << "|M|^2 = " << matrix_element << "     W :" << (pep+nu1).M() << " " << (pmum+nu2).M() << endl;
  return matrix_element;
  }


  Double_t ComputePdf(int pid, double x, double q2){
    // return the xf(pid,x,q2), be careful: this value must be divided by x to obtain f
    const double xf = pdf->xfxQ2(pid, x, q2);
    return xf;
  }
};


Double_t ME(TLorentzVector ep, TLorentzVector mum, TLorentzVector b, TLorentzVector bbar, TLorentzVector Met){

  TH1D *hst_Wm = new TH1D("test_mu", "test_1D", 150,0,150);
  TH1D *hst_We = new TH1D("test_ep", "test_1D", 150,0,150);
  TH1D *hst_t = new TH1D("test_t", "test_1D", 100,150,250);
  TH1D *hst_tbar = new TH1D("test_tbar", "test_1D", 100,150,250);


  TStopwatch chrono;
  std::cout << "initialising : " ;
  // TFoam Implementation
  Double_t *MCvect = new Double_t[4];
  TRandom3  *PseRan   = new TRandom3();
  PseRan->SetSeed(2245);  
  TFoam   *FoamX    = new TFoam("FoamX");
  FoamX->SetkDim(4);
  FoamX->SetnCells(20000);      // No. of cells, can be omitted, default=2000
  FoamX->SetnSampl(10);

  TString pdfname = "cteq6l1";
  //Int_t imem = 0;

  TFoamIntegrand *rho= new TFDISTR("/home/fynu/swertz/scratch/Madgraph/madgraph5/cpp_ttbar_epmum/Cards/param_card.dat", ep, mum, b, bbar, Met );

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

  for(Long_t loop=0; loop<50000; loop++){
    int count_old = mycount;
    FoamX->MakeEvent();          // generate MC event
    FoamX->GetMCvect( MCvect );   // get generated vector (x,y)
    
    /*Double_t px1=-500+1000*MCvect[0];
    Double_t py1=-500+1000*MCvect[1];
    Double_t pz1=-500+1000*MCvect[2];
    Double_t pz2=-500+1000*MCvect[3];*/
	Double_t px1=TMath::Tan(-TMath::Pi()/2. + TMath::Pi()*MCvect[0]);
  	Double_t py1=TMath::Tan(-TMath::Pi()/2. + TMath::Pi()*MCvect[1]);
  	Double_t pz1=TMath::Tan(-TMath::Pi()/2. + TMath::Pi()*MCvect[2]);
  	Double_t pz2=TMath::Tan(-TMath::Pi()/2. + TMath::Pi()*MCvect[3]);
    //if(loop<10) cout<<"(x,y) =  ( "<< x <<", "<< y <<" )"<<endl;
    TLorentzVector nu1,nu2;
    nu1.SetPxPyPzE(px1,py1,pz1,TMath::Sqrt(pow(px1,2)+pow(py1,2)+pow(pz1,2)));
    nu2.SetPxPyPzE(Met.Px()-px1,Met.Py()-py1,pz2,TMath::Sqrt(pow(Met.Px()-px1,2)+pow(Met.Py()-py1,2)+pow(pz2,2)));
    //cout << "W mass : " << (nu1+ep).M() << endl;  
	if(count_old == mycount){
		hst_We->Fill((nu1+ep).M());           // fill scattergram
		hst_Wm->Fill((nu2+mum).M());           // fill scattergram
		hst_t->Fill((nu1+ep+b).M());
		hst_tbar->Fill((nu2+mum+bbar).M());
	}
    
   }// loop



  Double_t mcResult, mcError;
  FoamX->GetIntegMC( mcResult, mcError);  // get MC integral, should be one
  
  std::cout << "CPU time : " << chrono.CpuTime() << "  Real-time : " << chrono.RealTime() << std::endl;
  cout << "nr. fails: " << mycount  << endl;

  cout << " mcResult= " << mcResult << " +- " << mcError <<endl;
  // now hst_xy will be plotted visualizing generated distribution
  TCanvas *c = new TCanvas("c","Canvas for plotting",600,600);
   c->cd();
   hst_We->Draw();
   c->Print(TString("plots/")+SSTR(count_wgt)+"_"+SSTR(count_perm)+"_Enu.png");
   delete hst_We;
   delete c;
  
   c = new TCanvas("c","Canvas for plotting",600,600);
   c->cd();
   hst_Wm->Draw();
   c->Print(TString("plots/")+SSTR(count_wgt)+"_"+SSTR(count_perm)+"_Munu.png");
   delete hst_Wm;
   delete c;


   c = new TCanvas("c","Canvas for plotting",600,600);
   c->cd();
   hst_t->Draw();
   c->Print(TString("plots/")+SSTR(count_wgt)+"_"+SSTR(count_perm)+"_t.png");
   delete hst_t;
   delete c;
  
   c = new TCanvas("c","Canvas for plotting",600,600);
   c->cd();
   hst_tbar->Draw();
   c->Print(TString("plots/")+SSTR(count_wgt)+"_"+SSTR(count_perm)+"_tbar.png");
   delete hst_tbar;
   delete c;

  delete FoamX;
  delete PseRan;
  delete MCvect;

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

  cout << "Entries:" << chain.GetEntries() << endl;

  double weight1 = 0, weight2 = 0;

  ofstream fout(outputFile);

  double madweight1[10] = {1.58495292058e-21, 2.09681384879e-21, 4.34399623629e-22, 1.68163897955e-22, 3.20350498956e-22, 5.22232034307e-22, 6.04738375743e-21, 9.55643564854e-22, 8.12425265344e-22, 5.81210532053e-23};
  double madweight2[10] = {1.02514966131e-21, 1.45375719248e-21, 1.65080839221e-22, 1.55653414654e-24, 5.60531044001e-25, 1.,  9.70526105314e-24, 3.89103636371e-22, 6.38206925825e-23, 9.37189585544e-26};

  for(Int_t entry = 0; entry < 10; ++entry)//numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    TLorentzVector gen_ep, gen_mum, gen_b, gen_bbar, gen_Met;

    GenParticle *gen;

	//int count_ep=0, count_mum=0;

    for (Int_t i = 0; i < branchGen->GetEntries(); i++)
    {
      gen = (GenParticle *) branchGen->At(i);
	  //cout << "Status=" << gen->Status << ", PID=" << gen->PID << ", E=" << gen->P4().E() << endl;
      if (gen->Status == 1)
      {
        if (gen->PID == -11){
			gen_ep = gen->P4();
			//count_ep++;
		}
        else if (gen->PID == 13){
			gen_mum = gen->P4();
			//count_mum++;
		}
		else if (gen->PID == 12) gen_Met += gen->P4();
		else if (gen->PID == -14) gen_Met += gen->P4();
		else if (gen->PID == 5) gen_b = gen->P4();
		else if (gen->PID == -5) gen_bbar = gen->P4();
      }
    }

	//if(count_ep != 1 || count_mum != 1)
	//	continue;
	//gen_Met.SetPz(0.);
  
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

  count_perm = 1;
  weight1 = ME(gen_ep, gen_mum, gen_b, gen_bbar, gen_Met);
  count_perm = 2;
  weight2 = ME(gen_ep, gen_mum, gen_bbar, gen_b, gen_Met);

  double weight3;
  if(weight1 < weight2){
	weight3 = weight2;
	weight2 = weight1;
	weight1 = weight3;
  }

  fout << entry << " " << weight1 << "   " << weight1 / (double) madweight1[entry] << ", " << weight2 << "    " << weight2/madweight2[entry] <<  endl;

  count_wgt++;

  }
 //delete treeReader; 
}
