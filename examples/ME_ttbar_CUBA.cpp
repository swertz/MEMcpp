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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TRandom3.h"

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"

#include "cuba.h"

#define M_T 173.
#define G_T 1.4915

#define M_W 80.419
#define G_W 2.0476

#define SQRT_S 13000

#define JAC_MIN 10e-20

#define SSTR( x ) dynamic_cast< std::ostringstream & > \
        ( std::ostringstream() << std::dec << x ).str()

using namespace LHAPDF;
using namespace std;

typedef complex<double> cd;

double solveQuadratic(double a, double b, double c, vector<double>& roots){
	if(a == 0){
		roots.push_back(-c/b);
		return true;
	}

	double rho = pow(b,2.) - 4*a*c;

	if(rho >= 0.){
		double x1 = (-b + sqrt(rho))/(2*a);
		double x2 = (-b - sqrt(rho))/(2*a);

		roots.push_back(x1);
		roots.push_back(x2);

		return true;
	}else
		return false;
}

double solveCubic(double a, double b, double c, double d, vector<double>& roots){
	if(a == 0)
		return solveQuadratic(b, c, d, roots);

	cd min1 = -1.;
	cd im = sqrt(min1);

	cd u1 = 1.;
	cd u2 = (-1. + im*sqrt(3.))/2.;
	cd u3 = (-1. - im*sqrt(3.))/2.;

	double delta0 = pow(b,2.) - 3*a*c;
	double delta1 = 2*pow(b,3.) - 9*a*b*c + 27*d*pow(a,2.);

	cd C = pow((delta1 + sqrt(pow(delta1,2.) - 4*pow(delta0,3.)))/2., 1./3.);

	cd x1 = -1./(3*a) * ( b + u1*C + delta0/(u1*C) );
	cd x2 = -1./(3*a) * ( b + u2*C + delta0/(u2*C) );
	cd x3 = -1./(3*a) * ( b + u3*C + delta0/(u3*C) );
	
	if( abs(real(x1)) > pow(10,5)*abs(imag(x1)) && abs(real(x1)) != 0 )
		roots.push_back(real(x1));
	
	if( abs(real(x2)) > pow(10,5)*abs(imag(x2)) && abs(real(x2)) != 0 )
		roots.push_back(real(x2));
	
	if( abs(real(x3)) > pow(10,5)*abs(imag(x3)) && abs(real(x3)) != 0 )
		roots.push_back(real(x3));

	if(roots.size())
		return true;
	else
		return false;
}

bool solveQuartic(double a, double b, double c, double d, double e, vector<double>& roots){
	// see http://en.wikipedia.org/wiki/Quartic_function#Solving_a_quartic_equation
	
	if(a == 0.)
		return solveCubic(b, c, d, e, roots);
	
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

	/*cout << "x1 = " << real(x1) << "+i*" << imag(x1) << endl;
	cout << "x2 = " << real(x2) << "+i*" << imag(x2) << endl;
	cout << "x2 = " << real(x3) << "+i*" << imag(x3) << endl;
	cout << "x4 = " << real(x4) << "+i*" << imag(x4) << endl;*/ 

	// checking which ones are real, with some tolerance
	
	if( abs(real(x1)) > pow(10,5)*abs(imag(x1)) && abs(real(x1)) != 0 )
		roots.push_back(real(x1));
	
	if( abs(real(x2)) > pow(10,5)*abs(imag(x2)) && abs(real(x2)) != 0 )
		roots.push_back(real(x2));
	
	if( abs(real(x3)) > pow(10,5)*abs(imag(x3)) && abs(real(x3)) != 0 )
		roots.push_back(real(x3));
	
	if( abs(real(x4)) > pow(10,5)*abs(imag(x4)) && abs(real(x4)) != 0 )
		roots.push_back(real(x4));

	if(roots.size())
		return true;
	else
		return false;
}

bool solve2Quads(double a11, double a22, double a12, double a10, double a01, double a00,
				double b11, double b22, double b12, double b10, double b01, double b00,
				vector<double>& E1, vector<double>& E2){
	double alpha = b11*a22-a11*b22;
	double beta = b11*a12-a11*b12;
	double gamma = b11*a10-a11*b10;
	double delta = b11*a01-a11*b01;
	double omega = b11*a00-a11*b00;

	double a = a11*pow(alpha,2.) + a22*pow(beta,2.) - a12*alpha*beta;
	double b = 2.*a11*alpha*delta - a12*alpha*gamma - a12*delta*beta - a10*alpha*beta + 2.*a22*beta*gamma + a01*pow(beta,2.);
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

double BreitWigner(double s, double m, double g){
  /*double ga = sqrt(m*m*(m*l+g*g));
  double k = 2*sqrt(2)*m*g*ga/(TMath::Pi()*sqrt(m*m+ga));*/
  double k = m*g;

  cout << "BW(" << s << "," << m << "," << g << ")=" << k/(pow(s-m*m,2.) + pow(m*g,2.)) << endl;

  return k/(pow(s-m*m,2.) + pow(m*g,2.));
}

double computeJacobian(vector<TLorentzVector> p){
	
	double E1  = p.at(0).E();
	double p1x = p.at(0).Px();
	double p1y = p.at(0).Py();
	double p1z = p.at(0).Pz();

	double E2  = p.at(1).E();
	double p2x = p.at(1).Px();
	double p2y = p.at(1).Py();
	double p2z = p.at(1).Pz();

	double E3  = p.at(2).E();
	double p3x = p.at(2).Px();
	double p3y = p.at(2).Py();
	double p3z = p.at(2).Pz();

	double E4  = p.at(3).E();
	double p4x = p.at(3).Px();
	double p4y = p.at(3).Py();
	double p4z = p.at(3).Pz();

	double E5  = p.at(4).E();
	double p5x = p.at(4).Px();
	double p5y = p.at(4).Py();
	double p5z = p.at(4).Pz();

	double E6  = p.at(5).E();
	double p6x = p.at(5).Px();
	double p6y = p.at(5).Py();
	double p6z = p.at(5).Pz();

	double E34  = E3 + E4;
	double p34x = p3x + p4x;
	double p34y = p3y + p4y;
	double p34z = p3z + p4z;

	double E56  = E5 + E6;
	double p56x = p5x + p6x;
	double p56y = p5y + p6y;
	double p56z = p5z + p6z;

    double jac = E3*(E5*
             (p34z*(p1y*p2z*p56x - p1x*p2z*p56y - p1y*p2x*p56z + 
                  p1x*p2y*p56z) + 
               p1z*(-(p2z*p34y*p56x) + p2z*p34x*p56y - 
                  p2y*p34x*p56z + p2x*p34y*p56z)) + 
            (E56*p2z - E2*p56z)*
             (p1z*p34y*p5x - p1y*p34z*p5x - p1z*p34x*p5y + 
               p1x*p34z*p5y) + 
            (E56*(p1z*p2y*p34x - p1z*p2x*p34y + p1y*p2x*p34z - 
                  p1x*p2y*p34z) + 
               E2*(p1z*p34y*p56x - p1y*p34z*p56x - p1z*p34x*p56y + 
                  p1x*p34z*p56y))*p5z) + 
         E34*(E5*p2z*(p1z*p3y*p56x - p1y*p3z*p56x - p1z*p3x*p56y + 
               p1x*p3z*p56y) + 
            E5*(p1z*p2y*p3x - p1z*p2x*p3y + p1y*p2x*p3z - 
               p1x*p2y*p3z)*p56z - 
            (E56*p2z - E2*p56z)*
             (p1z*p3y*p5x - p1y*p3z*p5x - p1z*p3x*p5y + p1x*p3z*p5y)
              - (E56*(p1z*p2y*p3x - p1z*p2x*p3y + p1y*p2x*p3z - 
                  p1x*p2y*p3z) + 
               E2*(p1z*p3y*p56x - p1y*p3z*p56x - p1z*p3x*p56y + 
                  p1x*p3z*p56y))*p5z) + 
         E1*(E5*(p2z*(-(p34z*p3y*p56x) + p34y*p3z*p56x + 
                  p34z*p3x*p56y - p34x*p3z*p56y) + 
               (-(p2y*p34z*p3x) + p2x*p34z*p3y + p2y*p34x*p3z - 
                  p2x*p34y*p3z)*p56z) + 
            (E56*p2z - E2*p56z)*
             (p34z*p3y*p5x - p34y*p3z*p5x - p34z*p3x*p5y + 
               p34x*p3z*p5y) + 
            (E56*(p2y*p34z*p3x - p2x*p34z*p3y - p2y*p34x*p3z + 
                  p2x*p34y*p3z) + 
               E2*(p34z*p3y*p56x - p34y*p3z*p56x - p34z*p3x*p56y + 
                  p34x*p3z*p56y))*p5z);

	//cout << "jac=" << abs(jac) << endl;
	
	if(abs(jac) < JAC_MIN)
		return -1.;
	else
		return 1./TMath::Abs(jac);
}

int mycount = 0, count_wgt = 0, count_perm=1;

class TFDISTR: public TFoamIntegrand {
private:
  CPPProcess process;
  TLorentzVector p3, p4, p5, p6;
  TLorentzVector Met;

  const PDF *pdf;

public:

  TFDISTR(std::string paramCardPath, TLorentzVector ep, TLorentzVector mum, TLorentzVector b, TLorentzVector bbar, TLorentzVector met);
  double Density(int nDim, double *Xarg);
  double ComputePdf(int pid, double x, double q2);
};

TFDISTR::TFDISTR(std::string paramCardPath, TLorentzVector ep, TLorentzVector mum, TLorentzVector b, TLorentzVector bbar, TLorentzVector met){
  process.initProc(paramCardPath);
  p3 = ep;
  p5 = mum;
  p4 = b;
  p4.SetE(p4.P());
  p6 = bbar;
  p6.SetE(p6.P());
  Met = met;
  pdf=LHAPDF::mkPDF("cteq6l1", 0);
}

double TFDISTR::Density(int nDim, double *Xarg){
  // Integrand for mFOAM
  //cout << endl << endl << endl << "########## Starting phase-space point ############" << endl << endl;

  //cout << "Inputs = [" << Xarg[0] << "," << Xarg[1] << "," << Xarg[2] << "," << Xarg[3] << endl;

  for(int i=0; i<nDim; i++){
	if(Xarg[i] == 1){
		mycount++;
		return 0.;
	}
  }

  // We flatten the Breit-Wigners by doing a change of variable for each resonance separately
  // s = M G tan(y) + M^2
  // jac = M G / cos^2(y)
  // ==> BW(s(y))*jac(y) is flat in the variable y, as BW(s) = 1/((s-M^2)^2 - (GM)^2)
  // Where y = -arctan(M/G) + (pi/2+arctan(M/G))*x_foam (x_foam between 0 and 1 => s between 0 and infinity)

  //double range1 = TMath::Pi();
  //double y1 = - TMath::Pi()/2. + range1 * Xarg[0];
  double range1 = TMath::Pi()/2. + TMath::ATan(M_W/G_W);
  double y1 = - TMath::ATan(M_W/G_W) + range1 * Xarg[0];
  double s13 = M_W * G_W * TMath::Tan(y1) + pow(M_W,2.);
  //double range1 = 2000;
  //double s13 = pow(M_W,2.)-1000 + range1*Xarg[0];

  //cout << "y1=" << y1 << ", m13=" << TMath::Sqrt(s13) << endl;

  //double range2 = TMath::Pi();
  //double y2 = - TMath::Pi()/2. + range2 * Xarg[1];
  double range2 = TMath::Pi()/2. + TMath::ATan(M_T/G_T);
  double y2 = - TMath::ATan(M_T/G_T) + range2 * Xarg[1];
  double s134 = M_T * G_T * TMath::Tan(y2) + pow(M_T,2.);
  //double range2 = 2000;
  //double s134 = pow(M_T,2.)-1000 + range2*Xarg[1];

  //cout << "y2=" << y2 << ", m134=" << TMath::Sqrt(s134) << endl;

  //double range3 = TMath::Pi();
  //double y3 = - TMath::Pi()/2. + range3 * Xarg[2];
  double range3 = TMath::Pi()/2. + TMath::ATan(M_W/G_W);
  double y3 = - TMath::ATan(M_W/G_W) + range3 * Xarg[2];
  double s25 = M_W * G_W * TMath::Tan(y3) + pow(M_W,2.);
  //double range3 = 2000;
  //double s25 = pow(M_W,2.)-1000 + range3*Xarg[2];

  //cout << "y3=" << y3 << ", m25=" << TMath::Sqrt(s25) << endl;

  //double range4 = TMath::Pi();
  //double y4 = - TMath::Pi()/2. + range4 * Xarg[3];
  double range4 = TMath::Pi()/2. + TMath::ATan(M_T/G_T);
  double y4 = - TMath::ATan(M_T/G_T) + range4 * Xarg[3];
  double s256 = M_T * G_T * TMath::Tan(y4) + pow(M_T,2.);
  //double range4 = 2000;
  //double s256 = pow(M_T,2.) -1000 + range4*Xarg[3];

  //cout << "y4=" << y4 << ", m256=" << TMath::Sqrt(s256) << endl;

  // #### FOR NWA
  /*s13 = pow(M_W,2.);
  s134 = pow(M_T,2.);
  s25 = pow(M_W,2.);
  s256 = pow(M_T,2.);*/

  // pT = transverse total momentum of the visible particles
  //TLorentzVector pT = p3 + p4 + p5 + p6;
  TLorentzVector pT = -Met;

  double p34 = p3*p4;
  double p56 = p5*p6;
  double p33 = p3.M2();
  double p44 = p4.M2();
  double p55 = p5.M2();
  double p66 = p6.M2();

  //cout << "p34=" << p34 << ",p56=" << p56 << ",p33=" << p33 << ",p44=" << p44 << ",p55" << p55 << ",p66" << p66;
  //cout << "Met=" << Met.Pt() << endl;

  // A1 p1x + B1 p1y + C1 = 0, with C1(E1,E2)
  // A2 p1y + B2 p2y + C2 = 0, with C2(E1,E2)
  // ==> express p1x and p1y as functions of E1, E2
  
  double A1 = 2.*( -p3.Px() + p3.Pz()*p4.Px()/p4.Pz() );
  double A2 = 2.*( p5.Px() - p5.Pz()*p6.Px()/p6.Pz() );
  
  double B1 = 2.*( -p3.Py() + p3.Pz()*p4.Py()/p4.Pz() );
  double B2 = 2.*( p5.Py() - p5.Pz()*p6.Py()/p6.Pz() );

  double Dx = B2*A1 - B1*A2;
  double Dy = A2*B1 - A1*B2;

  double X = 2*( pT.Px()*p5.Px() + pT.Py()*p5.Py() - p5.Pz()/p6.Pz()*( 0.5*(s25 - s256 + p66) + p56 + pT.Px()*p6.Px() + pT.Py()*p6.Py() ) ) + p55 - s25;
  double Y = p3.Pz()/p4.Pz()*( s13 - s134 + 2*p34 + p44 ) - p33 + s13;

  // p1x = alpha1 E1 + beta1 E2 + gamma1
  // p1y = ...(2)
  // p1z = ...(3)
  // p2z = ...(4)
  // p2x = ...(5)
  // p2y = ...(6)
  
  double alpha1 = -2*B2*(p3.E() - p4.E()*p3.Pz()/p4.Pz())/Dx;
  double beta1 = 2*B1*(p5.E() - p6.E()*p5.Pz()/p6.Pz())/Dx;
  double gamma1 = B1*X/Dx + B2*Y/Dx;

  double alpha2 = -2*A2*(p3.E() - p4.E()*p3.Pz()/p4.Pz())/Dy;
  double beta2 = 2*A1*(p5.E() - p6.E()*p5.Pz()/p6.Pz())/Dy;
  double gamma2 = A1*X/Dy + A2*Y/Dy;

  double alpha3 = (p4.E() - alpha1*p4.Px() - alpha2*p4.Py())/p4.Pz();
  double beta3 = -(beta1*p4.Px() + beta2*p4.Py())/p4.Pz();
  double gamma3 = ( 0.5*(s13 - s134 + p44) + p34 - gamma1*p4.Px() - gamma2*p4.Py() )/p4.Pz();

  double alpha4 = (alpha1*p6.Px() + alpha2*p6.Py())/p6.Pz();
  double beta4 = (p6.E() + beta1*p6.Px() + beta2*p6.Py())/p6.Pz();
  double gamma4 = ( 0.5*(s25 - s256 + p66) + p56 + (gamma1 + pT.Px())*p6.Px() + (gamma2 + pT.Py())*p6.Py() )/p6.Pz();

  double alpha5 = -alpha1;
  double beta5 = -beta1;
  double gamma5 = -pT.Px() - gamma1;

  double alpha6 = -alpha2;
  double beta6 = -beta2;
  double gamma6 = -pT.Py() - gamma2;

  // a11 E1^2 + a22 E2^2 + a12 E1E2 + a10 E1 + a01 E2 + a00 = 0
  // id. with bij

  double a11 = -1 + ( pow(alpha1,2.) + pow(alpha2,2.) + pow(alpha3,2.) );
  double a22 = pow(beta1,2.) + pow(beta2,2.) + pow(beta3,2.);
  double a12 = 2.*( alpha1*beta1 + alpha2*beta2 + alpha3*beta3 );
  double a10 = 2.*( alpha1*gamma1 + alpha2*gamma2 + alpha3*gamma3 );
  double a01 = 2.*( beta1*gamma1 + beta2*gamma2 + beta3*gamma3 );
  double a00 = pow(gamma1,2.) + pow(gamma2,2.) + pow(gamma3,2.);

  double b11 = pow(alpha5,2.) + pow(alpha6,2.) + pow(alpha4,2.);
  double b22 = -1 + ( pow(beta5,2.) + pow(beta6,2.) + pow(beta4,2.) );
  double b12 = 2.*( alpha5*beta5 + alpha6*beta6 + alpha4*beta4 );
  double b10 = 2.*( alpha5*gamma5 + alpha6*gamma6 + alpha4*gamma4 );
  double b01 = 2.*( beta5*gamma5 + beta6*gamma6 + beta4*gamma4 );
  double b00 = pow(gamma5,2.) + pow(gamma6,2.) + pow(gamma4,2.);

  // Find the intersection of the 2 ellipses (at most 4 real solutions for (E1,E2))
  vector<double> E1, E2;
  //cout << "coefs=" << a11 << "," << a22 << "," << a12 << "," << a10 << "," << a01 << "," << a00 << endl;
  //cout << "coefs=" << b11 << "," << b22 << "," << b12 << "," << b10 << "," << b01 << "," << b00 << endl;
  solve2Quads(a11, a22, a12, a10, a01, a00, b11, b22, b12, b10, b01, b00, E1, E2);

  // For each solution (E1,E2), find the neutrino 4-momenta p1,p2, find the initial quark momenta,
  // evaluate the matrix element and the jacobian

  if(E1.size() == 0){
	//cout << "No solutions!" << endl;
	mycount++;
	return 0.;
  }
  double integrand = 0.;

  //integrand = BreitWigner(s13,M_W,G_W) * BreitWigner(s25,M_W,G_W) * BreitWigner(s134,M_T,G_T) * BreitWigner(s256,M_T,G_T);

  /*unsigned int i = 0;
  if(E1.at(0) > E1.at(1))
	i = 0;
  else
    i = 1;*/

  for(unsigned int i=0; i<E1.size(); i++){
	//break;

	double e1 = E1.at(i);
	double e2 = E2.at(i);

	//cout << endl << "## Evaluating Matrix Element based on solutions e1 = " << e1 << ", e2 = " << e2 << endl << endl;

	if(e1 < 0. || e2 < 0.){
		mycount++;
		//cout << "Neg energies." << endl;
		//continue;
		return 0.;
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

	TLorentzVector p13 = p1 + p3;
	TLorentzVector p134 = p1 + p3 + p4;
	TLorentzVector p25 = p2 + p5;
	TLorentzVector p256 = p2 + p5 + p6;

	/*cout << "Input: W+ mass=" << TMath::Sqrt(s13) << ", Top mass=" << TMath::Sqrt(s134) << ", W- mass=" << TMath::Sqrt(s25) << ", Anti-top mass=" << TMath::Sqrt(s256) << endl;
	cout << "Output: W+ mass=" << p13.M() << ", Top mass=" << p134.M() << ", W- mass=" << p25.M() << ", Anti-top mass=" << p256.M() << endl << endl;
	
	cout << "Electron (E,Px,Py,Pz) = ";
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
  
	TLorentzVector tot = p1 + p2 + p3 + p4 + p5 + p6;

	double ETot = tot.E();
	double PzTot = tot.Pz();

	double q1Pz = (PzTot + ETot)/2.;
	double q2Pz = (PzTot - ETot)/2.;

	//cout << "===> Eext=" << ETot << ", Pzext=" << PzTot << ", q1Pz=" << q1Pz << ", q2Pz=" << q2Pz << endl << endl;
  
    if(q1Pz > SQRT_S/2. || q2Pz < -SQRT_S/2. || q1Pz < 0. || q2Pz > 0.){
      //cout << "Fail!" << endl;
      mycount++;
      //continue;
	  return 0.;
    }
	
	// momentum vector definition
  	vector<double*> p(1, new double[4]);
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

  	// Compute the Pdfs
  	//double pdf1_1 = ComputePdf(21,TMath::Abs(q1Pz/(SQRT_S/2.)), pow(tot.M(),2)) / TMath::Abs(q1Pz/(SQRT_S/2.));
  	//double pdf1_2 = ComputePdf(21,TMath::Abs(q2Pz/(SQRT_S/2.)), pow(tot.M(),2)) / TMath::Abs(q2Pz/(SQRT_S/2.));
  	double pdf1_1 = ComputePdf(21,TMath::Abs(q1Pz/(SQRT_S/2.)), pow(M_T,2)) / TMath::Abs(q1Pz/(SQRT_S/2.));
  	double pdf1_2 = ComputePdf(21,TMath::Abs(q2Pz/(SQRT_S/2.)), pow(M_T,2)) / TMath::Abs(q2Pz/(SQRT_S/2.));
  
	// Compute flux factor 1/(2*x1*x2*s)
	double PhaseSpaceIn = 1.0 / ( 2. * TMath::Abs(q1Pz/(SQRT_S/2.)) * TMath::Abs(q2Pz/(SQRT_S/2.0)) * pow(SQRT_S,2)); 

	// Compute finale Phase Space for observed particles (not concerned by the change of variable)
	// dPhi = |P|^2 sin(theta)/(2*E*(2pi)^3)
  	double dPhip3 = pow(p3.P(),2.)*TMath::Sin(p3.Theta())/(2.0*p3.E()*pow(2.*TMath::Pi(),3));
  	double dPhip4 = pow(p4.P(),2.)*TMath::Sin(p4.Theta())/(2.0*p4.E()*pow(2.*TMath::Pi(),3));
  	double dPhip5 = pow(p5.P(),2.)*TMath::Sin(p5.Theta())/(2.0*p5.E()*pow(2.*TMath::Pi(),3));
  	double dPhip6 = pow(p6.P(),2.)*TMath::Sin(p6.Theta())/(2.0*p6.E()*pow(2.*TMath::Pi(),3));
	double PhaseSpaceOut = dPhip5 * dPhip6 * dPhip3 * dPhip4;

	// Set momenta for this event
  	process.setMomenta(p);

	// Compute jacobian from change of variable:
	vector<TLorentzVector> momenta;
	momenta.push_back(p1);
	momenta.push_back(p2);
	momenta.push_back(p3);
	momenta.push_back(p4);
	momenta.push_back(p5);
	momenta.push_back(p6);
	double jac = computeJacobian(momenta);
	if(jac <= 0.){
		mycount++;
		continue;
	}
	jac /= 8.*16.*pow(TMath::Pi()*SQRT_S,2.);

  	// Evaluate matrix element
  	process.sigmaKin();
  	const double* matrix_elements1 = process.getMatrixElements();

	//cout << "Found PDF1 = " << pdf1_1 << ", PDF2 = " << pdf1_2 << ", PS in = " << PhaseSpaceIn << ", PS out = " << PhaseSpaceOut << ", jac = " << jac << endl;
	//cout << "===> Matrix element = " << matrix_elements1[0] << ", prod = " << PhaseSpaceIn * matrix_elements1[0] * pdf1_1 * pdf1_2 * PhaseSpaceOut * jac << endl << endl ;	
  
	integrand += PhaseSpaceIn * matrix_elements1[0] * pdf1_1 * pdf1_2 * PhaseSpaceOut * jac;

	// free up memory
	for(unsigned int i = 0; i < p.size(); ++i){
		delete p.at(i); p.at(i) = 0;
	}
  }

  if(integrand == 0.){
	mycount++;
    //cout << "Zero!" << endl;
	return integrand;
  }

  double flatterJac = range1 * range2 * range3 * range4;
  flatterJac *= M_W*G_W * M_T*G_T * M_W*G_W * M_T*G_T;
  flatterJac /= pow(TMath::Cos(y1) * TMath::Cos(y2) * TMath::Cos(y3) * TMath::Cos(y4), 2.);

  // ### FOR NWA
  //double flatterJac = pow(TMath::Pi(),4.) * (M_W*G_W * M_T*G_T * M_W*G_W * M_T*G_T);

  //cout << "## Phase Space point done. Integrand = " << integrand << ", flatterjac = " << flatterJac << ", prod = " << integrand*flatterJac <<  endl;

  integrand *= flatterJac;

  return integrand;
}


double TFDISTR::ComputePdf(int pid, double x, double q2){
    // return the xf(pid,x,q2), be careful: this value must be divided by x to obtain f
    const double xf = pdf->xfxQ2(pid, x, q2);
    return xf;
}

double ME(double *error, TLorentzVector ep, TLorentzVector mum, TLorentzVector b, TLorentzVector bbar, TLorentzVector Met, int nCells = 2000, int nSampl = 100, int nPoints = 50000){
  
  TH1D *hst_Wm = new TH1D("test_mu", "test_1D", 150,0,150);
  TH1D *hst_We = new TH1D("test_ep", "test_1D", 150,0,150);
  TH1D *hst_t = new TH1D("test_t", "test_1D", 100,150,250);
  TH1D *hst_tbar = new TH1D("test_tbar", "test_1D", 100,150,250);

  TStopwatch chrono;
  std::cout << "initialising : " ;
  // TFoam Implementation
  double *MCvect = new double[4];
  TRandom3  *PseRan   = new TRandom3();
  //PseRan->SetSeed(514638543);  
  TFoam   *FoamX    = new TFoam("FoamX");
  FoamX->SetkDim(4);
  FoamX->SetnCells(nCells);      // No. of cells, can be omitted, default=2000
  FoamX->SetnSampl(nSampl);

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

  for(Long_t loop=0; loop<nPoints; loop++){
    int count_old = mycount;
    FoamX->MakeEvent();          // generate MC event
    FoamX->GetMCvect( MCvect );   // get generated vector (x,y)
	
	double range1 = TMath::Pi()/2. + TMath::ATan(M_W/G_W);
  	double y1 = - TMath::ATan(M_W/G_W) + range1 * MCvect[0];
  	double s13 = M_W * G_W * TMath::Tan(y1) + pow(M_W,2.);
	//double range1 = 2000;
	//double s13 = pow(M_W,2.) -1000+ range1*MCvect[0];

	double range2 = TMath::Pi()/2. + TMath::ATan(M_T/G_T);
  	double y2 = - TMath::ATan(M_T/G_T) + range2 * MCvect[1];
  	double s134 = M_T * G_T * TMath::Tan(y2) + pow(M_T,2.);
	//double range2 = 2000;
	//double s134 = pow(M_T,2.) -1000+ range2*MCvect[1];

	double range3 = TMath::Pi()/2. + TMath::ATan(M_W/G_W);
  	double y3 = - TMath::ATan(M_W/G_W) + range3 * MCvect[2];
  	double s25 = M_W * G_W * TMath::Tan(y3) + pow(M_W,2.);
	//double range3 = 2000;
	//double s25 = pow(M_W,2.)-1000 + range3*MCvect[2];

	double range4 = TMath::Pi()/2. + TMath::ATan(M_T/G_T);
  	double y4 = - TMath::ATan(M_T/G_T) + range4 * MCvect[3];
  	double s256 = M_T * G_T * TMath::Tan(y4) + pow(M_T,2.);
	//double range4 = 2000;
	//double s256 = pow(M_T,2.)-1000 + range4*MCvect[3];

	if(count_old == mycount){
		hst_We->Fill(TMath::Sqrt(s13));
		hst_Wm->Fill(TMath::Sqrt(s25));
		hst_t->Fill(TMath::Sqrt(s134));
		hst_tbar->Fill(TMath::Sqrt(s256));
	}
  }

  double mcResult, mcError;
  FoamX->GetIntegMC( mcResult, mcError);  // get MC integral, should be one
  
  std::cout << "CPU time : " << chrono.CpuTime() << "  Real-time : " << chrono.RealTime() << std::endl;
  cout << "nr. fails: " << mycount  << endl;
  mycount = 0;

  cout << " mcResult= " << mcResult << " +- " << mcError <<endl;

  
  TCanvas *c = new TCanvas("c","Canvas for plotting",600,600);
  c->cd();
  hst_We->Draw();
  //c->Print(TString("plots/")+SSTR(count_wgt)+"_"+SSTR(count_perm)+"_Enu.png");
  delete hst_We; hst_We = 0;
  delete c; c = 0;
  
  c = new TCanvas("c","Canvas for plotting",600,600);
  c->cd();
  hst_Wm->Draw();
  //c->Print(TString("plots/")+SSTR(count_wgt)+"_"+SSTR(count_perm)+"_Munu.png");
  delete hst_Wm; hst_Wm = 0;
  delete c; c = 0;


  c = new TCanvas("c","Canvas for plotting",600,600);
  c->cd();
  hst_t->Draw();
  //c->Print(TString("plots/")+SSTR(count_wgt)+"_"+SSTR(count_perm)+"_t.png");
  delete hst_t; hst_t = 0;
  delete c; c = 0;
  
  c = new TCanvas("c","Canvas for plotting",600,600);
  c->cd();
  hst_tbar->Draw();
  //c->Print(TString("plots/")+SSTR(count_wgt)+"_"+SSTR(count_perm)+"_tbar.png");
  delete hst_tbar; hst_tbar = 0;
  delete c; c = 0;

  delete FoamX; FoamX = 0;
  delete PseRan; PseRan = 0;
  delete MCvect; MCvect = 0;

  *error = mcError;
  if(std::isnan(*error))
	*error = 0.;
  if(std::isnan(mcResult))
	mcResult = 0.;
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

  // Get pointers to branches used in this analysis
  TClonesArray *branchGen = treeReader->UseBranch("Particle");

  cout << "Entries:" << chain.GetEntries() << endl;

  ofstream fout(outputFile);

  // Full weights
  //double madweight1[10] = {1.58495292058e-21, 2.09681384879e-21, 4.34399623629e-22, 1.68163897955e-22, 3.20350498956e-22, 5.22232034307e-22, 6.04738375743e-21, 9.55643564854e-22, 8.12425265344e-22, 5.81210532053e-23};
  //double madweight2[10] = {1.02514966131e-21, 1.45375719248e-21, 1.65080839221e-22, 1.55653414654e-24, 5.60531044001e-25, 1., 9.70526105314e-24, 3.89103636371e-22, 6.38206925825e-23, 9.37189585544e-26};
  double madweight1[10] = {1.48990458909e-21,2.00433822978e-21,4.08998078881e-22,1.56339237714e-22,2.98606743727e-22,4.79498317117e-22,5.63645701583e-21,8.99177777775e-22,7.68316733666e-22,5.42606461617e-23};
  double madweight1Err[10] = {8.63813589113e-24,1.08426062115e-23,2.5750146827e-24,7.0506407196e-25,1.10554655068e-24,2.31140842678e-24,2.71677566322e-23,4.8290429288e-24,1.69718762833e-24,2.66346844676e-25};
  double madweight2[10] = {9.62646303545e-22,1.38143123163e-21,1.54526017444e-22,1.45628835295e-24,6.80263123625e-25,1.,1.07797730384e-23,3.61278172744e-22,6.19087950579e-23,7.20276231557e-26};
  double madweight2Err[10] = {2.96180414077e-24,4.8856162625e-24,1.0218999515e-24,1.29754825587e-25,2.72733072519e-25,1.,4.03010515215e-24,4.29592188061e-24,1.67765665953e-24,8.06569780018e-27};

  // NWA weights
  //double madweight1[10] = {1.26069381322e-21, 2.85437676736e-21, 4.8136572269e-22, 1., 3.99894656854e-22, 5.7603822256e-22, 6.99323258475e-21, 1.0892124248e-21, 8.28291668972e-22, 1.};
  //double madweight2[10] = {1.46272073513e-21, 1.51733772927e-21, 1.61193875253e-22, 1., 1., 1., 1., 1., 1., 1.};
  // NWA weights, first sol
  //double madweight1[3] = {8.93501179418e-22, 7.42359826601e-22, 1.49577468649e-22};
  //double madweight2[3] = {1.04113131882e-23, 7.04643552065e-22, 4.3214935529e-23};
  // NWA weights, first sol, ME=1
  //double madweight1[3] = {1.9968889994e-17, 1.10734832869e-17, 2.17966664756e-18};
  //double madweight2[3] = {1.2718723458e-19, 2.38734853175e-17, 6.27800021816e-19};

  for(int entry = 0; entry < 10 ; ++entry){//numberOfEntries; ++entry)
  for(int permutation = 1; permutation <= 2; permutation ++){

	//entry = 3; permutation = 2;

  	count_perm = permutation;
  
	double weight1 = 0/*, weight2 = 0*/;

  	vector<double> nbr_points;
  	vector<double> wgts;
  	vector<double> wgtsErr;
  	vector<double> wgtsErrX;
  	double max_wgt = 0.;
  	double min_wgt = 1.;

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    TLorentzVector gen_ep, gen_mum, gen_b, gen_bbar, gen_Met;

    GenParticle *gen;

	//int count_ep=0, count_mum=0;

    for (Int_t i = 0; i < branchGen->GetEntries(); i++){
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

  	/*weight1 = ME(gen_ep, gen_mum, gen_b, gen_bbar, gen_Met, 20000, 10, 50000)/2.;
  	weight2 = ME(gen_ep, gen_mum, gen_bbar, gen_b, gen_Met, 20000, 10, 50000)/2.;
  	
  	double weight3;
  	if(weight1 < weight2){
  	  weight3 = weight2;
  	  weight2 = weight1;
  	  weight1 = weight3;
  	}
  	if(madweight1[entry] < madweight2[entry]){
  	  weight3 = madweight2[entry];
  	  madweight2[entry] = madweight1[entry];
  	  madweight1[entry] = weight3;
  	  cout << "entry " << madweight1[entry] << "," << madweight2[entry] << endl;
  	}

  	fout << entry << " weight1 = " << weight1 << " (" << weight1 / (double) madweight1[entry] << "), weight2 = " << weight2 << " (" << weight2/madweight2[entry] << ")" << endl;
  	
  	count_wgt++;*/
  	
  	for(int k = 1; k <= 51; k+=5){
		/*int nCells = k*100;
		int nSampl = 50;
		int nPoints = k*2000;*/
		int nCells = k*40;
		int nSampl = 200;
		int nPoints = 50000;
		/*int nCells = 1000;
		int nSampl = 20;
		int nPoints = 40000;*/

		double error = 0;

  		if(permutation == 1)
  		  weight1 = ME(&error, gen_ep, gen_mum, gen_b, gen_bbar, gen_Met, nCells, nSampl, nPoints)/2;
  		if(permutation == 2)
  		  weight1 = ME(&error, gen_ep, gen_mum, gen_bbar, gen_b, gen_Met, nCells, nSampl, nPoints)/2;
  		
  		nbr_points.push_back(nCells);
  		//nbr_points.push_back(nPoints);
  		wgts.push_back(weight1);
  		wgtsErr.push_back(error);
  		wgtsErrX.push_back(0.);
  		if(weight1 > max_wgt) max_wgt = weight1;
  		if(weight1 < min_wgt) min_wgt = weight1;
  		
  		if(abs(weight1-madweight1[entry]) > abs(madweight2[entry]-weight1)){
  		  double weight3 = madweight2[entry];
  		  madweight2[entry] = madweight1[entry];
  		  madweight1[entry] = weight3;
  		  
  		  weight3 = madweight2Err[entry];
  		  madweight2Err[entry] = madweight1Err[entry];
  		  madweight1Err[entry] = weight3;
  		}

  		fout << entry << "/" << permutation << ", points=" << nPoints << ", cells=" << nCells << ", evts/cell=" << nSampl << ": weight = " << weight1 << " +- " << error << " (" << weight1 / (double) madweight1[entry] << ")" <<  endl;

		//break;

  	}

  	double mwX[2] = {nbr_points.at(0), nbr_points.at(nbr_points.size()-1)};
  	double mwErrX[2] = {0.,0.};
  	double correct = madweight1[entry];
  	double correctErr = madweight1Err[entry];
  	if(abs(correct-wgts.at(wgts.size()-1)) > abs(madweight2[entry]-wgts.at(wgts.size()-1))){
  	  correct = madweight2[entry];
  	  correctErr = madweight2Err[entry];
  	}
  	double madwgts[2] = {correct, correct};
  	double madwgtsErr[2] = {correctErr, correctErr};
  	if(correct < min_wgt)
  	  min_wgt = correct;
  	if(correct > max_wgt)
  	  max_wgt = correct;
  	
  	TGraphErrors* wgt_vs_points = new TGraphErrors(wgts.size(), &nbr_points[0], &wgts[0], &wgtsErrX[0], &wgtsErr[0]);
  	TGraphErrors* madwgt_vs_points = new TGraphErrors(2, mwX, madwgts, mwErrX, madwgtsErr);
  	TCanvas* c = new TCanvas("c","Canvas for plotting",600,600);
  	c->cd();
  	wgt_vs_points->Draw("AC*");
  	wgt_vs_points->GetHistogram()->SetMaximum(max_wgt+wgtsErr[0]);
  	wgt_vs_points->GetHistogram()->SetMinimum(min_wgt-wgtsErr[0]);
  	wgt_vs_points->SetMarkerColor(kRed);
  	wgt_vs_points->SetLineColor(kRed);
  	wgt_vs_points->SetMarkerStyle(21);
  	madwgt_vs_points->Draw("CP3");
  	madwgt_vs_points->SetMarkerColor(kBlue);
  	madwgt_vs_points->SetLineColor(kBlue);
  	madwgt_vs_points->SetFillColor(kBlue);
  	madwgt_vs_points->SetFillStyle(3005);
  	c->Print(TString("plots/wgt_vs_cells_50000p_2000c_200s_")+SSTR(entry)+"_"+SSTR(permutation)+".png");
  	//c->Print(TString("plots/wgt_vs_points_102000p_5100c_50s_")+SSTR(entry)+"_"+SSTR(permutation)+".png");
  	delete wgt_vs_points; wgt_vs_points = 0;
  	delete madwgt_vs_points; madwgt_vs_points = 0;
  	delete c;
  
  //break;
  }
  count_wgt++;
  //break;
  }

  delete treeReader; 
}

