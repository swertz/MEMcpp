#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdlib.h>

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

#define SQRT_S 8000 

#define SSTR( x ) dynamic_cast< std::ostringstream & > \
				( std::ostringstream() << std::dec << x ).str()

using namespace LHAPDF;
using namespace std;

int mycount = 0, count_wgt = 0, count_perm=1;

unsigned int setFlags(char verbosity = 0, bool subregion = false, bool retainStateFile = false, unsigned int level = 0, bool smoothing = false, bool takeOnlyGridFromFile = true){
	/* Set option flags for CUBA integrator
	 * 
	 * smoothing only used by Suave
	 * smoothing and takeOnlyGridFromFile only used by Vegas
	 *
	 * verbosity = 0-3
	 * subregion true = only last set of samples is used for final evaluation of integral
	 * smoothing true = smoothe importance function
	 * retainStateFile true => retain state file when integration ends
	 * takeOnlyGridFromFile false => full state taken from file (if present), true => only grid is taken (e.g. to use it for another integrand)
	 * level determines random-number generator:
	 *	seed = 0 => Sobol quasi-random
	 *	seed != 0 => level is used:
	 *		level = 0 => Mersenne Twister pseudo-random
	 *		level != 0 => Ranlux pseudo-random, level = determines period of generator
	*/

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

class MEWeight{
	public:

	MEWeight(const string paramCardPath, const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met);
	inline double ComputePdf(const int pid, const double x, const double q2);

	inline TLorentzVector GetP3(void) const { return p3; }
	inline TLorentzVector GetP4(void) const { return p4; }
	inline TLorentzVector GetP5(void) const { return p5; }
	inline TLorentzVector GetP6(void) const { return p6; }
	inline TLorentzVector GetMet(void) const { return Met; }

	inline void setProcessMomenta(vector<double*> &p){ process.setMomenta(p); }
	inline void computeMatrixElements(){ process.sigmaKin(); }
	inline const double* const getMatrixElements() const { return process.getMatrixElements(); }

	private:

	TLorentzVector p3, p4, p5, p6, Met;

	CPPProcess process;
	PDF *pdf;
};

MEWeight::MEWeight(const string paramCardPath, const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met){
	p3 = ep;
	p5 = mum;
	p4 = b;
	p6 = bbar;
	Met = met;

	process.initProc(paramCardPath);
	pdf = mkPDF("cteq6l1", 0);
}

double MEWeight::ComputePdf(const int pid, const double x, const double q2){
	// return f(pid,x,q2)
	if(x <= 0 || x >= 1)
		return 0.;
	else
		return pdf->xfxQ2(pid, x, q2)/x;
}

int BWTest(const int *nDim, const double* Xarg, const int *nComp, double *Value, void *inputs){
	*Value = 0.;

	for(int i=0; i<*nDim; ++i){
		if(Xarg[i] == 1. || Xarg[i] == 0.){
			mycount++;
			return 0;
		}
	}

	double range1 = TMath::Pi();
	double y1 = - TMath::Pi()/2. + range1 * Xarg[0];
	const double s13 = M_W * G_W * TMath::Tan(y1) + pow(M_W,2.);

	double range2 = TMath::Pi();
	double y2 = - TMath::Pi()/2. + range2 * Xarg[1];
	const double s134 = M_T * G_T * TMath::Tan(y2) + pow(M_T,2.);

	double range3 = TMath::Pi();
	double y3 = - TMath::Pi()/2. + range3 * Xarg[2];
	const double s25 = M_W * G_W * TMath::Tan(y3) + pow(M_W,2.);

	double range4 = TMath::Pi();
	double y4 = - TMath::Pi()/2. + range4 * Xarg[3];
	const double s256 = M_T * G_T * TMath::Tan(y4) + pow(M_T,2.);
	
	*Value = BreitWigner(s13,M_W,G_W) * BreitWigner(s25,M_W,G_W) * BreitWigner(s134,M_T,G_T) * BreitWigner(s256,M_T,G_T);
	
	double flatterJac = range1 * range2 * range3 * range4;
	flatterJac *= M_W*G_W * M_T*G_T * M_W*G_W * M_T*G_T;
	flatterJac /= pow(TMath::Cos(y1) * TMath::Cos(y2) * TMath::Cos(y3) * TMath::Cos(y4), 2.);
	
	flatterJac /= pow(TMath::Pi(),4);

	*Value *= flatterJac;

	return 0;
}

int MEFunct(const int *nDim, const double* Xarg, const int *nComp, double *Value, void *inputs){

	//cout << endl << endl << endl << "########## Starting phase-space point ############" << endl << endl;
	//cout << "Inputs = [" << Xarg[0] << "," << Xarg[1] << "," << Xarg[2] << "," << Xarg[3] << endl;

	MEWeight* myWeight = (MEWeight*) inputs;

	TLorentzVector p3 = myWeight->GetP3();
	TLorentzVector p4 = myWeight->GetP5();
	TLorentzVector Met = myWeight->GetMet();

	*Value = 0.;

	for(int i=0; i<*nDim; ++i){
		if(Xarg[i] == 1. || Xarg[i] == 0.){
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

	const double range2 = TMath::Pi()/2. + TMath::ATan(M_W/G_W);
	const double y2 = - TMath::ATan(M_W/G_W) + range2 * Xarg[1];
	const double s24 = M_W * G_W * TMath::Tan(y2) + pow(M_W,2.);


	//const double q1 = Xarg[2]/50.0;
	//const double q2 = -Xarg[3]/50.0;
	//
	const double q1 = 0.0615737;
	const double q2 = 0.0082855;

	//cout << "y2=" << y2 << ", M24=" << TMath::Sqrt(s24) << endl;

	//if(s13 > s134 || s25 > s256 || s13 < p3.M() || s25 < p5.M()){
		//cout << "Masses too small!" << endl;
	//	mycount++;
	//	return 0;
	//}


	// pb = transverse total momentum of the visible particles
	const TLorentzVector pb = p3+p4;

	// P2x = a1 E2 + a2 P2y + a3
	// P2z = b1 E2 + b2 P2y + b3

	const double Qm = SQRT_S*(q1-q2)/2;
	const double Qp = SQRT_S*(q1+q2)/2;
	
        const double ka = -4*p4.Px()*p3.Pz()+4*p3.Px()+p4.Pz();
	const double kb = 2*(p3.Pz()*p4.Px()-p3.Px()*p4.Pz());

	const double b1 = -(1/kb) * (2*p4.E()*p3.Px()-2*p3.E()*p4.Px());
	const double b2 = -(1/kb) * (2*p3.Py()*p4.Px()-2*p3.Px()*p4.Py());
	const double b3 = -(1/kb) * (pow(p4.M(),2)*p3.Px()-2*p3.E()*pb.E()*p4.Px()+pow(p3.M(),2)*p4.Px()+2*p3.Px()*p4.Px()*pb.Px()+2*p3.Py()*p4.Px()*pb.Py()+2*p3.Pz()*p4.Px()*pb.Pz()-2*p3.Pz()*p4.Px()*Qm+2*p3.E()*p4.Px()*Qp-p4.Px()*s13-p3.Px()*s24);

	const double a1 = (1/ka)*(-2*p4.Pz()*(-2*p3.E())-2*p3.Pz()*2*p4.E()) ;
	const double a2 = (1/ka)*(-4*p3.Py()*p4.Pz()-2*p3.Pz()*(-2*p4.Py()));
	const double a3 = -pb.Px()+(1/ka)*(-4*p4.Px()*pb.Px()*p3.Pz()-4*p3.Py()*pb.Py()*p4.Pz()-2*p4.Pz()*(-2*p3.E()*(pb.E()-Qp)+pow(p3.M(),2)+2*p3.Pz()*(pb.Pz()-Qm)-s13)-2*p3.Pz()*(pow(p4.M(),2)-s24)) ;


	cout << "p2x : -51.7897 , p2y : 23.0622 , p2z : 112.36 , E2 : 125.852" << endl;
	cout << "p2x : 50.8451 , p2y : -21.9069 , p2z : 46.6542 , E2 : 72.3999" << endl;
	
	const double test = a1*125.852+a2*(23.0622)+a3;
	const double test2 = b1*125.852+b2*(23.0622)+b3;
        const double test3 = a1*72.3999+a2*(-21.9069)+a3;
        const double test4 = b1*72.3999+b2*(-21.9069)+b3;

	cout << " test p2x : " << test << endl;
        cout << " test p2z : " << test2 << endl;
	cout << " test p2x : " << test3 << endl;
	cout << " test p2z : " << test4 << endl;


	//cout << " ka : " << ka << endl;
	//cout << " kb : " << kb << endl; 

	//cout << " a : " << a1 << " " << a2 << " " << a3 << endl;
	//cout << " b : " << b1 << " " << b2 << " " << b3 << endl;

	// 0 = c1 E2 + c2 P2y + c3 
	// E2 = d1 P2y + d2
	
	const double c1 = 2*(Qp-pb.E())+2*pb.Px()*a1-2*(Qm-pb.Pz())*b1;
	const double c2 = 2*pb.Px()*a2+2*pb.Py()-2*(Qm-pb.Pz())*b2;
	const double c3 = -pow(Qp-pb.E(),2)+pow(pb.Px(),2)+2*pb.Px()*a3+pow(pb.Py(),2)+pow(Qm-pb.Pz(),2)-2*(Qm-pb.Pz())*b3;

	//cout << " c : " << c1 << " " << c2 << " " << c3 << endl;
	
	const double d1 = -c2/c1;
	const double d2 = -c3/c1;

	//cout << " d : " << d1 << " " << d2 << endl; 

	// alpha*P2y^2 + beta*P2y + gamma = 0
	
	const double alpha = pow(a1,2)*pow(d1,2)+pow(a2,2)+2*a1*a2*d1+1+pow(b1,2)*pow(d1,2)+pow(b2,2)+2*b1*b2*d1-pow(d1,2);
	const double beta  = 2*pow(a1,2)*d1*d2+2*a1*a3*d1+2*a1*a2*d2+2*a2*a3+2*pow(b1,2)*d1*d2+2*b1*b3*d1+2*b1*b2*d2+2*b2*b3-2*d1*d2;
	const double gamma = pow(a1,2)*pow(d2,2)+pow(a3,2)+2*a1*a3*d2+pow(b1,2)*pow(d2,2)+pow(b3,2)+2*b1*b3*d2-pow(d2,2);
	
	// Find P2Y
	vector<double> P2Y;

	solveQuadratic(alpha, beta, gamma, P2Y, false);

	//cout << "coefs=" << a11 << "," << a22 << "," << a12 << "," << a10 << "," << a01 << "," << a00 << endl;
	//cout << "coefs=" << b11 << "," << b22 << "," << b12 << "," << b10 << "," << b01 << "," << b00 << endl;

	// For each solution (E1,E2), find the neutrino 4-momenta p1,p2, find the initial quark momenta,
	// evaluate the matrix element and the jacobian
	
	if(P2Y.size() == 0){
		//cout << "No solutions!" << endl;
		mycount++;
		return 0;
	}

	cout << "Checking W Mass... " << P2Y.size() << " entries" << endl;
	for(unsigned int i=0; i<P2Y.size(); i++){
		const double P2X = a1*d1*P2Y.at(i)+a1*d2+a2*P2Y.at(i)+a3;
		const double P2Z = b1*d1*P2Y.at(i)+b1*d2+b2*P2Y.at(i)+b3;

		TLorentzVector P2, W;
		const double P2E = sqrt(pow(P2X,2)+pow(P2Y.at(i),2)+pow(P2Z,2));
		P2.SetPxPyPzE(P2X,P2Y.at(i),P2Z,P2E);
		W = P2+p4;

		cout << "  W24 mass : " << sqrt(s24) << " " <<  W.M();
		//cout << "W13 mass : " << s13 << " " <<  W.M() << endl;


	}
	cout << endl;


	/*
	for(unsigned int i=0; i<E1.size(); i++){
		if(E1.size() >= 0){
			if(E1.at(0) > E1.at(1))
				i = 0;
			else
				i = 1;
		}

		const double e1 = E1.at(i);
		const double e2 = E2.at(i);

		//cout << endl << "## Evaluating Matrix Element based on solutions e1 = " << e1 << ", e2 = " << e2 << endl << endl;

		if(e1 < 0. || e2 < 0.){
			//mycount++;
			//cout << "Neg energies." << endl;
			//continue;
			break;
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

	
		const TLorentzVector tot = p1 + p2 + p3 + p4 + p5 + p6;

		const double ETot = tot.E();
		const double PzTot = tot.Pz();

		const double q1Pz = (PzTot + ETot)/2.;
		const double q2Pz = (PzTot - ETot)/2.;

		//cout << "===> Eext=" << ETot << ", Pzext=" << PzTot << ", q1Pz=" << q1Pz << ", q2Pz=" << q2Pz << endl << endl;
	
		if(q1Pz > SQRT_S/2. || q2Pz < -SQRT_S/2. || q1Pz < 0. || q2Pz > 0.){
			//cout << "Fail!" << endl;
			mycount++;
			//continue;
			break;
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
		const double pdf1_1 = myWeight->ComputePdf(21,TMath::Abs(q1Pz/(SQRT_S/2.)), pow(M_T,2));
		const double pdf1_2 = myWeight->ComputePdf(21,TMath::Abs(q2Pz/(SQRT_S/2.)), pow(M_T,2));
	
		// Compute flux factor 1/(2*x1*x2*s)
		const double PhaseSpaceIn = 1.0 / ( 2. * TMath::Abs(q1Pz/(SQRT_S/2.)) * TMath::Abs(q2Pz/(SQRT_S/2.0)) * pow(SQRT_S,2)); 

		// Compute finale Phase Space for observed particles (not concerned by the change of variable)
		// dPhi = |P|^2 sin(theta)/(2*E*(2pi)^3)
		const double dPhip3 = pow(p3.P(),2.)*TMath::Sin(p3.Theta())/(2.0*p3.E()*pow(2.*TMath::Pi(),3));
		const double dPhip4 = pow(p4.P(),2.)*TMath::Sin(p4.Theta())/(2.0*p4.E()*pow(2.*TMath::Pi(),3));

		const double PhaseSpaceOut = dPhip3 * dPhip4;

		// Set momenta for this event
		myWeight->setProcessMomenta(p);

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
			//cout << "Jac infinite!" << endl;
			mycount++;
			//continue;
			break;
		}

		// Evaluate matrix element
		myWeight->computeMatrixElements();
		const double* const matrix_elements1 = myWeight->getMatrixElements();

		//cout << "Found PDF1 = " << pdf1_1 << ", PDF2 = " << pdf1_2 << ", PS in = " << PhaseSpaceIn << ", PS out = " << PhaseSpaceOut << ", jac = " << jac << endl;
		//cout << "===> Matrix element = " << matrix_elements1[0] << ", prod = " << PhaseSpaceIn * matrix_elements1[0] * pdf1_1 * pdf1_2 * PhaseSpaceOut * jac << endl << endl ;	
		
		*Value += PhaseSpaceIn * matrix_elements1[0] * pdf1_1 * pdf1_2 * PhaseSpaceOut * jac;

		// free up memory
		for(unsigned int i = 0; i < p.size(); ++i){
			delete p.at(i); p.at(i) = 0;
		}

		break;
	}

	if(*Value == 0.){
		mycount++;
		//cout << "Zero!" << endl;
		return 0;
	}

	double flatterJac = range1 * range2 * range3 * range4;
	flatterJac *= M_W*G_W * M_T*G_T * M_W*G_W * M_T*G_T;
	flatterJac /= pow(TMath::Cos(y1) * TMath::Cos(y2) * TMath::Cos(y3) * TMath::Cos(y4), 2.);

	// ### FOR NWA
	//double flatterJac = pow(TMath::Pi(),4.) * (M_W*G_W * M_T*G_T * M_W*G_W * M_T*G_T);

	//cout << "## Phase Space point done. Integrand = " << integrand << ", flatterjac = " << flatterJac << ", prod = " << integrand*flatterJac <<	endl;

	*Value *= flatterJac;
	*/
	return 0;
}

double ME(double *error, TLorentzVector ep, TLorentzVector mum, TLorentzVector Met, double *time){
	
	/*TH1D *hst_Wm = new TH1D("test_mu", "test_1D", 150,0,150);
	TH1D *hst_We = new TH1D("test_ep", "test_1D", 150,0,150);
	TH1D *hst_t = new TH1D("test_t", "test_1D", 100,150,250);
	TH1D *hst_tbar = new TH1D("test_tbar", "test_1D", 100,150,250);*/

	cout << "Initializing integration..." << endl;
	
	TStopwatch chrono;
	chrono.Start();

	TLorentzVector b, bbar;

	MEWeight myWeight("/home/fynu/swertz/scratch/Madgraph/madgraph5/cpp_ttbar_epmum/Cards/param_card.dat", ep, mum, b, bbar, Met);

	int neval, nfail;
	double mcResult=0, mcError=0, prob=0;
	
	char verbosity = 0; // 0-3
	bool subregion = false; // true = only last set of samples is used for final evaluation of integral
	bool smoothing = false;
	bool retainStateFile = false; // false => delete state file when integration ends
	bool takeOnlyGridFromFile = true; // false => full state taken from file (if present), true => only grid is taken (e.g. to use it for another integrand)
	unsigned int level = 0; 

	unsigned int flags = setFlags(verbosity, subregion, retainStateFile, level, smoothing, takeOnlyGridFromFile);

	cout << "Starting integration..." << endl;
	cubacores(0,0);	
	Vegas(
		4,						// (int) dimensions of the integrated volume
		1,						// (int) dimensions of the integrand
		(integrand_t) MEFunct,	// (integrand_t) integrand (cast to integrand_t)
		//(integrand_t) BWTest,	// (integrand_t) integrand (cast to integrand_t)
		(void*) &myWeight,		// (void*) pointer to additional arguments passed to integrand
		1,						// (int) maximum number of points given the integrand in each invocation (=> SIMD) ==> PS points = vector of sets of points (x[ndim][nvec]), integrand returns vector of vector values (f[ncomp][nvec])
		0.005,					// (double) requested relative accuracy |-> error < max(rel*value,abs)
		0.,						// (double) requested absolute accuracy |
		flags,					// (int) various control flags in binary format, see setFlags function
		8945,						// (int) seed (seed==0 && no control flag => SOBOL; seed!=0 && control flag==0 => Mersenne Twister)
		0,					// (int) minimum number of integrand evaluations
		500,					// (int) maximum number of integrand evaluations (approx.!)
		1000,					// (int) number of integrand evaluations per interations (to start)
		0,						// (int) increase in number of integrand evaluations per interations
		1000,					// (int) batch size for sampling
		0,						// (int) grid number, 1-10 => up to 10 grids can be stored, and re-used for other integrands (provided they are not too different)
		"", 					// (char*) name of state file => state can be stored and retrieved for further refinement
		NULL,					// (int*) "spinning cores": -1 || NULL <=> integrator takes care of starting & stopping child processes (other value => keep or retrieve child processes, probably not useful here)
		&neval,					// (int*) actual number of evaluations done
		&nfail,					// 0=desired accuracy was reached; -1=dimensions out of range; >0=accuracy was not reached
		&mcResult,				// (double*) integration result ([ncomp])
		&mcError,				// (double*) integration error ([ncomp])
		&prob					// (double*) Chi-square p-value that error is not reliable (ie should be <0.95) ([ncomp])
	);
	
	cout << "Integration done." << endl;

	/*for(Long_t loop=0; loop<nPoints; loop++){
		int count_old = mycount;
		FoamX->MakeEvent();					// generate MC event
		FoamX->GetMCvect( MCvect );	 // get generated vector (x,y)
	
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
	}*/

	*time = chrono.CpuTime();

	cout << "CPU time : " << chrono.CpuTime() << "	Real-time : " << chrono.RealTime() << endl;
	cout << "Status: " << nfail << ", nr. fails: " << mycount << endl;
	mycount = 0;

	cout << " mcResult= " << mcResult << " +- " << mcError << " in " << neval << " evaluations. Chi-square prob. = " << prob << endl;
	
	/*TCanvas *c = new TCanvas("c","Canvas for plotting",600,600);
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
	delete c; c = 0;*/

	*error = mcError;
	if(std::isnan(*error))
	*error = 0.;
	if(std::isnan(mcResult))
	mcResult = 0.;
	return mcResult;
}


int main(int argc, char *argv[])
{
	TString inputFile(argv[1]);
	TString outputFile(argv[2]);
	int start_evt = atoi(argv[3]);
	int end_evt = atoi(argv[4]);

	//gSystem->Load("libDelphes");

	// Create chain of root trees
	TChain chain("Delphes");
	chain.Add(inputFile);

	// Create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

	// Get pointers to branches used in this analysis
	TClonesArray *branchGen = treeReader->UseBranch("Particle");

	cout << "Entries:" << chain.GetEntries() << endl;

	TFile* outFile = new TFile(outputFile, "RECREATE");
	TTree* outTree = chain.CloneTree(0);

	double Weight_TT_cpp, Weight_TT_Error_cpp;
	bool Weighted_TT_cpp;
	double time = 0;
	outTree->Branch("Weight_TT_cpp", &Weight_TT_cpp);
	outTree->Branch("Weight_TT_Error_cpp", &Weight_TT_Error_cpp);
	outTree->Branch("Weighted_TT_cpp", &Weighted_TT_cpp);
	outTree->Branch("Weight_TT_cpp_time", &time);

	//ofstream fout(outputFile);

	// Full weights
	//double madweight1[10] = {1.58495292058e-21, 2.09681384879e-21, 4.34399623629e-22, 1.68163897955e-22, 3.20350498956e-22, 5.22232034307e-22, 6.04738375743e-21, 9.55643564854e-22, 8.12425265344e-22, 5.81210532053e-23};
	//double madweight2[10] = {1.02514966131e-21, 1.45375719248e-21, 1.65080839221e-22, 1.55653414654e-24, 5.60531044001e-25, 1., 9.70526105314e-24, 3.89103636371e-22, 6.38206925825e-23, 9.37189585544e-26};
	/*double madweight1[10] = {1.48990458909e-21,2.00433822978e-21,4.08998078881e-22,1.56339237714e-22,2.98606743727e-22,4.79498317117e-22,5.63645701583e-21,8.99177777775e-22,7.68316733666e-22,5.42606461617e-23};
	double madweight1Err[10] = {8.63813589113e-24,1.08426062115e-23,2.5750146827e-24,7.0506407196e-25,1.10554655068e-24,2.31140842678e-24,2.71677566322e-23,4.8290429288e-24,1.69718762833e-24,2.66346844676e-25};
	double madweight2[10] = {9.62646303545e-22,1.38143123163e-21,1.54526017444e-22,1.45628835295e-24,6.80263123625e-25,0.,1.07797730384e-23,3.61278172744e-22,6.19087950579e-23,7.20276231557e-26};
	double madweight2Err[10] = {2.96180414077e-24,4.8856162625e-24,1.0218999515e-24,1.29754825587e-25,2.72733072519e-25,0.,4.03010515215e-24,4.29592188061e-24,1.67765665953e-24,8.06569780018e-27};*/

	// NWA weights
	//double madweight1[10] = {1.26069381322e-21, 2.85437676736e-21, 4.8136572269e-22, 1., 3.99894656854e-22, 5.7603822256e-22, 6.99323258475e-21, 1.0892124248e-21, 8.28291668972e-22, 1.};
	//double madweight2[10] = {1.46272073513e-21, 1.51733772927e-21, 1.61193875253e-22, 1., 1., 1., 1., 1., 1., 1.};
	// NWA weights, first sol
	//double madweight1[3] = {8.93501179418e-22, 7.42359826601e-22, 1.49577468649e-22};
	//double madweight2[3] = {1.04113131882e-23, 7.04643552065e-22, 4.3214935529e-23};
	// NWA weights, first sol, ME=1
	//double madweight1[3] = {1.9968889994e-17, 1.10734832869e-17, 2.17966664756e-18};
	//double madweight2[3] = {1.2718723458e-19, 2.38734853175e-17, 6.27800021816e-19};

	if(end_evt >= chain.GetEntries())
		end_evt = chain.GetEntries()-1;

	for(int entry = start_evt; entry <= end_evt ; ++entry){
		
		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);
		chain.GetEntry(entry);

		TLorentzVector gen_ep, gen_mum, gen_Met, gen_nm, gen_ne;

		GenParticle *gen;

		//int count_ep=0, count_mum=0;

		for (Int_t i = 0; i < branchGen->GetEntries(); i++){
			gen = (GenParticle *) branchGen->At(i);
			//cout << "Status=" << gen->Status << ", PID=" << gen->PID << ", E=" << gen->P4().E() << endl;
			if (gen->Status == 3){
				if (gen->PID == -11){
					gen_ep = gen->P4();
					//count_ep++;
				}else if (gen->PID == 13){
					gen_mum = gen->P4();
					//count_mum++;
				}
				else if (gen->PID == 12) {gen_Met += gen->P4(); gen_ne = gen->P4();}
				else if (gen->PID == -14) {gen_Met += gen->P4(); gen_nm = gen->P4();}
			}
		}

		cout << "p2x : " << gen_ne.Px() << " , p2y : " << gen_ne.Py() << " , p2z : " << gen_ne.Pz() << " , E2 : " << gen_ne.E() << endl;
		cout << "p2x : " << gen_nm.Px() << " , p2y : " << gen_nm.Py() << " , p2z : " << gen_nm.Pz() << " , E2 : " << gen_nm.E() << endl;

		double Pzext = (gen_Met+gen_ep+gen_mum).Pz();
		double Eext  = (gen_Met+gen_ep+gen_mum).E();

		cout << "pz   : " << (gen_Met+gen_ep+gen_mum).Pz() << endl;
		cout << "etot : " << (gen_Met+gen_ep+gen_mum).E() << endl;
		cout << "q1 : " << (Pzext+Eext)/(2*4000) << endl;;
		cout << "q2 : " << (Pzext-Eext)/(2*4000) << endl;;
	


		//if(count_ep != 1 || count_mum != 1)
		//	continue;
		//gen_Met.SetPz(0.);
	
		cout << "From MadGraph:" << endl;
		cout << "Electron" << endl;
		cout << gen_ep.E() << "," << gen_ep.Px() << "," << gen_ep.Py() << "," << gen_ep.Pz() << endl;
		cout << "Muon" << endl;
		cout << gen_mum.E() << "," << gen_mum.Px() << "," << gen_mum.Py() << "," << gen_mum.Pz() << endl;
		cout << "MET" << endl;
		cout << gen_Met.E() << "," << gen_Met.Px() << "," << gen_Met.Py() << "," << gen_Met.Pz() << endl;

		Weight_TT_cpp = 0.;
		Weight_TT_Error_cpp = 0.;

		for(int permutation = 1; permutation <= 2; permutation ++){

			//entry = 3; permutation = 2;

			count_perm = permutation;
		
			double weight = 0, temp_time;

			vector<double> nbr_points;
			vector<double> wgts;
			vector<double> wgtsErr;
			vector<double> wgtsErrX;
			/*double max_wgt = 0.;
			double min_wgt = 1.;*/


			for(int k = 1; k <= 51; k+=5){
				/*int nCells = k*10;
				int nSampl = 100;
				int nPoints = 50000;*/
				/*int nCells = k*40;
				int nSampl = 200;
				int nPoints = 50000;*/
				/*int nCells = 750;
				int nSampl = 50;
				int nPoints = 50000;*/

				double error = 0;

				if(permutation == 1)
					weight = ME(&error, gen_ep, gen_mum, gen_Met, &temp_time)/2;
				if(permutation == 2)
					weight = ME(&error, gen_ep, gen_mum, gen_Met, &temp_time)/2;

				Weight_TT_cpp += weight;
				Weight_TT_Error_cpp += pow(error/2,2.);
				time += temp_time;
				
				/*nbr_points.push_back(nCells);
				//nbr_points.push_back(nPoints);
				wgts.push_back(weight);
				wgtsErr.push_back(error);
				wgtsErrX.push_back(0.);
				if(weight > max_wgt) max_wgt = weight;
				if(weight < min_wgt) min_wgt = weight;
				
				if(abs(weight-madweight1[entry]) > abs(madweight2[entry]-weight)){
					double weight3 = madweight2[entry];
					madweight2[entry] = madweight1[entry];
					madweight1[entry] = weight3;
					
					weight3 = madweight2Err[entry];
					madweight2Err[entry] = madweight1Err[entry];
					madweight1Err[entry] = weight3;
				}
				
				if(madweight1[entry] == 0.)
					fout << entry << "/" << permutation << ", points=" << nPoints << ", cells=" << nCells << ", evts/cell=" << nSampl << ": weight = " << weight << " +- " << error << " (1)" << endl;
				else{
					double ratErr;
					if(weight != 0)
						ratErr = weight/madweight1[entry] * TMath::Sqrt( pow(error/weight,2.) + pow(madweight1Err[entry]/madweight1[entry],2.) );
					else
						ratErr = 0;
					fout << entry << "/" << permutation << ", points=" << nPoints << ", cells=" << nCells << ", evts/cell=" << nSampl << ": weight = " << weight << " +- " << error << " (" << weight / (double) madweight1[entry] << " +- " << ratErr << ")" << endl;
				}*/

				break;

			}

			/*double mwX[2] = {nbr_points.at(0), nbr_points.at(nbr_points.size()-1)};
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
			//c->Print(TString("plots/test3_")+SSTR(entry)+"_"+SSTR(permutation)+".png");
			c->Print(TString("plots/wgt_vs_cells_50000p_500c_100s_")+SSTR(entry)+"_"+SSTR(permutation)+".png");
			//c->Print(TString("plots/wgt_vs_points_102000p_5100c_50s_")+SSTR(entry)+"_"+SSTR(permutation)+".png");
			delete wgt_vs_points; wgt_vs_points = 0;
			delete madwgt_vs_points; madwgt_vs_points = 0;
			delete c;*/
		
		//break;
		}
		
		Weight_TT_Error_cpp = TMath::Sqrt(Weight_TT_Error_cpp);
		Weighted_TT_cpp = true;

		outTree->Fill();

		count_wgt++;
		//fout << endl;
		//break;
	}

	outTree->Write();
	outFile->Close();

	delete treeReader; 
}

