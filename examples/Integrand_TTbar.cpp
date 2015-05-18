#include <iostream>
#include <vector>

#include "TLorentzVector.h"
#include "TMath.h"

#include "MEWeight.h"
#include "jacobianD.h"

#define M_T 173.
#define G_T 1.4915

#define M_W 80.419
#define G_W 2.0476

#define SQRT_S 13000

using namespace std;

double MEWeight::Integrand(const double* Xarg, const double *weight){
  double Value = 0.;

  for(int i=0; i<4; ++i){
    if(Xarg[i] == 1.){
      //mycount++;
      return 0;
    }
  }

  TLorentzVector p3rec = myEvent->GetP3();
  TLorentzVector p4rec = myEvent->GetP4();
  TLorentzVector p5rec = myEvent->GetP5();
  TLorentzVector p6rec = myEvent->GetP6();
  TLorentzVector Met = myEvent->GetMet();

  double TFValue = 1.;

  TLorentzVector p3 = p3rec;
  const double p3DeltaRange = myTF->GetDeltaRange("electron");
  const double E3gen = p3rec.E() - myTF->GetDeltaMax("electron") + p3DeltaRange * Xarg[4];
  p3.SetE(E3gen);
  if(p3DeltaRange != 0.)
    TFValue *= myTF->Evaluate("electron", p3rec.E(), E3gen) * p3DeltaRange;

  TLorentzVector p4 = p4rec;
  const double p4DeltaRange = myTF->GetDeltaRange("jet");
  const double E4gen = p4rec.E() - myTF->GetDeltaMax("jet") + p4DeltaRange * Xarg[5];
  p4.SetE(E4gen);
  if(p4DeltaRange != 0.)
    TFValue *= myTF->Evaluate("jet", p4rec.E(), E4gen) * p4DeltaRange;

  TLorentzVector p5 = p5rec;
  const double p5DeltaRange = myTF->GetDeltaRange("muon");
  const double E5gen = p5rec.E() - myTF->GetDeltaMax("muon") + p5DeltaRange * Xarg[6];
  p5.SetE(E5gen);
  if(p5DeltaRange != 0.)
    TFValue *= myTF->Evaluate("muon", p5rec.E(), E5gen) * p5DeltaRange;

  TLorentzVector p6 = p6rec;
  const double p6DeltaRange = myTF->GetDeltaRange("jet");
  const double E6gen = p6rec.E() - myTF->GetDeltaMax("jet") + p6DeltaRange * Xarg[7];
  p6.SetE(E6gen);
  if(p6DeltaRange != 0.)
    TFValue *= myTF->Evaluate("jet", p6rec.E(), E6gen) * p6DeltaRange;

  //cout << "Final TF = " << TFValue << endl;


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

  if(s13 > s134 || s25 > s256 || s13 < p3.M() || s25 < p5.M() || s134 < p4.M() || s256 < p6.M()){
    //cout << "Masses too small!" << endl;
    //mycount++;
    return 0;
  }

  //cout << "weight = " << *weight << endl;
  
  std::vector<TLorentzVector> p1vec, p2vec;

  ComputeTransformD(s13, s134, s25, s256,
                    p3, p4, p5, p6, Met,
                    p1vec, p2vec);

  int countSol = 0;

  for(unsigned short i = 0; i < p1vec.size(); ++i){

    const TLorentzVector p1 = p1vec.at(i);
    const TLorentzVector p2 = p2vec.at(i);

    const TLorentzVector p13 = p1 + p3;
    const TLorentzVector p134 = p1 + p3 + p4;
    const TLorentzVector p25 = p2 + p5;
    const TLorentzVector p256 = p2 + p5 + p6;

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
      //mycount++;
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
      //mycount++;
      continue;
      //break;
    }
  
    // Compute the Pdfs
    const double x1 = TMath::Abs(q1Pz/(SQRT_S/2.));
    const double x2 = TMath::Abs(q2Pz/(SQRT_S/2.));

    const double pdf1_1 = ComputePdf(21, x1, pow(M_T,2));
    const double pdf1_2 = ComputePdf(21, x2, pow(M_T,2));
  
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
    setProcessMomenta(p);

    // Evaluate matrix element
    computeMatrixElements();
    const double* const matrix_elements1 = getMatrixElements();
    
    // free up memory
    for(unsigned int j = 0; j < p.size(); ++j){
      delete [] p.at(j); p.at(j) = NULL;
    }

    const double thisSolResult = PhaseSpaceIn * matrix_elements1[0] * pdf1_1 * pdf1_2 * PhaseSpaceOut * jac * flatterJac * TFValue;
    Value += thisSolResult; 
    
    // Check whether the next solutions for the neutrinos are the same => don't redo all this!
    int countEqualSol = 1;
    for(unsigned int j = i+1; j<p1vec.size(); j++){
      if(p1 == p1vec.at(j) && p2 == p2vec.at(j)){
        Value += thisSolResult;
        countEqualSol++;
      }
    }

    //cout << "Found PDF1 = " << pdf1_1 << ", PDF2 = " << pdf1_2 << ", PS in = " << PhaseSpaceIn << ", PS out = " << PhaseSpaceOut << ", jac = " << jac << endl;
    //cout << "===> Matrix element = " << matrix_elements1[0] << ", prod = " << thisSolResult << ", multiplicity = " << countEqualSol << endl << endl; 

    // If we have included the next solutions already, skip them!
    i += countEqualSol - 1;
    countSol += countEqualSol;
    
    //myEvent->GetTTbar()->Fill(tot.M(), *weight * (double)countEqualSol * thisSolResult);
  }

  if(Value == 0.){
    //mycount++;
    //cout << "Zero!" << endl;
    //return 0;
  }

  // ### FOR NWA
  //double flatterJac = pow(TMath::Pi(),4.) * (M_W*G_W * M_T*G_T * M_W*G_W * M_T*G_T);

  //cout << "## Phase Space point done. Integrand = " << integrand << ", flatterjac = " << flatterJac << ", prod = " << integrand*flatterJac <<  endl;

  return Value;
}