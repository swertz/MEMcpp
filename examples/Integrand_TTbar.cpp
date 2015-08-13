#include <iostream>
#include <vector>
#include <cmath>

#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/Boost.h"
#include "TMath.h"

#include "MEWeight.h"
#include "jacobianD.h"
#include "utils.h"

#define M_T 173.
#define G_T 1.4915

#define M_W 80.419
#define G_W 2.0476

#define SQRT_S 13000

using namespace std;

double MEWeight::Integrand(const double* Xarg, const double *weight){
  double returnValue = 0.;

  for(int i=0; i<4; ++i){
    if(Xarg[i] == 1.)
      return 0;
  }

  const ROOT::Math::PtEtaPhiEVector p3rec( _recEvent->GetP3() );
  const ROOT::Math::PtEtaPhiEVector p4rec( _recEvent->GetP4() );
  const ROOT::Math::PtEtaPhiEVector p5rec( _recEvent->GetP5() );
  const ROOT::Math::PtEtaPhiEVector p6rec( _recEvent->GetP6() );
  const ROOT::Math::PtEtaPhiEVector RecMet( _recEvent->GetMet() );

  // Define ISR vector from the observed particles and MET
  // Do this in PxPyPzE basis as this will be more practical for the next stages (transverse boost)
  ROOT::Math::PxPyPzEVector ISR( -(p3rec + p4rec + p5rec + p6rec + RecMet) );

  ///// Transfer functions

  double TFreturnValue = 1.;

  ROOT::Math::PtEtaPhiEVector p3gen = p3rec;
  const double E3rec = p3rec.E();
  const double p3DeltaRange = _TF->GetDeltaRange("electron", E3rec);
  const double E3gen = E3rec - _TF->GetDeltaMax("electron", E3rec) + p3DeltaRange * Xarg[4];
  const double pt3gen = sqrt( SQ(E3gen) - SQ(p3rec.M()) ) / cosh(p3rec.Eta());
  p3gen.SetCoordinates(pt3gen, p3rec.Eta(), p3rec.Phi(), E3gen);
  if(p3DeltaRange != 0.)
    TFreturnValue *= _TF->Evaluate("electron", E3rec, E3gen) * p3DeltaRange * dEoverdP(E3gen, p3gen.M());

  ROOT::Math::PtEtaPhiEVector p5gen = p5rec;
  const double E5rec = p5rec.E();
  const double p5DeltaRange = _TF->GetDeltaRange("muon", E5rec);
  const double E5gen = E5rec - _TF->GetDeltaMax("muon", E5rec) + p5DeltaRange * Xarg[6];
  const double pt5gen = sqrt( SQ(E5gen) - SQ(p5rec.M()) ) / cosh(p5rec.Eta());
  p5gen.SetCoordinates(pt5gen, p5rec.Eta(), p5rec.Phi(), E5gen);
  if(p5DeltaRange != 0.)
    TFreturnValue *= _TF->Evaluate("muon", E5rec, E5gen) * p5DeltaRange * dEoverdP(E5gen, p5gen.M());

  ROOT::Math::PtEtaPhiEVector p4gen = p4rec;
  const double E4rec = p4rec.E();
  const double p4DeltaRange = _TF->GetDeltaRange("jet", E4rec);
  const double E4gen = E4rec - _TF->GetDeltaMax("jet", E4rec) + p4DeltaRange * Xarg[5];
  const double pt4gen = sqrt( SQ(E4gen) - SQ(p4rec.M()) ) / cosh(p4rec.Eta());
  p4gen.SetCoordinates(pt4gen, p4rec.Eta(), p4rec.Phi(), E4gen);
  if(p4DeltaRange != 0.)
    TFreturnValue *= _TF->Evaluate("jet", E4rec, E4gen) * p4DeltaRange * dEoverdP(E4gen, p4gen.M());

  ROOT::Math::PtEtaPhiEVector p6gen = p6rec;
  const double E6rec = p6rec.E();
  const double p6DeltaRange = _TF->GetDeltaRange("jet", E6rec);
  const double E6gen = p6rec.E() - _TF->GetDeltaMax("jet", E6rec) + p6DeltaRange * Xarg[7];
  const double pt6gen = sqrt( SQ(E6gen) - SQ(p6rec.M()) ) / cosh(p6rec.Eta());
  p6gen.SetCoordinates(pt6gen, p6rec.Eta(), p6rec.Phi(), E6gen);
  if(p6DeltaRange != 0.)
    TFreturnValue *= _TF->Evaluate("jet", E6rec, E6gen) * p6DeltaRange * dEoverdP(E6gen, p6gen.M());

  // In the following, we want to use PxPyPzE vectors, since the change of variables is done over those variables => we are already in the right basis, no need to recompute quantities every time
  ROOT::Math::PxPyPzEVector p3(p3gen);
  ROOT::Math::PxPyPzEVector p4(p4gen);
  ROOT::Math::PxPyPzEVector p5(p5gen);
  ROOT::Math::PxPyPzEVector p6(p6gen);
  ROOT::Math::PxPyPzEVector Met(RecMet);

  //cout << "Final TF = " << TFreturnValue << endl;

  // We flatten the Breit-Wigners by doing a change of variable for each resonance separately
  // s = M G tan(y) + M^2
  // jac = M G / cos^2(y)
  // ==> BW(s(y))*jac(y) is flat in the variable y, as BW(s) = 1/((s-M^2)^2 - (GM)^2)
  // Where y = -arctan(M/G) + (pi/2+arctan(M/G))*x_foam (x_foam between 0 and 1 => s between 0 and infinity)

  const double range1 = TMath::Pi()/2. + atan(M_W/G_W);
  const double y1 = - atan(M_W/G_W) + range1 * Xarg[0];
  const double s13 = M_W * G_W * tan(y1) + SQ(M_W);

  //cout << "y1=" << y1 << ", m13=" << sqrt(s13) << endl;

  const double range2 = TMath::Pi()/2. + atan(M_T/G_T);
  const double y2 = - atan(M_T/G_T) + range2 * Xarg[1];
  const double s134 = M_T * G_T * tan(y2) + SQ(M_T);

  //cout << "y2=" << y2 << ", m134=" << sqrt(s134) << endl;

  const double range3 = TMath::Pi()/2. + atan(M_W/G_W);
  const double y3 = - atan(M_W/G_W) + range3 * Xarg[2];
  const double s25 = M_W * G_W * tan(y3) + SQ(M_W);

  //cout << "y3=" << y3 << ", m25=" << sqrt(s25) << endl;

  const double range4 = TMath::Pi()/2. + atan(M_T/G_T);
  const double y4 = - atan(M_T/G_T) + range4 * Xarg[3];
  const double s256 = M_T * G_T * tan(y4) + SQ(M_T);

  //cout << "y4=" << y4 << ", m256=" << sqrt(s256) << endl;
  
  double flatterJac = range1 * range2 * range3 * range4;
  flatterJac *= M_W*G_W * M_T*G_T * M_W*G_W * M_T*G_T;
  flatterJac /= pow(cos(y1) * cos(y2) * cos(y3) * cos(y4), 2.);

  if(s13 > s134 || s25 > s256 || s13 < p3.M() || s25 < p5.M() || s134 < p4.M() || s256 < p6.M())
    return 0;

  //cout << "weight = " << *weight << endl;
  
  std::vector<ROOT::Math::PxPyPzEVector> p1vec, p2vec;

  ComputeTransformD(s13, s134, s25, s256,
                    p3, p4, p5, p6, Met, ISR,
                    p1vec, p2vec);

  int countSol = 0;

  for(unsigned short i = 0; i < p1vec.size(); ++i){

    const ROOT::Math::PxPyPzEVector p1 = p1vec.at(i);
    const ROOT::Math::PxPyPzEVector p2 = p2vec.at(i);

    /*const ROOT::Math::PxPyPzEVector p13 = p1 + p3;
    const ROOT::Math::PxPyPzEVector p134 = p1 + p3 + p4;
    const ROOT::Math::PxPyPzEVector p25 = p2 + p5;
    const ROOT::Math::PxPyPzEVector p256 = p2 + p5 + p6;

    cout << "Solution " << i << ":" << endl;
    cout << "Input: W+ mass=" << sqrt(s13) << ", Top mass=" << sqrt(s134) << ", W- mass=" << sqrt(s25) << ", Anti-top mass=" << sqrt(s256) << endl;
    cout << "Output: W+ mass=" << p13.M() << ", Top mass=" << p134.M() << ", W- mass=" << p25.M() << ", Anti-top mass=" << p256.M() << endl << endl;
    //cout << "Differences: W+ mass=" << sqrt(s13)-p13.M() << ", Top mass=" << sqrt(s134)-p134.M() << ", W- mass=" << sqrt(s25)-p25.M() << ", Anti-top mass=" << sqrt(s256)-p256.M() << endl << endl;
    cout << "Accuracies: " << endl;
    cout << "W+ mass=" << (sqrt(s13)-p13.M())/p13.M() << endl;
    cout << "Top mass=" << (sqrt(s134)-p134.M())/p134.M() << endl;
    cout << "W- mass=" << (sqrt(s25)-p25.M())/p25.M() << endl;
    cout << "Anti-top mass=" << (sqrt(s256)-p256.M())/p256.M() << endl << endl;*/
    
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
  
    const ROOT::Math::PxPyPzEVector tot = p1 + p2 + p3 + p4 + p5 + p6;

    /*cout << "**********************" << endl;
    cout << "Total: " << tot  << endl;
    cout << "ISR: " << ISR << endl;*/
    
    const double ETot = tot.E();
    const double PzTot = tot.Pz();

    const double q1Pz = (PzTot + ETot)/2.;
    const double q2Pz = (PzTot - ETot)/2.;

    //cout << "===> Eext=" << ETot << ", Pzext=" << PzTot << ", q1Pz=" << q1Pz << ", q2Pz=" << q2Pz << endl << endl;
  
    if(q1Pz > SQRT_S/2. || q2Pz < -SQRT_S/2. || q1Pz < 0. || q2Pz > 0.)
      continue;
    
    // Compute jacobian from change of variable:
    vector<ROOT::Math::PxPyPzEVector> momenta( { p1, p2, p3, p4, p5, p6 } );
    const double jacobian = computeJacobianD(momenta, SQRT_S);
    if(jacobian <= 0.){
      cout << "Jac infinite!" << endl;
      continue;
    }

    // Compute the Pdfs
    const double x1 = abs(q1Pz/(SQRT_S/2.));
    const double x2 = abs(q2Pz/(SQRT_S/2.));

    // Compute flux factor 1/(2*x1*x2*s)
    const double phaseSpaceIn = 1.0 / ( 2. * x1 * x2 * SQ(SQRT_S) ); 

    // Compute phase space density for observed particles (not concerned by the change of variable)
    // dPhi = |P|^2 sin(theta)/(2*E*(2pi)^3)
    const double dPhip3 = pow(p3.P(),2.)*TMath::Sin(p3.Theta())/(2.0*p3.E()*pow(2.*TMath::Pi(),3));
    const double dPhip4 = pow(p4.P(),2.)*TMath::Sin(p4.Theta())/(2.0*p4.E()*pow(2.*TMath::Pi(),3));
    const double dPhip5 = pow(p5.P(),2.)*TMath::Sin(p5.Theta())/(2.0*p5.E()*pow(2.*TMath::Pi(),3));
    const double dPhip6 = pow(p6.P(),2.)*TMath::Sin(p6.Theta())/(2.0*p6.E()*pow(2.*TMath::Pi(),3));
    const double phaseSpaceOut = dPhip5 * dPhip6 * dPhip3 * dPhip4;

    // Boost the initial particle 4-momenta in order to match parton1 + parton2 = - ISR
    ROOT::Math::PxPyPzEVector parton1(0., 0., q1Pz, q1Pz);
    ROOT::Math::PxPyPzEVector parton2(0., 0., q2Pz, abs(q2Pz));
   
    /*cout << "Before:" << endl;
    cout << " Parton1: " << parton1 << endl;
    cout << " Parton2: " << parton2 << endl;*/
 
    // This gives reasonable results (test transverse < 1 GeV)
    ROOT::Math::XYZVector isrBoostVector = -(tot.BoostToCM()); // WTF are -tot.BoostToCM() and -(tot.BoostToCM()) giving different results??
    
    // this does not give the same result as above, since beta_x(boost) = x/E, and while x_ISR = -x_tot, E_ISR != E_tot
    //ROOT::Math::XYZVector isrBoostVector = ISR.BoostToCM(); 
    
    //isrBoostVector.SetZ(0.); // We want a transverse boost only
    //cout << "Boost: " << isrBoostVector << endl;
    ROOT::Math::Boost isrBoost(isrBoostVector);
    parton1 = isrBoost*parton1;
    parton2 = isrBoost*parton2;

    /*cout << "After:" << endl;
    cout << " Parton1: " << parton1 << endl;
    cout << " Parton2: " << parton2 << endl;*/
    
    //ROOT::Math::PxPyPzEVector testT = parton1 + parton2 + ISR;
    //ROOT::Math::PxPyPzEVector testT = -parton1 - parton2 + p1 + p2 + p3 + p4 + p5 + p6;
    //cout << "Test transverse: " << testT << endl;

    // Define initial momenta to be passed to matrix element
    std::vector< std::vector<double> > initialMomenta = 
    {
      { parton1.E(), parton1.Px(), parton1.Py(), parton1.Pz() },
      { parton2.E(), parton2.Px(), parton2.Py(), parton2.Pz() },
    };
    
    // Define final PID and momenta to be passed to matrix element
    std::vector< std::pair<int, std::vector<double> > > finalState = 
    {
      std::make_pair<int, std::vector<double>>(-11, { p3.E(), p3.Px(), p3.Py(), p3.Pz() }),
      std::make_pair<int, std::vector<double>>(12, { p1.E(), p1.Px(), p1.Py(), p1.Pz() }),
      std::make_pair<int, std::vector<double>>(5, { p4.E(), p4.Px(), p4.Py(), p4.Pz() }),
      std::make_pair<int, std::vector<double>>(13, { p5.E(), p5.Px(), p5.Py(), p5.Pz() }),
      std::make_pair<int, std::vector<double>>(-14, { p2.E(), p2.Px(), p2.Py(), p2.Pz() }),
      std::make_pair<int, std::vector<double>>(-5, { p6.E(), p6.Px(), p6.Py(), p6.Pz() }),
    };

    // Evaluate matrix element
    std::map< std::pair<int, int>, double > matrixElements = getMatrixElements(initialMomenta, finalState);

    double thisSolResult = phaseSpaceIn * phaseSpaceOut * jacobian * flatterJac * TFreturnValue;

    double pdfMESum = 0.;
    // If no initial states have been defined explicitly, loop over all states returned by the matrix element
    if(!_initialStates.size()){
      for(auto const &me: matrixElements){
        const double pdf1 = ComputePdf(me.first.first, x1, SQ(M_T));
        const double pdf2 = ComputePdf(me.first.second, x2, SQ(M_T));
        pdfMESum += me.second * pdf1 * pdf2;
        //cout << "Initial state (" << me.first.first << ", " << me.first.second << "): " << me.second << endl;
      }
    }else{
      // Otherwise, loop over all states defined by user
      for(auto const &initialState: _initialStates){
        const double pdf1 = ComputePdf(initialState.first, x1, SQ(M_T));
        const double pdf2 = ComputePdf(initialState.second, x2, SQ(M_T));
        pdfMESum += matrixElements[initialState] * pdf1 * pdf2;
        //cout << "Initial state (" << initialState.first << ", " << initialState.second << "): " << matrixElements[initialState] << endl;
      }
    }

    thisSolResult *= pdfMESum;
    returnValue += thisSolResult; 
    
    // Check whether the next solutions for the neutrinos are the same => don't redo all this!
    int countEqualSol = 1;
    for(unsigned int j = i+1; j<p1vec.size(); j++){
      if(p1 == p1vec.at(j) && p2 == p2vec.at(j)){
        returnValue += thisSolResult;
        countEqualSol++;
      }
    }

    //cout << "Found PDF1 = " << pdf1_1 << ", PDF2 = " << pdf1_2 << ", PS in = " << PhaseSpaceIn << ", PS out = " << PhaseSpaceOut << ", jac = " << jac << endl;
    //cout << "===> Matrix element = " << matrix_elements1[0] << ", prod = " << thisSolResult << ", multiplicity = " << countEqualSol << endl << endl; 

    // If we have included the next solutions already, skip them!
    i += countEqualSol - 1;
    countSol += countEqualSol;
  }

  //cout << "## Phase Space point done. Integrand = " << integrand << ", flatterjac = " << flatterJac << ", prod = " << integrand*flatterJac <<  endl;

  return returnValue;
}
