#include <iostream>
#include <vector>

#include "TLorentzVector.h"
#include "TMath.h"

#include "utils.h"
#include "jacobianD.h"

using namespace std;

int ComputeTransformD(const double s13, const double s134, const double s25, const double s256,
                      const TLorentzVector p3, const TLorentzVector p4, const TLorentzVector p5, const TLorentzVector p6, const TLorentzVector Met,
                      std::vector<TLorentzVector> &p1, std::vector<TLorentzVector> &p2){
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

  // For each solution (E1,E2), find the neutrino 4-momenta p1,p2

  if(E1.size() == 0){
    //cout << "No solutions!" << endl;
    return 0;
  }

  for(unsigned int i=0; i<E1.size(); i++){
    const double e1 = E1.at(i);
    const double e2 = E2.at(i);

    //cout << endl << "## Evaluating Matrix Element based on solutions e1 = " << e1 << ", e2 = " << e2 << endl << endl;

    if(e1 < 0. || e2 < 0.){
      //cout << "Neg energies." << endl;
      continue;
    }

    TLorentzVector tempp1,tempp2;

    tempp1.SetPx( alpha1*e1 + beta1*e2 + gamma1 );
    tempp1.SetPy( alpha2*e1 + beta2*e2 + gamma2 );
    tempp1.SetPz( alpha3*e1 + beta3*e2 + gamma3 );
    tempp1.SetE(e1);

    tempp2.SetPx( alpha5*e1 + beta5*e2 + gamma5 );
    tempp2.SetPy( alpha6*e1 + beta6*e2 + gamma6 );
    tempp2.SetPz( alpha4*e1 + beta4*e2 + gamma4 ); 
    tempp2.SetE(e2);

    p1.push_back(tempp1);
    p2.push_back(tempp2);
  }

  return p1.size();
}

double computeJacobianD(const std::vector<TLorentzVector> p, const double sqrt_s){
	
	const double E1  = p.at(0).E();
	const double p1x = p.at(0).Px();
	const double p1y = p.at(0).Py();
	const double p1z = p.at(0).Pz();

	const double E2  = p.at(1).E();
	const double p2x = p.at(1).Px();
	const double p2y = p.at(1).Py();
	const double p2z = p.at(1).Pz();

	const double E3  = p.at(2).E();
	const double p3x = p.at(2).Px();
	const double p3y = p.at(2).Py();
	const double p3z = p.at(2).Pz();

	const double E4  = p.at(3).E();
	const double p4x = p.at(3).Px();
	const double p4y = p.at(3).Py();
	const double p4z = p.at(3).Pz();

	const double E5  = p.at(4).E();
	const double p5x = p.at(4).Px();
	const double p5y = p.at(4).Py();
	const double p5z = p.at(4).Pz();

	const double E6  = p.at(5).E();
	const double p6x = p.at(5).Px();
	const double p6y = p.at(5).Py();
	const double p6z = p.at(5).Pz();

	const double E34  = E3 + E4;
	const double p34x = p3x + p4x;
	const double p34y = p3y + p4y;
	const double p34z = p3z + p4z;

	const double E56  = E5 + E6;
	const double p56x = p5x + p6x;
	const double p56y = p5y + p6y;
	const double p56z = p5z + p6z;

    double inv_jac = E3*(E5*
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
                  
    inv_jac *= 8.*16.*pow(TMath::Pi()*sqrt_s,2.);

	//std::cout << "jac=" << abs(jac) << std::endl;
	
	if(TMath::Abs(inv_jac) < INV_JAC_MIN){
        std::cout << "Warning: jacobian is close to zero!" << std::endl;
		return -1.;
	}else
		return 1./TMath::Abs(inv_jac);
}

