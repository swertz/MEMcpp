#include <vector>
#include <iostream>
#include <complex>

#include "utils.h"

using namespace std;

typedef complex<double> cd;

bool solveQuadratic(const double a, const double b, const double c, vector<double>& roots){
	if(a == 0){
		roots.push_back(-c/b);
		return true;
	}

	const double rho = pow(b,2.) - 4*a*c;

	if(rho >= 0.){
		const double x1 = (-b + sqrt(rho))/(2*a);
		const double x2 = (-b - sqrt(rho))/(2*a);

		roots.push_back(x1);
		roots.push_back(x2);

		return true;
	}else
		return false;
}

bool solveCubic(const double a, const double b, const double c, const double d, vector<double>& roots){
	if(a == 0)
		return solveQuadratic(b, c, d, roots);

	const cd min1 = -1.;
	const cd im = sqrt(min1);

	const cd u1 = 1.;
	const cd u2 = (-1. + im*sqrt(3.))/2.;
	const cd u3 = (-1. - im*sqrt(3.))/2.;

	const double delta0 = pow(b,2.) - 3*a*c;
	const double delta1 = 2*pow(b,3.) - 9*a*b*c + 27*d*pow(a,2.);

	const cd C = pow((delta1 + sqrt(pow(delta1,2.) - 4*pow(delta0,3.)))/2., 1./3.);

	const cd x1 = -1./(3*a) * ( b + u1*C + delta0/(u1*C) );
	const cd x2 = -1./(3*a) * ( b + u2*C + delta0/(u2*C) );
	const cd x3 = -1./(3*a) * ( b + u3*C + delta0/(u3*C) );
	
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

bool solveQuartic(const double a, const double b, const double c, const double d, const double e, vector<double>& roots){
	// see http://en.wikipedia.org/wiki/Quartic_function#Solving_a_quartic_equation
	
	if(a == 0.)
		return solveCubic(b, c, d, e, roots);
	
	/*double delta = 256.*pow(a,3.)*pow(e,3.) - 192.*pow(a,2.)*b*d*pow(e,2.) - 128.*pow(a*c*e,2.)
		+ 144.*c*e*pow(a*d,2.) - 27.*pow(a,2.)*pow(d,4.) + 144.*a*c*pow(b*e,2.) - 6.*a*e*pow(b*d,2.)
		- 80.*a*b*d*e*pow(c,2.) + 18.*a*b*c*pow(d,3.) + 16*a*e*pow(c,4.) - 4.a*pow(c,3.)*pow(d,2.)
		- 27.*pow(b,4.)*pow(e,2.) + 18*c*d*e*pow(b,3.) - 4.*pow(d*b,3.) - 4.*e*pow(b,2.)*pow(c,3.)
		+ pow(b*c*d,2.);*/
	
	const cd delta0 = pow(c,2.) - 3.*b*d + 12.*a*e;
	const cd delta1 = 2.*pow(c,3.) - 9.*b*c*d + 27.*pow(b,2)*e + 27.*a*pow(d,2.) - 72.*a*c*e;
	const cd p = (8.*a*c - 3.*pow(b,2.))/(8.*pow(a,2.));
	const cd q = (pow(b,3.) - 4.*a*b*c + 8.*pow(a,2.)*d)/(8*pow(a,3.));

	const cd Q = pow( ( delta1 + sqrt( pow(delta1,2.) - 4.*pow(delta0,3.) ) )/2., 1./3.);
	const cd S = 0.5 * sqrt( -2./3. * p + 1./(3.*a) * ( Q + delta0/Q ) );

	// four solutions
	
	const cd x1 = -b/(4.*a) - S + 0.5 * sqrt( -4.*pow(S,2.) - 2.*p + q/S );
	const cd x2 = -b/(4.*a) - S - 0.5 * sqrt( -4.*pow(S,2.) - 2.*p + q/S );

	const cd x3 = -b/(4.*a) + S + 0.5 * sqrt( -4.*pow(S,2.) - 2.*p - q/S );
	const cd x4 = -b/(4.*a) + S - 0.5 * sqrt( -4.*pow(S,2.) - 2.*p - q/S );

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

bool solve2Quads(const double a11, const double a22, const double a12, const double a10, const double a01, const double a00,
				const double b11, const double b22, const double b12, const double b10, const double b01, const double b00,
				vector<double>& E1, vector<double>& E2){
	const double alpha = b11*a22-a11*b22;
	const double beta = b11*a12-a11*b12;
	const double gamma = b11*a10-a11*b10;
	const double delta = b11*a01-a11*b01;
	const double omega = b11*a00-a11*b00;

	const double a = a11*pow(alpha,2.) + a22*pow(beta,2.) - a12*alpha*beta;
	const double b = 2.*a11*alpha*delta - a12*alpha*gamma - a12*delta*beta - a10*alpha*beta + 2.*a22*beta*gamma + a01*pow(beta,2.);
	const double c = a11*pow(delta,2.) + 2.*a11*alpha*omega - a12*delta*gamma - a12*omega*beta - a10*alpha*gamma - a10*delta*beta
		+ a22*pow(gamma,2.) + 2.*a01*beta*gamma + a00*pow(beta,2.);
	const double d = 2.*a11*delta*omega - a12*omega*gamma - a10*delta*gamma - a10*omega*beta + a01*pow(gamma,2.) + 2.*a00*beta*gamma;
	const double e = a11*pow(omega,2.) - a10*omega*gamma + a00*pow(gamma,2.);

	solveQuartic(a, b, c, d, e, E2);

	for(unsigned int i =0; i < E2.size(); ++i){
		const double e2 = E2.at(i);
		const double e1 = -(alpha * pow(e2,2.) + delta*e2 + omega)/(beta*e2 + gamma);
		E1.push_back(e1);
	}

	return true;
}

double BreitWigner(const double s, const double m, const double g){
	/*double ga = sqrt(m*m*(m*l+g*g));
	double k = 2*sqrt(2)*m*g*ga/(TMath::Pi()*sqrt(m*m+ga));*/
	double k = m*g;
	
	//cout << "BW(" << s << "," << m << "," << g << ")=" << k/(pow(s-m*m,2.) + pow(m*g,2.)) << endl;
	
	return k/(pow(s-m*m,2.) + pow(m*g,2.));
}

