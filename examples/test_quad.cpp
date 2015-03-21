#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <complex>

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

	cout << "x1 = " << real(x1) << " + i*" << imag(x1) << endl;
	cout << "x2 = " << real(x2) << " + i*" << imag(x2) << endl;
	cout << "x3 = " << real(x3) << " + i*" << imag(x3) << endl;
	cout << "x4 = " << real(x4) << " + i*" << imag(x4) << endl;

	// checking which ones are real, with some tolerance
	
	if( abs(real(x1)) > pow(10,5)*abs(imag(x1)) && abs(real(x1)) != 0 )
		roots.push_back(real(x1));
	
	if( abs(real(x2)) > pow(10,5)*abs(imag(x2)) && abs(real(x2)) != 0 )
		roots.push_back(real(x2));
	
	if( abs(real(x3)) > pow(10,5)*abs(imag(x3)) && abs(real(x3)) != 0 )
		roots.push_back(real(x3));
	
	if( abs(real(x4)) > pow(10,5)*abs(imag(x4)) && abs(real(x4)) != 0 )
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

int main(int argv, char** argc){
	
	vector<double> roots;
	double a=2, b=-1, c=4, d=8, e=-5;
	solveQuartic(a,b,c,d,e,roots);
	for(unsigned int i=0; i<roots.size(); i++){
		double x = roots.at(i);
		double eval = a*pow(x,4)+b*pow(x,3)+c*pow(x,2)+d*x+e;
		cout << "Root " << x << ", eval = " << eval << endl;
	}
	if(roots.size() == 0)
		cout << "No real root found!" << endl;

	cout << endl << endl;

	double a11=-23,a22=14,a12=52,a10=-45,a01=0.1,a00=1.6;
	double b11=55,b22=-898,b12=3,b10=-59,b01=43,b00=0.5;

	vector<double> E1,E2;
	solve2Quads(a11,a22,a12,a10,a01,a00,
				b11,b22,b12,b10,b01,b00,
				E1, E2);
	
	
	for(unsigned int i=0; i<E1.size(); i++){
		double x1 = E1.at(i);
		double x2 = E2.at(i);
		double eval1 = a11*pow(x1,2)+a22*pow(x2,2)+a12*x1*x2+a10*x1+a01*x2+a00;
		double eval2 = b11*pow(x1,2)+b22*pow(x2,2)+b12*x1*x2+b10*x1+b01*x2+b00;
		cout << "Roots " << x1 << ", " << x2 << ", eval1 = " << eval1 << ", eval2 = " << eval2 << endl;
	}
	

	return 0;
}
