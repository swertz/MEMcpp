#include <vector>
#include <iostream>
#include <complex>
#include <cmath>
#include <stdlib.h>

#include "utils.h"

using namespace std;

typedef complex<double> cd;

// Overrides double std::cbrt(double)
// Returns principal cubic root of complex number
// Caution: the principal cubic root of a negative real number is complex
cd cbrt(cd x){ 
    double cbrt_abs = cbrt(abs(x));
    double cbrt_arg = arg(x)/3.;
    return cbrt_abs*cd(cos(cbrt_arg), sin(cbrt_arg));
}

double sign(double x){
    if(x > 0)
        return 1.;
    else if(!x)
        return 0.;
    else
        return -1.;
}

// Compute cos(x +- 2*pi/3) in a more "analytical" way (pm = +- 1)
double cosXpm2PI3(double x, double pm){
    return -0.5*( cos(x) + pm * sin(x) * sqrt(3.) );
}

bool solveQuadratic(const double a, const double b, const double c, vector<double>& roots, bool verbose){
	if(!a){
        if(!b){
            if(verbose)
                cout << "No solution to equation " << a << " x^2 + " << b << " x + " << c << endl << endl;
            return false;
        }
		roots.push_back(-c/b);
        if(verbose)
            cout << "Solution of " << b << " x + " << c << ": " << roots[0] << ", test = " << b*roots[0] + c << endl << endl;
		return true;
	}

	const double rho = SQ(b) - 4.*a*c;

	if(rho >= 0.){
        if(b == 0.){
            roots.push_back( sqrt(rho)/(2.*a) );
            roots.push_back( -sqrt(rho)/(2.*a) );
        }else{
		    const double x = -0.5*(b + sign(b)*sqrt(rho));
		    roots.push_back(x/a);
		    roots.push_back(c/x);
        }
        if(verbose){
            cout << "Solutions of " << a << " x^2 + " << b << " x + " << c << ":" << endl;
            for(unsigned short i=0; i<roots.size(); i++)
	            cout << "x" << i << " = " << roots[i] << ", test = " << a*SQ(roots[i]) + b*roots[i] + c << endl;
            cout << endl;
        }
		return true;
	}else{
        if(verbose)
            cout << "No real solutions to " << a << " x^2 + " << b << " x + " << c << endl << endl;
		return false;
    }
}

bool solveCubic(const double a, const double b, const double c, const double d, vector<double>& roots, bool verbose){
	if(a == 0)
		return solveQuadratic(b, c, d, roots, verbose);

    const double an = b/a;
    const double bn = c/a;
    const double cn = d/a;

    const double Q = SQ(an)/9. - bn/3.;
    const double R = CB(an)/27. - an*bn/6. + cn/2.;

    if( SQ(R) < CB(Q) ){
        const double theta = acos( R/sqrt(CB(Q)) );

        roots.push_back( -2. * sqrt(Q) * cos( theta/3. ) - an/3. );
        roots.push_back( -2. * sqrt(Q) * cosXpm2PI3(theta/3., 1.) - an/3. );
        roots.push_back( -2. * sqrt(Q) * cosXpm2PI3(theta/3., -1.) - an/3. );
    }else{
        const double A = - sign(R) * cbrt( abs(R) + sqrt( SQ(R) - CB(Q) ) );

        double B;

        if(A == 0.)
            B = 0.;
        else
            B = Q/A;

        roots.push_back( A + B - an/3.);
        roots.push_back( A + B - an/3.);
        roots.push_back( A + B - an/3.);
    }
    
    if(verbose){
        cout << "Solutions of " << a << " x^3 + " << b << " x^2 + " << c << " x + " << d << ":" << endl;
        for(unsigned short i=0; i<roots.size(); i++)
	        cout << "x" << i << " = " << roots[i] << ", test = " << a*CB(roots[i]) + b*SQ(roots[i]) + c*roots[i] + d << endl;
        cout << endl;
    }

    return true;

	/*const cd u1 = 1.;
	const cd u2 = (-1. + 1i*sqrt(3.))/2.;
	const cd u3 = (-1. - 1i*sqrt(3.))/2.;

    const cd delta = 18*a*b*c*d - 4*CB(b)*d + SQ(b*c) - 4*a*CB(c) - 27*SQ(a*d);
	const cd delta0 = SQ(b) - 3.*a*c;
	const cd delta1 = 2*CB(b) - 9.*a*b*c + 27.*d*SQ(a);

    cd temp_deltas;
    if(delta != 0.){
        if(delta0 == 0.)
            temp_deltas = delta1;
        else
            temp_deltas = sqrt(SQ(delta1) - 4.*CB(delta0));
    }else{
        if(delta0 == 0.){
            roots.push_back(-b/3.*a);
            roots.push_back(-b/3.*a);
            roots.push_back(-b/3.*a);
            return true;
        }else{
            roots.push_back( (9.*a*d - b*c)/(2.*real(delta0)) );
            roots.push_back( (9.*a*d - b*c)/(2.*real(delta0)) );
            roots.push_back( (4.*a*b*c - 9*SQ(a)*d - CB(b))/(a*real(delta0)) );
            return true;
        }
    }

	const cd C = cbrt( ( delta1 + temp_deltas )/2. );

	cd x1 = -1./(3.*a) * ( b + u1*C + delta0/(u1*C) );
	cd x2 = -1./(3.*a) * ( b + u2*C + delta0/(u2*C) );
	cd x3 = -1./(3.*a) * ( b + u3*C + delta0/(u3*C) );
    
    // reorganize them s.t., if there are double roots, they are grouped (0,1)
    if(real(x1) == real(x3)){
        cd temp = x2;
        x2 = x3;
        x3 = temp;
    }
    if(real(x2) == real(x3)){
        cd temp = x1;
        x1 = x3;
        x3 = temp;
    }*/
    
    /*cout << "Solutions of " << a << " x^3 + " << b << " x^2 + " << c << " x + " << d << endl;
	cout << "x1 = " << real(x1) << "+i*" << imag(x1) << endl;
	cout << "x2 = " << real(x2) << "+i*" << imag(x2) << endl;
	cout << "x3 = " << real(x3) << "+i*" << imag(x3) << endl;*/
	
	/*if( abs(real(x1)) > pow(10, REAL_THRES)*abs(imag(x1)) )
		roots.push_back(real(x1));
	
	if( abs(real(x2)) > pow(10, REAL_THRES)*abs(imag(x2)) )
		roots.push_back(real(x2));
	
	if( abs(real(x3)) > pow(10, REAL_THRES)*abs(imag(x3)) )
		roots.push_back(real(x3));

	if(roots.size())
		return true;
	else
		return false;*/
}

bool solveQuartic(const double a, const double b, const double c, const double d, const double e, vector<double>& roots, bool verbose){
	// See https://en.wikipedia.org/wiki/Quartic_function#Solving_by_factoring_into_quadratics 
	
	if(!a)
		return solveCubic(b, c, d, e, roots, verbose);

    if(!b && !c && !d){
        roots.push_back(0.);
        roots.push_back(0.);
        roots.push_back(0.);
        roots.push_back(0.);
    }else{
        const double an = b/a;
        const double bn = c/a - (3./8.) * SQ(b/a);
        const double cn = CB(0.5*b/a) - 0.5*b*c/SQ(a) + d/a;
        const double dn = -3.*QU(0.25*b/a) + e/a - 0.25*b*d/SQ(a) + c*SQ(b/4.)/CB(a);

        vector<double> res;
        solveCubic(1., 2.*bn, SQ(bn) - 4.*dn, -SQ(cn), res, verbose);
        unsigned short pChoice = -1;

        for(unsigned short i = 0; i<res.size(); ++i){
            if(res[i] > 0){
                pChoice = i;
                break;
            }
        }

        if(pChoice < 0){
            if(verbose)
                cout << "No real solution to " << a << " x^4 + " << b << " x^3 + " << c << " x^2 + " << d << " x + " << e << " (no positive root for the resolvent cubic)." << endl << endl;
            return false;
        }

        const double p = sqrt(res[pChoice]);
        solveQuadratic(p, SQ(p), 0.5*(p*(bn + res[pChoice]) - cn), roots, verbose);
        solveQuadratic(p, -SQ(p), 0.5*(p*(bn + res[pChoice]) + cn), roots, verbose);

        for(unsigned short i = 0; i<roots.size(); ++i)
            roots[i] -= an/4.;
    }

    const unsigned short nRoots = roots.size();

    if(verbose){
        if(nRoots){
            cout << "Solutions of " << a << " x^4 + " << b << " x^3 + " << c << " x^2 + " << d << " x + " << e << ":" << endl;
            for(unsigned short i=0; i<nRoots; i++)
	            cout << "x" << i << " = " << roots[i] << ", test = " << a*QU(roots[i]) + b*CB(roots[i]) + c*SQ(roots[i]) + d*roots[i] + e << endl;
            cout << endl;
        }else
            cout << "No real solution to " << a << " x^4 + " << b << " x^3 + " << c << " x^2 + " << d << " x + " << e << endl << endl;
    }

    return nRoots;
	
	/*const cd delta = 256.*pow(a,3.)*pow(e,3.) - 192.*pow(a,2.)*b*d*pow(e,2.) - 128.*pow(a*c*e,2.)
		+ 144.*c*e*pow(a*d,2.) - 27.*pow(a,2.)*pow(d,4.) + 144.*a*c*pow(b*e,2.) - 6.*a*e*pow(b*d,2.)
		- 80.*a*b*d*e*pow(c,2.) + 18.*a*b*c*pow(d,3.) + 16.*a*e*pow(c,4.) - 4.*a*pow(c,3.)*pow(d,2.)
		- 27.*pow(b,4.)*pow(e,2.) + 18.*c*d*e*pow(b,3.) - 4.*pow(d*b,3.) - 4.*e*pow(b,2.)*pow(c,3.)
		+ pow(b*c*d,2.);
	
	const cd delta0 = SQ(c) - 3.*b*d + 12.*a*e;
	const cd delta1 = 2.*CB(c) - 9.*b*c*d + 27.*SQ(b)*e + 27.*a*SQ(d) - 72.*a*c*e;

    cd temp_deltas;
    if(delta != 0. && delta0 == 0.)
        temp_deltas = delta1;
    else
        temp_deltas = sqrt( SQ(delta1) - 4.*CB(delta0) );

    if(delta == 0. && delta0 == 0.){
        cout << "Warning: special case in solveQuartic, not covered!" << endl;
        return false;
    }

	cd Q = cbrt( ( delta1 + temp_deltas )/2. );
	
    const cd p = (8.*a*c - 3.*SQ(b))/(8.*SQ(a));
	const cd q = (CB(b) - 4.*a*b*c + 8.*SQ(a)*d)/(8*CB(a));

	cd S = 0.5 * sqrt( -2./3. * p + 1./(3.*a) * ( Q + delta0/Q ) );

    while(S == 0.){
        // Q is not zero, so this should not be an inifite loop
        Q *= exp( cd(0.,2.*M_PI/3.) );
        S = 0.5 * sqrt( -2./3. * p + 1./(3.*a) * ( Q + delta0/Q ) );
    }

	// four solutions
	
	cd x1 = -b/(4.*a) - S + 0.5 * sqrt( -4.*SQ(S) - 2.*p + q/S );
	cd x2 = -b/(4.*a) - S - 0.5 * sqrt( -4.*SQ(S) - 2.*p + q/S );

	cd x3 = -b/(4.*a) + S + 0.5 * sqrt( -4.*SQ(S) - 2.*p - q/S );
	cd x4 = -b/(4.*a) + S - 0.5 * sqrt( -4.*SQ(S) - 2.*p - q/S );

    // reorganize them s.t., if there are double roots, they are grouped (0,1) and (2,3)
    //                                    triple roots, they are grouped (0,1,2) and (3)
    if(real(x1) == real(x3)){
        cd temp = x2;
        x2 = x3;
        x3 = temp;
    }
    if(real(x1) == real(x4)){
        cd temp = x2;
        x2 = x4;
        x4 = temp;
    }
    if(real(x2) == real(x4)){
        cd temp = x3;
        x3 = x4;
        x4 = temp;
    }

    //cout << "Solutions of " << a << " x^4 + " << b << " x^3 + " << c << " x^2 + " << d << " x + " << e << endl;
	//if( abs(real(x1)) > pow(10,5)*abs(imag(x1)) && abs(real(x1)) < pow(10,8)*abs(imag(x1)) && abs(real(x1)) != 0 )
	    //cout << "x1 = " << real(x1) << "+i*" << imag(x1) << endl;
	//if( abs(real(x2)) > pow(10,5)*abs(imag(x2)) && abs(real(x2)) < pow(10,8)*abs(imag(x2)) && abs(real(x2)) != 0 )
	    //cout << "x2 = " << real(x2) << "+i*" << imag(x2) << endl;
	//if( abs(real(x3)) > pow(10,5)*abs(imag(x3)) && abs(real(x3)) < pow(10,8)*abs(imag(x3)) && abs(real(x3)) != 0 )
	    //cout << "x3 = " << real(x3) << "+i*" << imag(x3) << endl;
	//if( abs(real(x4)) > pow(10,5)*abs(imag(x4)) && abs(real(x4)) < pow(10,8)*abs(imag(x4)) && abs(real(x4)) != 0 )
	    //cout << "x4 = " << real(x4) << "+i*" << imag(x4) << endl;

	// checking which ones are real, with some tolerance
	
	if( abs(real(x1)) > pow(10, REAL_THRES)*abs(imag(x1)) )
		roots.push_back(real(x1));
	
	if( abs(real(x2)) > pow(10, REAL_THRES)*abs(imag(x2)) )
		roots.push_back(real(x2));
	
	if( abs(real(x3)) > pow(10, REAL_THRES)*abs(imag(x3)) )
		roots.push_back(real(x3));
	
	if( abs(real(x4)) > pow(10, REAL_THRES)*abs(imag(x4)) )
		roots.push_back(real(x4));

	if(roots.size())
		return true;
	else
		return false;*/
}

bool solve2Quads(const double a11, const double a22, const double a12, const double a10, const double a01, const double a00,
				const double b11, const double b22, const double b12, const double b10, const double b01, const double b00,
				vector<double>& E1, vector<double>& E2,
                bool verbose){

	const double alpha = b11*a22-a11*b22;
	const double beta = b11*a12-a11*b12;
	const double gamma = b11*a10-a11*b10;
	const double delta = b11*a01-a11*b01;
	const double omega = b11*a00-a11*b00;

	const double a = a11*SQ(alpha) + a22*SQ(beta) - a12*alpha*beta;
	const double b = 2.*a11*alpha*delta - a12*alpha*gamma - a12*delta*beta - a10*alpha*beta + 2.*a22*beta*gamma + a01*SQ(beta);
	const double c = a11*SQ(delta) + 2.*a11*alpha*omega - a12*delta*gamma - a12*omega*beta - a10*alpha*gamma - a10*delta*beta
		+ a22*SQ(gamma) + 2.*a01*beta*gamma + a00*SQ(beta);
	const double d = 2.*a11*delta*omega - a12*omega*gamma - a10*delta*gamma - a10*omega*beta + a01*SQ(gamma) + 2.*a00*beta*gamma;
	const double e = a11*SQ(omega) - a10*omega*gamma + a00*SQ(gamma);

	solveQuartic(a, b, c, d, e, E2, verbose);

	for(unsigned int i = 0; i < E2.size(); ++i){
		const double e2 = E2[i];
        if(beta*e2 + gamma != 0){        // Everything OK
		    const double e1 = -(alpha * SQ(e2) + delta*e2 + omega)/(beta*e2 + gamma);
		    E1.push_back(e1);
        }else if(alpha*SQ(e2) + delta*e2 + omega == 0){
            // Up to two solutions for e1 (aligned along y-axis)
            vector<double> e1;
            if( !solveQuadratic(a11, a12*e2 + a10, a22*SQ(e2) + a01*e2 + a00, e1, verbose) ){
                cout << "Error in solve2Quads: there should be at least one solution for e1!" << endl;
                exit(1);
            }
            if(e1[0] == e1[1]){ // no ambiguity
                E1.push_back(e1[0]);
                if(i < E2.size() - 1){
                    if(e2 == E2[i+1]){
                        E1.push_back(e1[0]);
                        E1.push_back(e1[0]);
                        i++;
                        continue;
                    }
                }
            }else{ // in that case, e2 SHOULD be twice degenerate!
                if(i < E2.size() - 1){
                    if(e2 != E2[i+1]){
                        cout << "Error in solve2Quads: if there are two solutions for e1, e2 should be degenerate!" << endl;
                        exit(1);
                    }
                    if(i < E2.size() - 2){
                        if(e2 == E2[i+2]){
                            cout << "Error in solve2Quads: if there are two solutions for e1, e2 cannot be thrice degenerate!" << endl;
                            exit(1);
                        }
                    }
                    E1.push_back(e1[0]);
                    E1.push_back(e1[1]);
                    i++;
                    continue;
                }else{
                    cout << "Error in solve2Quads: if there are two solutions for e1, e2 should be degenerate!" << endl;
                    exit(1);
                }
            }
        }else{        // There is no solution given this e2
            E2.erase(E2.begin() + i);
            cout << "Error in solve2Quads: no solution!" << endl;
        }
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

