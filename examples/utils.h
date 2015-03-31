#ifndef _INC_UTILS
#define _INC_UTILS

#include <vector>

double solveQuadratic(double a, double b, double c, std::vector<double>& roots);
double solveCubic(double a, double b, double c, double d, std::vector<double>& roots);
bool solveQuartic(double a, double b, double c, double d, double e, std::vector<double>& roots);
bool solve2Quads(double a11, double a22, double a12, double a10, double a01, double a00,
				double b11, double b22, double b12, double b10, double b01, double b00,
				std::vector<double>& E1, std::vector<double>& E2);
double BreitWigner(double s, double m, double g);

#endif
