#ifndef _INC_UTILS
#define _INC_UTILS

#include <vector>

bool solveQuadratic(const double a, const double b, const double c, std::vector<double>& roots);
bool solveCubic(const double a, const double b, const double c, const double d, std::vector<double>& roots);
bool solveQuartic(const double a, const double b, const double c, const double d, const double e, std::vector<double>& roots);
bool solve2Quads(const double a11, const double a22, const double a12, const double a10, const double a01, const double a00,
				const double b11, const double b22, const double b12, const double b10, const double b01, const double b00,
				std::vector<double>& E1, std::vector<double>& E2);
double BreitWigner(const double s, const double m, const double g);

#endif
