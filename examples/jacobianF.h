#ifndef _INC_JACOBIAND
#define _INC_JACOBIAND

#include <vector>
#include "TLorentzVector.h"

#define JAC_MIN 1.e-20

double computeJacobianF(const std::vector<TLorentzVector> p, const double sqrt_s);

#endif
