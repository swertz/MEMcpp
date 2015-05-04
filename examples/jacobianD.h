#ifndef _INC_JACOBIAND
#define _INC_JACOBIAND

#include <vector>
#include "TLorentzVector.h"

#define INV_JAC_MIN 1e3 // Just as in MW

double computeJacobianD(const std::vector<TLorentzVector> p, const double sqrt_s);

#endif
