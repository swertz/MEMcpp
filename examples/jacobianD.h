#ifndef _INC_JACOBIAND
#define _INC_JACOBIAND

#include <vector>
#include "TLorentzVector.h"

#define JAC_MIN 10e-20

double computeJacobianD(std::vector<TLorentzVector> p, double sqrt_s);

#endif
