#ifndef _INC_JACOBIAND
#define _INC_JACOBIAND

#include <vector>
#include "TLorentzVector.h"

#define INV_JAC_MIN 1e3 // Just as in MW

int ComputeTransformD(const double s13, const double s134, const double s25, const double s256,
                      const TLorentzVector p3, const TLorentzVector p4, const TLorentzVector p5, const TLorentzVector p6, const TLorentzVector Met,
                      std::vector<TLorentzVector> &p1, std::vector<TLorentzVector> &p2);

  double computeJacobianD(const std::vector<TLorentzVector> p, const double sqrt_s);

#endif
