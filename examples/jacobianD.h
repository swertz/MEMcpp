#ifndef _INC_JACOBIAND
#define _INC_JACOBIAND

#include <vector>
#include "Math/Vector4D.h"

#define INV_JAC_MIN 1e3 // Just as in MW

int ComputeTransformD(const double s13, const double s134, const double s25, const double s256,
                      const ROOT::Math::PxPyPzEVector p3, const ROOT::Math::PxPyPzEVector p4, const ROOT::Math::PxPyPzEVector p5, const ROOT::Math::PxPyPzEVector p6, const ROOT::Math::PxPyPzEVector Met,
                      std::vector<ROOT::Math::PxPyPzEVector> &p1, std::vector<ROOT::Math::PxPyPzEVector> &p2);

double computeJacobianD(const std::vector<ROOT::Math::PxPyPzEVector> p, const double sqrt_s);

#endif
