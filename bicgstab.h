/*
 * BiCGStab subroutine.
 *
 * copyright@ Yimin Zhong. yzhong@math.utexas.edu. All Rights Reserved.
 *
 */


#ifndef BICGSTAB_H
#define BICGSTAB_H

#include "blas_wrapper.h"

void bicgstab(std::function<Vector(Vector &)> f, Vector &x, Vector &rhs,
              const int _maxIter, const double _tol);

#endif //BICGSTAB_H
