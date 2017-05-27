/*
 * bicgstab.h
 *
 *  Created on: Oct 9, 2016
 *      Author: Yimin Zhong
 */

#ifndef BICGSTAB_H
#define BICGSTAB_H

#include "blas_wrapper.h"

void bicgstab(std::function<Vector(Vector &)> f, Vector &x, Vector &rhs,
              const int _maxIter, const double _tol);

#endif //BICGSTAB_H
