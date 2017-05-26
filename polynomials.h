//
// Created by lurker on 5/14/17.
//

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "utils.h"
#include "linalg.h"

using namespace bbfmm;

/*
 * can be used for quad precision. todo
 * each value is computed with O(order).
 */

scalar_t jacobi(int order, scalar_t alpha, scalar_t beta, scalar_t x);

scalar_t koornwinder(int k, int n, scalar_t x, scalar_t y);

void makeKoornwinderMatrix(Matrix &K, int N, vector<scalar_t> &x, vector<scalar_t> &y, vector<scalar_t> &w);


#endif //ISORTES_JACOBI_POLYNOMIAL_H
