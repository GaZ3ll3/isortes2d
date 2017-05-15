#include "polynomials.h"

scalar_t jacobi(int order, scalar_t alpha, scalar_t beta, scalar_t x) {
    scalar_t prev, cur, next;
    scalar_t a1, a2, a3, a4;
    prev = cur = 0.;
    next = 1.0;
    if (order > 0) {
        cur = next;
        next = (alpha - beta + (alpha + beta + 2.0) * x) / 2.0;
    }
    if (order > 1) {
        for (int n = 1; n < order; n++) {
            prev = cur;
            cur = next;
            a1 = 2 * (n + 1) * (n + alpha + beta + 1) * (2 * n + alpha + beta);
            a2 = (2 * n + alpha + beta + 1) * (alpha * alpha - beta * beta);
            a3 = (2 * n + alpha + beta) * (2 * n + alpha + beta + 1) * (2 * n + alpha + beta + 2);
            a4 = 2 * (n + alpha) * (n + beta) * (2 * n + alpha + beta + 2);
            next = ((a2 + a3 * x) * cur - a4 * prev) / (a1);
        }
    }
    return next;
}

scalar_t koornwinder(int n, int k, scalar_t x, scalar_t y) {
    return jacobi(n - k, 2 * k + 1, 0, 2 * x - 1) * jacobi(k, 0, 0, 2 * y / (1 - x) - 1.0) * pow(1 - x, k);
}
