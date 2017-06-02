/*
 * BiCGStab subroutine.
 *
 * copyright@ Yimin Zhong. yzhong@math.utexas.edu. All Rights Reserved.
 *
 */



#include "bicgstab.h"

using namespace bbfmm;

void bicgstab(std::function<Vector(Vector &)> A, Vector &x, Vector &b,
              const int _maxIter, const double _tol) {

    /*
       * for performance use
       */
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    int iter = 0;
    int flag = 0;

    scalar_t bnrm2 = nrm2(b);

    if (bnrm2 == 0.) bnrm2 = 1.0;

    Vector r = b;
    Vector c = A(x);
    daxpy(-1.0, c, r);

    scalar_t error = nrm2(r) / bnrm2;
    if (error < _tol) return;

    scalar_t omega = 1.0;
    Vector r_tld = r;

    scalar_t rho_1 = 0.;
    scalar_t alpha = 0.;
    scalar_t resid = 0.;
    scalar_t rho = 0.;
    Vector p(r.row());
    setValue(p, 0.);
    Vector v(r.row());
    setValue(v, 0.);
    Vector s(r.row());
    setValue(s, 0.);

    std::cout << "=============== BICGSTAB =============" << std::endl;
    std::cout << "    iter    |  rel error   |   time   " << std::endl;

    begin = std::chrono::steady_clock::now();

    for (iter = 0; iter < _maxIter; ++iter) {

        end = std::chrono::steady_clock::now();
        std::cout << std::setw(6) << iter << std::setw(20) << std::scientific << sqrt(nrm2(r) / bnrm2) << std::setw(12)
                  << std::fixed
                  << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0
                  << std::fixed << std::endl;
        begin = std::chrono::steady_clock::now();

        rho = ddot(r_tld, r);
        if (rho == 0.) break;

        if (iter > 0) {
            scalar_t beta = (rho / rho_1) * (alpha / omega);
            daxpy(-omega, v, p);
            dscal(beta, p);
            daxpy(1.0, r, p);
        } else {
            p = r;
        }

        v = A(p);

        alpha = rho / (ddot(r_tld, v));

        s = r;
        daxpy(-alpha, v, s);

        if (nrm2(s) < _tol) {
            daxpy(alpha, p, x);
            resid = nrm2(s) / bnrm2;
            break;
        }

        Vector t = A(s);

        omega = (ddot(t, s)) / (ddot(t, t));

        daxpy(alpha, p, x);
        daxpy(omega, s, x);

        r = s;
        daxpy(-omega, t, r);
        error = nrm2(r) / bnrm2;

        if (error <= _tol) {
            break;
        }
        if (omega == 0.) {
            std::cout << "!!!" << std::endl;
            break;
        }

        rho_1 = rho;

    }

    if (error <= _tol || nrm2(s) <= _tol) {
        if (nrm2(s) <= _tol) {
            error = nrm2(s) / bnrm2;
        }
        flag = 0;
        std::cout << "=============== CONVERGED =============" << std::endl;

    } else if (omega == 0.) {
        flag = -2;
        std::cout << "============= BREAK DOWN 1  ===========" << std::endl;
    } else if (rho == 0.) {
        flag = -1;
        std::cout << "============= BREAK DOWN 2  ===========" << std::endl;
    } else {
        flag = 1;
        std::cout << "============= NOT CONVERGED ===========" << std::endl;
    }
}