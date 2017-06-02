
/*
 *  Steady-State Radiative Transport Equation Solver for Isotropic Scattering in 2D unit square.
 *
 *  copyright@ Yimin Zhong. yzhong@math.utexas.edu.
 *
 *  About:
 *
 *  1. This program solves RTE by volumetric integral formulation, where integral is evaluated by FMM.
 *
 *  2. The integral equation is second kind, solved by Krylov subspace method e.g. GMRES, MINRES, BICGSTAB.
 *
 *  3. Integral part is based on the paper "A NYSTRÖM METHOD FOR WEAKLY SINGULAR INTEGRAL OPERATORS ON SURFACES" of
 *  J.Bremer and Z.Gimbutas.
 *
 *  4. FMM method uses BlackBox Fast Multipole Method by Darve and Fong.
 *
 *
 *  Limitation:
 *
 *  Support int32 particle number only.
 *
 */
#include <iostream>
#include "bbfmm.h"
#include "modes.h"
#include "polynomials.h"
#include "unit_square_mesh.h"
#include "profiler.h"
#include "gmres.h"
#include "Config.h"
#include "bicgstab.h"
#include "matlab_io.h"

using namespace bbfmm;

extern void
duffy_transform(int deg, vector<scalar_t> points, vector<scalar_t> &X, vector<scalar_t> &Y, vector<scalar_t> &W);

int main(int argc, char *argv[]) {
#ifdef RUN_OMP
    omp_set_num_threads(omp_get_max_threads());
#endif
    // configuration.

    if (argc <= 1) {
        std::cout << "USE " << argv[0] << " PATH_OF_CONFIG_FILE " << std::endl;
        exit(0);
    }

    Config cfg;
    std::ifstream cfgFile;
    cfgFile.open(argv[1], std::ifstream::in);
    cfg.parse(cfgFile);
    cfgFile.close();

    cfg.print();

    int M = atoi(cfg.options["M"].c_str());
    int qrule_tri_deg = atoi(cfg.options["q_tri"].c_str());
    int qrule_lin_deg = atoi(cfg.options["q_lin"].c_str());
    index_t np = atoi(cfg.options["np"].c_str());

    scalar_t mu_t = 2.2;
    scalar_t mu_s = 2.0;

    int nRefineLevel = 2;


    auto source_function = [&](scalar_t x, scalar_t y) {
        scalar_t d = sqrt(SQR(x - 0.6) + SQR(y - 0.4));
        if (fabs(d - 0.2) <= 0.05) {
            return (cos((d - 0.2) * M_PI / 0.05) + 1.0) / 2.0;
        } else {
            return 0.;
        }
    };

    vector<int> sample(8);

    profiler prof;


    /*
     * build mesh and related quantities.
     */

    prof.tic("domain generation");
    geometry_triangle *domain = build_unit_square_mesh(M);
    prof.toc();

    /*
     * get near field for near interaction.
     */

    prof.tic("near-field assignment");
    vector<vector<int>> nearField;
    domain->getNearField(nearField);
    prof.toc();

    /*
     * get standard quadrature rule for equilateral triangle
     * and affine onto reference triangle [0, 1] -- [1, 0] -- [0, 0].
     */
    prof.tic("quadrature rule");
    QRule_tri quadratue_rule;
    get_vr_data((size_t) qrule_tri_deg, quadratue_rule);
    affine(quadratue_rule);
    assert(quadratue_rule.weights.size() > 0);
    prof.toc();


    int coarseQuadratureSize = (int) quadratue_rule.weights.size();
    int refineQuadratureSize = coarseQuadratureSize;

    int nC = int((-3.0 + sqrt(9 + 8 * coarseQuadratureSize)) / 2.0);
    assert((nC + 1) * (nC + 2) == 2 * coarseQuadratureSize);


    Matrix interpolate(coarseQuadratureSize, coarseQuadratureSize);
    makeKoornwinderMatrix(interpolate, nC, quadratue_rule.points_x, quadratue_rule.points_y, quadratue_rule.weights);

    Vector sqrtW((int) quadratue_rule.weights.size());
    for (int wi = 0; wi < quadratue_rule.weights.size(); ++wi) {
        sqrtW(wi) = (scalar_t) sqrt(quadratue_rule.weights[wi]);
    }

    vector<scalar_t> refine_quad_x = quadratue_rule.points_x;
    vector<scalar_t> refine_quad_y = quadratue_rule.points_y;
    vector<scalar_t> refine_weight = quadratue_rule.weights;

    vector<scalar_t> refine_quad_x_tmp, refine_quad_y_tmp, refine_weight_tmp;

    prof.tic("refinement");
    for (int level = 0; level < nRefineLevel; ++level) {
        for (int id = 0; id < refineQuadratureSize; ++id) {
            refine_quad_x_tmp.push_back(refine_quad_x[id] / 2.0);
            refine_quad_y_tmp.push_back(refine_quad_y[id] / 2.0);
            refine_weight_tmp.push_back(refine_weight[id] / 4.0);


            refine_quad_x_tmp.push_back(refine_quad_x[id] / 2.0);
            refine_quad_y_tmp.push_back(refine_quad_y[id] / 2.0 + 0.50);
            refine_weight_tmp.push_back(refine_weight[id] / 4.0);

            refine_quad_x_tmp.push_back(refine_quad_x[id] / 2.0 + 0.50);
            refine_quad_y_tmp.push_back(refine_quad_y[id] / 2.0);
            refine_weight_tmp.push_back(refine_weight[id] / 4.0);

            refine_quad_x_tmp.push_back(0.50 - refine_quad_x[id] / 2.0);
            refine_quad_y_tmp.push_back(0.50 - refine_quad_y[id] / 2.0);
            refine_weight_tmp.push_back(refine_weight[id] / 4.0);

        }

        refine_quad_x = refine_quad_x_tmp;
        refine_quad_y = refine_quad_y_tmp;
        refine_weight = refine_weight_tmp;
        refineQuadratureSize *= 4;


        refine_quad_x_tmp.clear();
        refine_quad_y_tmp.clear();
        refine_weight_tmp.clear();
    }


    Matrix refinements(coarseQuadratureSize, refineQuadratureSize);
    makeKoornwinderMatrix(refinements, nC, refine_quad_x, refine_quad_y, refine_weight);

    Matrix mapping(refineQuadratureSize, coarseQuadratureSize);
    t_dgemm(1.0, refinements, interpolate, 0., mapping);
    prof.toc();


    vector<scalar_t> X1, Y1, W1, X2, Y2, W2, X3, Y3, W3;



    Matrix basis(coarseQuadratureSize, coarseQuadratureSize);
    setValue(basis, 0.);

    vector<scalar_t> norm_of_koornwinder(quadratue_rule.weights.size());


    prof.tic("singular precomputing");
#pragma omp parallel for num_threads(omp_get_max_threads()) private(X1, X2, X3, Y1, Y2, Y3, W1, W2, W3)
    for (int row = 0; row <= nC * (nC + 3) / 2; ++row) {
        int n = int(sqrt(2 * row + 0.25) - 0.5);
        int k = row - (n + 1) * n / 2;
        for (int I = 0; I < coarseQuadratureSize; ++I) {
            norm_of_koornwinder[row] += SQR(
                    koornwinder(n, k, quadratue_rule.points_x[I], quadratue_rule.points_y[I]) *
                    sqrt(quadratue_rule.weights[I]));


            duffy_transform(qrule_lin_deg, {quadratue_rule.points_x[I], quadratue_rule.points_y[I], 0, 1, 0, 0}, X1,
                            Y1,
                            W1);
            duffy_transform(qrule_lin_deg, {quadratue_rule.points_x[I], quadratue_rule.points_y[I], 1, 0, 0, 1}, X2,
                            Y2,
                            W2);
            duffy_transform(qrule_lin_deg, {quadratue_rule.points_x[I], quadratue_rule.points_y[I], 0, 0, 1, 0}, X3,
                            Y3,
                            W3);


            for (int q = 0; q < W1.size(); ++q) {
                scalar_t dist = sqrt(SQR(X1[q] - quadratue_rule.points_x[I]) +
                                     SQR(Y1[q] - quadratue_rule.points_y[I]));
                basis(row, I) += koornwinder(n, k,
                                             X1[q], Y1[q]) * W1[q] / (dist);
            }


            for (int q = 0; q < W2.size(); ++q) {
                scalar_t dist = sqrt(SQR(X2[q] - quadratue_rule.points_x[I]) +
                                     SQR(Y2[q] - quadratue_rule.points_y[I]));
                basis(row, I) += koornwinder(n, k,
                                             X2[q], Y2[q]) * W2[q] / (dist);

            }

            for (int q = 0; q < W3.size(); ++q) {
                scalar_t dist = sqrt(SQR(X3[q] - quadratue_rule.points_x[I]) +
                                     SQR(Y3[q] - quadratue_rule.points_y[I]));
                basis(row, I) += koornwinder(n, k,
                                             X3[q], Y3[q]) * W3[q] / (dist);
            }
        }
        norm_of_koornwinder[row] = sqrt(norm_of_koornwinder[row]);
    }

#pragma omp parallel for collapse(2)
    for (int _row = 0; _row < basis.row(); ++_row) {
        for (int _col = 0; _col < basis.col(); ++_col) {
            basis(_row, _col) /= norm_of_koornwinder[_row];
        }
    }


    Matrix local(coarseQuadratureSize, coarseQuadratureSize);
    t_dgemm(1.0, basis, interpolate, 0., local);
    prof.toc(false);


    /*
     * get source and target points for second kind integral equation.
     * source equals to target.
     */
    vector<point> source, target;
    vector<scalar_t> weights;

    prof.tic("source/target allocation");

    for (index_t i = 0; i < domain->numberoftriangles; ++i) {

        auto _u = domain->getTriangle(3 * i);
        auto _v = domain->getTriangle(3 * i + 1);
        auto _w = domain->getTriangle(3 * i + 2);

        auto _area = domain->arealist[i];

        for (index_t j = 0; j < coarseQuadratureSize; ++j) {
            auto lambda = (scalar_t) quadratue_rule.points_x[j];
            auto mu = (scalar_t) quadratue_rule.points_y[j];
            point p = {
                    domain->getX(_u) * lambda +
                    domain->getX(_v) * mu +
                    domain->getX(_w) * (1.0 - lambda - mu),
                    domain->getY(_u) * lambda +
                    domain->getY(_v) * mu +
                    domain->getY(_w) * (1.0 - lambda - mu)
            };

            source.push_back(p);
            target.push_back(p);

            weights.push_back(2.0 * _area * quadratue_rule.weights[j]);
        }
    }
    prof.toc();


    for (int sample_id = 0; sample_id < source.size(); ++sample_id) {
        if (fabs(source[sample_id].x - 1. / 6.) + fabs(source[sample_id].y - 1. / 6.) < 1E-9) {
            sample[0] = sample_id;
        }
        if (fabs(source[sample_id].x - 1. / 6.) + fabs(source[sample_id].y - 2. / 3.) < 1E-9) {
            sample[1] = sample_id;
        }
        if (fabs(source[sample_id].x - 1. / 3.) + fabs(source[sample_id].y - 1. / 3.) < 1E-9) {
            sample[2] = sample_id;
        }
        if (fabs(source[sample_id].x - 2. / 3.) + fabs(source[sample_id].y - 1. / 6.) < 1E-9) {
            sample[3] = sample_id;
        }
        if (fabs(source[sample_id].x - 5. / 6.) + fabs(source[sample_id].y - 5. / 6.) < 1E-9) {
            sample[4] = sample_id;
        }
        if (fabs(source[sample_id].x - 5. / 6.) + fabs(source[sample_id].y - 1. / 3.) < 1E-9) {
            sample[5] = sample_id;
        }
        if (fabs(source[sample_id].x - 2. / 3.) + fabs(source[sample_id].y - 2. / 3.) < 1E-9) {
            sample[6] = sample_id;
        }
        if (fabs(source[sample_id].x - 1. / 3.) + fabs(source[sample_id].y - 5. / 6.) < 1E-9) {
            sample[7] = sample_id;
        }

    }

    /*
     * use BBFMM for evaluation of particle interaction.
     * the kernel in 2D has smooth part and determinstic singular part.
     * for smooth part, there is nothing extra to do.
     * for singular part, we adopt Bremer & Gimbutas method.
     *
     * kernel function signature is G.eval (source ,target).
     *
     *
     * mu_t can be variable coefficient in (x,y).
     */

    auto smooth_eval = [&](point &src, point &trg) {
        scalar_t dist = sqrt(SQR(src.x - trg.x) + SQR(src.y - trg.y));
        if (dist == 0.) return -mu_t;
        else {
            return (exp(-mu_t * dist) - 1.0) / (dist);
        }
    };


    auto singular_eval = [&](point &src, point &trg) {
        scalar_t dist = sqrt(SQR(src.x - trg.x) + SQR(src.y - trg.y));
        if (dist == 0.) return 0.;
        return (1.0) / (dist);
    };


    auto smooth_mapping = [&](Vector &charge) {

        kernel smooth_kernel;
        smooth_kernel.eval = smooth_eval;

        assert(charge.row() == source.size());
        Vector weightedCharge(charge.row());
        setValue(weightedCharge, 0.);

        for (auto chargeId = 0; chargeId < weightedCharge.row(); ++chargeId) {
            weightedCharge(chargeId) = charge(chargeId) * weights[chargeId];
        }

        prof.tic("quadtree generation");
        smooth_kernel.initialize(np, source, target, weightedCharge, (index_t) source.size(), (index_t) target.size(),
                                 np * np, 10);
        prof.toc();

        Vector ret;

        prof.tic("fast multipole method");
        smooth_kernel.run(ret);
        prof.toc();

        return ret;
    };


    auto singular_mapping = [&](Vector &charge) {

        kernel singular_kernel;
        singular_kernel.eval = singular_eval;

        assert(charge.row() == source.size());
        Vector weightedCharge(charge.row());
        setValue(weightedCharge, 0.);

        for (auto chargeId = 0; chargeId < weightedCharge.row(); ++chargeId) {
            weightedCharge(chargeId) = charge(chargeId) * weights[chargeId];
        }

        prof.tic("quadtree generation");
        singular_kernel.initialize(np, source, target, weightedCharge, (index_t) source.size(), (index_t) target.size(),
                                   np * np, 10);
        prof.toc();

        Vector ret;

        prof.tic("fast multipole method");
        singular_kernel.run(ret);
        prof.toc();

        Vector coarseNearInteraction((int) target.size());
        setValue(coarseNearInteraction, 0.);

        prof.tic("near-field removal");
        for (int targetTriangleId = 0; targetTriangleId < domain->numberoftriangles; ++targetTriangleId) {
            for (int targetQuadratureId = 0; targetQuadratureId < quadratue_rule.weights.size(); ++targetQuadratureId) {

                int targetId = (int) (coarseQuadratureSize * targetTriangleId + targetQuadratureId);

                for (int nearSourceTriangleId = 0;
                     nearSourceTriangleId < nearField[targetTriangleId].size();
                     ++nearSourceTriangleId) {
                    for (int nearSourceQuadratureId = 0;
                         nearSourceQuadratureId < quadratue_rule.weights.size(); ++nearSourceQuadratureId) {

                        int nearSourceId =
                                (int) (coarseQuadratureSize * nearField[targetTriangleId][nearSourceTriangleId] +
                                       nearSourceQuadratureId);
                        coarseNearInteraction(targetId) +=
                                singular_kernel.eval(source[nearSourceId], target[targetId]) *
                                weightedCharge(nearSourceId);
                    }
                }
            }
        }
        daxpy(-1.0, coarseNearInteraction, ret);
        prof.toc();

        Vector refineNearInteraction((int) source.size());
        setValue(refineNearInteraction, 0.);

        prof.tic("near interaction addon");
#pragma omp parallel for schedule(static, CHUNKSIZE) collapse(2) num_threads(omp_get_max_threads())
        for (int targetTriangleId = 0; targetTriangleId < domain->numberoftriangles; ++targetTriangleId) {
            for (int targetQuadratureId = 0; targetQuadratureId < coarseQuadratureSize; ++targetQuadratureId) {

                int targetId = coarseQuadratureSize * targetTriangleId + targetQuadratureId;

                for (int nearSourceTriangleId = 0;
                     nearSourceTriangleId < nearField[targetTriangleId].size();
                     ++nearSourceTriangleId) {

                    if (nearField[targetTriangleId][nearSourceTriangleId] == targetTriangleId) continue;

                    auto near_u = domain->getTriangle(3 * nearField[targetTriangleId][nearSourceTriangleId]);
                    auto near_v = domain->getTriangle(3 * nearField[targetTriangleId][nearSourceTriangleId] + 1);
                    auto near_w = domain->getTriangle(3 * nearField[targetTriangleId][nearSourceTriangleId] + 2);

                    Vector newValues(refineQuadratureSize);
                    Vector oldValues(coarseQuadratureSize);

                    for (int sourceQuadratureId = 0; sourceQuadratureId < coarseQuadratureSize; ++sourceQuadratureId) {
                        oldValues(sourceQuadratureId) =
                                weightedCharge(
                                        coarseQuadratureSize * nearField[targetTriangleId][nearSourceTriangleId] +
                                        sourceQuadratureId) / sqrtW(sourceQuadratureId);
                    }

                    setValue(newValues, 0.);

                    dgemv(1.0, mapping, oldValues, 0., newValues);

                    for (int nearSourceQuadratureId = 0;
                         nearSourceQuadratureId < refineQuadratureSize; ++nearSourceQuadratureId) {

                        scalar_t lambda = (scalar_t) refine_quad_x[nearSourceQuadratureId];
                        scalar_t mu = (scalar_t) refine_quad_y[nearSourceQuadratureId];
                        scalar_t w = (scalar_t) refine_weight[nearSourceQuadratureId];


                        point cur_nearSourcePoint = {
                                domain->getX(near_u) * lambda +
                                domain->getX(near_v) * mu +
                                domain->getX(near_w) * (1.0 - lambda - mu),
                                domain->getY(near_u) * lambda +
                                domain->getY(near_v) * mu +
                                domain->getY(near_w) * (1.0 - lambda - mu)
                        };

                        refineNearInteraction(targetId) +=
                                singular_kernel.eval(cur_nearSourcePoint, target[targetId]) * sqrt(w) *
                                newValues(nearSourceQuadratureId);

                    }
                }
            }
        }
        prof.toc();

        daxpy(1.0, refineNearInteraction, ret);

        prof.tic("singular calculation");
#pragma omp parallel for  schedule(static, CHUNKSIZE) collapse(2) num_threads(omp_get_max_threads())
        for (int tid = 0; tid < domain->numberoftriangles; ++tid) {
            for (int qid = 0; qid < coarseQuadratureSize; ++qid) {
                for (int rid = 0; rid < coarseQuadratureSize; ++rid) {
                    ret(tid * coarseQuadratureSize + qid) +=
                            charge(tid * coarseQuadratureSize + rid) * sqrtW(rid) * local(qid, rid) / M;
                }
            }
        }
        prof.toc();

        return ret;

    };


    auto apply_integral = [&](Vector &in) {
        Vector singular_part = singular_mapping(in);
        Vector out = smooth_mapping(in);

        daxpy(1.0, singular_part, out);
        dscal(M_1_PI / 2.0, out);
        return out;
    };


    auto forward_mapping = [&](Vector &in) {
        Vector out = in;
        Vector tmp = apply_integral(in);
        dscal(mu_s, tmp);
        daxpy(-1.0, tmp, out);
        return out;
    };

    Vector rhs((int) source.size());
    setValue(rhs, 0.);
    for (int i = 0; i < rhs.row(); ++i) {
        rhs(i) = source_function(source[i].x, source[i].y);
    }

    Vector RHS = apply_integral(rhs);
    Vector x(RHS.row());
    setValue(x, 0.);


    if (strcmp(cfg.options["Krylov"].c_str(), "GMRES") == 0) {
        GMRES(forward_mapping, x, RHS, 20, 400, 1e-14);
    } else {
        bicgstab(forward_mapping, x, RHS, 400, 1e-14);
    }

    if (atoi(cfg.options["IO"].c_str())) {
        write_to_csv(source, "points.csv", " ");
        write_to_csv(x, "result.csv");
    }


    for (auto sample_id : sample) {
        std::cout << std::scientific << std::setprecision(16) << x(sample_id) << std::endl;
    }

    std::cout << "number of points: " << source.size() << std::endl;
    std::cout << "triangles: " << domain->numberoftriangles << std::endl;



    delete (domain);
    return 0;
}

extern void
duffy_transform(int deg, vector<scalar_t> points, vector<scalar_t> &X, vector<scalar_t> &Y, vector<scalar_t> &W) {
    QRule_lin ql;
    get_legendre_data(deg, ql);
    affine(ql);


    scalar_t a = points[0];
    scalar_t b = points[1];

    scalar_t a11 = points[2] - points[0];
    scalar_t a12 = points[4] - points[2];
    scalar_t a21 = points[3] - points[1];
    scalar_t a22 = points[5] - points[3];

    Matrix A(2, 2);
    A(0, 0) = a11;
    A(0, 1) = a12;
    A(1, 0) = a21;
    A(1, 1) = a22;

    scalar_t detA = a11 * a22 - a12 * a21;

    size_t nsp = deg;

    vector<scalar_t> sx(nsp * nsp);
    vector<scalar_t> sy(nsp * nsp);
    vector<scalar_t> sw(nsp * nsp);


    int id = 0;
    for (int ri = 0; ri < nsp; ++ri) {
        for (int ci = 0; ci < nsp; ++ci) {
            sx[id] = ql.points_x[ri];
            sy[id] = ql.points_x[ci];
            sw[id++] = ql.weights[ri] * ql.weights[ci];
        }
    }

    vector<scalar_t> Z1x(nsp * nsp);
    vector<scalar_t> Z1y(nsp * nsp);
    vector<scalar_t> Z1w(nsp * nsp);

    X.resize(nsp * nsp);
    Y.resize(nsp * nsp);
    W.resize(nsp * nsp);


    for (int i = 0; i < nsp * nsp; ++i) {
        scalar_t u = sx[i];
        scalar_t v = sy[i];
        scalar_t w = sw[i];

        Z1x[i] = u;
        Z1y[i] = u * v;
        Z1w[i] = w * u;

        u = Z1x[i];
        v = Z1y[i];
        w = Z1w[i];

        scalar_t x = a11 * u + a12 * v + a;
        scalar_t y = a21 * u + a22 * v + b;
        scalar_t wr = detA * w;

        X[i] = x;
        Y[i] = y;
        W[i] = wr;
    }
}