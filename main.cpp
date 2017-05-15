
/*
 *  Steady-State Radiative Transport Equation Solver for Isotropic Scattering in 2D.
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

#define DISP

#include <iostream>
#include "bbfmm.h"
#include "geometry.h"
#include "modes.h"
#include "polynomials.h"

using namespace bbfmm;

extern void generate_domain(Array<double_t> &, Array<int> &, std::string, geometry *);

extern void domain_information(geometry *);

extern void modes_information(QRule_tri &);

extern void
makeKoornwinderMatrix(Matrix &K, int N, vector<long double> &x, vector<long double> &y, vector<long double> &w);

const int S = 6;

int main() {
    // 1. generate domain mesh.
    scalar_t area[8] = {0.5, 0.125, 0.03125, 0.0078125, 0.001953125, 0.00048828125, 0.0001220703125,
                        0.0001220703125 / 4.0};

    geometry *domain = new geometry();
    Array<scalar_t> boundary_points = {0., 0., 0., 1., 1., 1., 1., 0.};
    Array<int> edge_connection = {0, 1, 1, 2, 2, 3, 3, 0};
    generate_domain(boundary_points, edge_connection, "prq34.0a" + std::to_string(area[S]) + "nzjQ", domain);
    domain->getAllArea();
    domain_information(domain);

    // 1-1. get near field.
    vector<vector<int>> nearField;
    domain->getNearField(nearField);

    // 2. get high order quadrature rule on standard reference triangle T1.
    QRule_tri quadratue_rule;
    get_vr_data(10, quadratue_rule);
    affine(quadratue_rule);
    modes_information(quadratue_rule);
    index_t coarseQuadratureSize = quadratue_rule.weights.size();


    // 3. generate source, target nodes.
    vector<point> source, target;
    vector<scalar_t> weights;

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

    // 4. FMM
    // calculate particle interactions, near field is included.
    index_t np = 9;
    /*
     * G.eval (source ,target).
     */
    kernel G;
    G.eval = [&](point &a, point &b) {
        scalar_t dist = sqrt(SQR(b.x - a.x) + SQR(b.y - a.y));
        if (dist == 0.) return 0.;
        return 1.0 / dist;
    };

    Vector p((int) source.size());
    setValue(p, 1.0);
    for (int i = 0; i < p.row(); ++i) {

        p(i) *= source[i].x * source[i].y * weights[i];
    }

    RUN("quadtree", G.initialize(np, source, target, p, (index_t) source.size(), (index_t) target.size(), np * np, 10));
    Vector ret;
    G.run(ret);

    std::cout << std::setprecision(15) << std::scientific << ret(0) << std::endl;

    std::cout << std::setprecision(15) << std::scientific << source[0].x << " " << source[0].y << std::endl;

    // 5. Remove near interaction.

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    Vector coarseNearInteraction((int) target.size());
    setValue(coarseNearInteraction, 0.);

    for (int targetTriangleId = 0; targetTriangleId < domain->numberoftriangles; ++targetTriangleId) {
        for (int targetQuadratureId = 0; targetQuadratureId < quadratue_rule.weights.size(); ++targetQuadratureId) {

            int targetId = (int) (coarseQuadratureSize * targetTriangleId + targetQuadratureId);

            for (int nearSourceTriangleId = 0;
                 nearSourceTriangleId < nearField[targetTriangleId].size();
                 ++nearSourceTriangleId) {
                for (int nearSourceQuadratureId = 0;
                     nearSourceQuadratureId < quadratue_rule.weights.size(); ++nearSourceQuadratureId) {

                    int nearSourceId = (int) (coarseQuadratureSize * nearField[targetTriangleId][nearSourceTriangleId] +
                                              nearSourceQuadratureId);
                    coarseNearInteraction(targetId) += G.eval(source[nearSourceId], target[targetId]) * p(nearSourceId);

                }
            }
        }
    }

    daxpy(-1.0, coarseNearInteraction, ret);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << std::setw(15) << "REMOVE" << " " << std::setprecision(5) << std::setw(8)
              << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0 << " seconds"
              << std::endl;

    std::cout << std::setprecision(15) << std::scientific << ret(0) << " " << coarseNearInteraction(0) << std::endl;


    // 5. Near field re-calculation.
    // use interpolation for near field.
    // There will be T * 12 * M^2 * (4^l) flops.
    // storage: (4^l * M) x M, precomputed.


    // 5-1 build interpolation matrix, by producting the same matrix.
    // The matrix is of size (4^level * M) x M
    // Time cost is linear, but constant will be huge.
    // The interpolation matrix is uniform.
    Matrix interpolate(coarseQuadratureSize, coarseQuadratureSize);

    int nC = int((-3.0 + sqrt(9 + 8 * coarseQuadratureSize)) / 2.0);
    assert((nC + 1) * (nC + 2) == 2 * coarseQuadratureSize);

    makeKoornwinderMatrix(interpolate, nC, quadratue_rule.points_x, quadratue_rule.points_y, quadratue_rule.weights);

    Vector sqrtW((int) quadratue_rule.weights.size());
    for (int wi = 0; wi < quadratue_rule.weights.size(); ++wi) {
        sqrtW(wi) = (scalar_t) sqrt(quadratue_rule.weights[wi]);
    }


    int nRefineLevel = 0;
    int refineQuadratureSize = coarseQuadratureSize;

    vector<long double> refine_quad_x, refine_quad_y, refine_weight;

    for (int level = 0; level < nRefineLevel; ++level) {
        for (int id = 0; id < refineQuadratureSize; ++id) {
            refine_quad_x.push_back(quadratue_rule.points_x[id] / 2.0);
            refine_quad_y.push_back(quadratue_rule.points_y[id] / 2.0);
            refine_weight.push_back(quadratue_rule.weights[id] / 4.0);

            refine_quad_x.push_back(quadratue_rule.points_x[id] / 2.0);
            refine_quad_y.push_back(quadratue_rule.points_y[id] / 2.0 + 0.50);
            refine_weight.push_back(quadratue_rule.weights[id] / 4.0);

            refine_quad_x.push_back(quadratue_rule.points_x[id] / 2.0 + 0.50);
            refine_quad_y.push_back(quadratue_rule.points_y[id] / 2.0);
            refine_weight.push_back(quadratue_rule.weights[id] / 4.0);

            refine_quad_x.push_back(0.50 - quadratue_rule.points_x[id] / 2.0);
            refine_quad_y.push_back(0.50 - quadratue_rule.points_y[id] / 2.0);
            refine_weight.push_back(quadratue_rule.weights[id] / 4.0);
        }

        quadratue_rule.points_x = refine_quad_x;
        quadratue_rule.points_y = refine_quad_y;
        quadratue_rule.weights = refine_weight;
        refineQuadratureSize *= 4;


        refine_quad_x.clear();
        refine_quad_y.clear();
        refine_weight.clear();
    }

    Vector sqrtWR((int) quadratue_rule.weights.size());
    for (int wi = 0; wi < quadratue_rule.weights.size(); ++wi) {
        sqrtWR(wi) = (scalar_t) sqrt(quadratue_rule.weights[wi]);
    }


    Matrix refinements(coarseQuadratureSize, refineQuadratureSize);
    makeKoornwinderMatrix(refinements, nC, quadratue_rule.points_x, quadratue_rule.points_y, quadratue_rule.weights);

    Matrix mapping(refineQuadratureSize, coarseQuadratureSize);
    t_dgemm(1.0, refinements, interpolate, 0., mapping);


    Vector Rout(refineQuadratureSize);
    setValue(Rout, 0.);
    dgemv(1.0, mapping, sqrtW, 0., Rout);

    daxpy(-1.0, sqrtWR, Rout);

    std::cout << nrm2(Rout) << std::endl;

    // 5-2. Evaluate refined quadrature.
    // todo: improve the application of interpolation matrix. should be linear to depth, not exponential.

    begin = std::chrono::steady_clock::now();

    Vector refineNearInteraction((int) target.size());
    setValue(refineNearInteraction, 0.);

#pragma omp parallel for schedule(static, 80) collapse(2) num_threads(4)
    for (int targetTriangleId = 0; targetTriangleId < domain->numberoftriangles; ++targetTriangleId) {
        for (int targetQuadratureId = 0; targetQuadratureId < coarseQuadratureSize; ++targetQuadratureId) {

            int targetId = (int) (coarseQuadratureSize * targetTriangleId + targetQuadratureId);

            for (int nearSourceTriangleId = 0;
                 nearSourceTriangleId < nearField[targetTriangleId].size() &&
                 nearField[targetTriangleId][nearSourceTriangleId] != targetTriangleId;
                 ++nearSourceTriangleId) {

                auto near_u = domain->getTriangle(3 * nearField[targetTriangleId][nearSourceTriangleId]);
                auto near_v = domain->getTriangle(3 * nearField[targetTriangleId][nearSourceTriangleId] + 1);
                auto near_w = domain->getTriangle(3 * nearField[targetTriangleId][nearSourceTriangleId] + 2);

//                scalar_t cur_area =  domain->arealist[nearField[targetTriangleId][nearSourceTriangleId]];

                Vector newValues(refineQuadratureSize);
                Vector oldValues(coarseQuadratureSize);

                for (int sourceQuadratureId = 0; sourceQuadratureId < coarseQuadratureSize; ++sourceQuadratureId) {
                    oldValues(sourceQuadratureId) =
                            p(coarseQuadratureSize * nearField[targetTriangleId][nearSourceTriangleId] +
                              sourceQuadratureId) / sqrtW(sourceQuadratureId);
                }

                setValue(newValues, 0.);

                dgemv(1.0, mapping, oldValues, 0., newValues);

                for (int nearSourceQuadratureId = 0;
                     nearSourceQuadratureId < refineQuadratureSize; ++nearSourceQuadratureId) {

                    scalar_t lambda = (scalar_t) quadratue_rule.points_x[nearSourceQuadratureId];
                    scalar_t mu = (scalar_t) quadratue_rule.points_y[nearSourceQuadratureId];
                    scalar_t w = (scalar_t) quadratue_rule.weights[nearSourceQuadratureId];


                    point cur_nearSourcePoint = {
                            domain->getX(near_u) * lambda +
                            domain->getX(near_v) * mu +
                            domain->getX(near_w) * (1.0 - lambda - mu),
                            domain->getY(near_u) * lambda +
                            domain->getY(near_v) * mu +
                            domain->getY(near_w) * (1.0 - lambda - mu)
                    };

                    refineNearInteraction(targetId) +=
                            G.eval(cur_nearSourcePoint, target[targetId]) * sqrt(w) * newValues(nearSourceQuadratureId);

                }
            }
        }
    }

    end = std::chrono::steady_clock::now();
    std::cout << std::setw(15) << "NEAR" << " " << std::setprecision(5) << std::setw(8)
              << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0 << " seconds"
              << std::endl;

    std::cout << std::setprecision(15) << std::scientific << ret(0) << " " << refineNearInteraction(0) << std::endl;


    // 6. All singular integral.
    // on each triangle, integrate with high order rule with tensored L^2 points
    // There will be ~ T * M^2 * 3L^2 flops.
    // plus storage: (M * 3L^2) * M has to be precomputed.


    // 6-1. get local quadrature rule.




    delete (domain);
    return 0;
}

extern void generate_domain(Array<double_t> &b, Array<int> &e, std::string switches, geometry *out) {
    geometry *setting = new geometry();
    setting->set_points(b);
    setting->set_facets(e);
    geometry *coarse = new geometry();
    setting->build("pczjQ", *coarse);
    delete setting;
    coarse->refine(switches, *out);
    delete coarse;
}

extern void domain_information(geometry *in) {
#ifdef DISP
    std::cout << std::setw(15) << "vertices" << std::setw(8) << in->_meta.numberofpoints << std::endl;
    std::cout << std::setw(15) << "triangles" << std::setw(8) << in->_meta.numberoftriangles << std::endl;
    std::cout << std::setw(15) << "edges" << std::setw(8) << in->_meta.numberofedges << std::endl;
    std::cout << std::setw(15) << "segments" << std::setw(8) << in->_meta.numberofsegments << std::endl;
    std::cout << std::endl;
#endif
}

extern void modes_information(QRule_tri &rule) {
#ifdef DISP
    std::cout << std::setw(15) << "degree" << std::setw(8) << rule.degree << std::endl;
    std::cout << std::setw(15) << "#points" << std::setw(8) << rule.points_x.size() << std::endl;
    std::cout << std::endl;
#endif
}

extern void
makeKoornwinderMatrix(Matrix &K, int N, vector<long double> &x, vector<long double> &y, vector<long double> &w) {
    // assume matrix size is correct.
    int row = 0;
    for (int n = 0; n <= N; ++n) {
        for (int k = 0; k <= n; ++k) {
            for (int I = 0; I < x.size(); ++I) {
                K(row, I) = (scalar_t) (koornwinder(k, n, x[I], y[I]) * sqrt(w[I]));
            }
            row += 1;
        }
    }

    for (int r = 0; r < K.row(); ++r) {
        long double nrm = 0.;
        for (int c = 0; c < K.col(); ++c) {
            nrm += SQR(K(r, c));
        }
        nrm = sqrt(nrm);
        for (int c = 0; c < K.col(); ++c) {
            K(r, c) /= nrm;
        }
    }
}