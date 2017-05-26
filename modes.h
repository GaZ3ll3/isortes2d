//
// Created by lurker on 5/12/17.
//

#ifndef MODES_H
#define MODES_H


#include <vector>
#include <quadmath.h>
#include "utils.h"


using std::vector;

class QRule_tri {
public:
    vector<scalar_t> points_x;
    vector<scalar_t> points_y;
    vector<scalar_t> weights;
    size_t degree;

    void resize(size_t n) {
        points_x.resize(n);
        points_y.resize(n);
        weights.resize(n);
    }
};


class QRule_lin {
public:
    vector<scalar_t> points_x;
    vector<scalar_t> weights;

    void resize(size_t n) {
        points_x.resize(n);
        weights.resize(n);
    }
};

void get_legendre_data(size_t deg, QRule_lin &table);

void affine(QRule_lin &table);

void get_vr_data(size_t deg, QRule_tri &table);

void affine(QRule_tri &table);


#endif //ISORTES_MODES_H
