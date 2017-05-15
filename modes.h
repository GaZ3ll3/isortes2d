//
// Created by lurker on 5/12/17.
//

#ifndef MODES_H
#define MODES_H


#include <vector>
#include <quadmath.h>


using std::vector;

class QRule_tri {
public:
    vector<long double> points_x;
    vector<long double> points_y;
    vector<long double> weights;
    size_t degree;

    void resize(size_t n) {
        points_x.resize(n);
        points_y.resize(n);
        weights.resize(n);
    }
};

void get_vr_data(size_t deg, QRule_tri &table);

void affine(QRule_tri &table);


#endif //ISORTES_MODES_H
