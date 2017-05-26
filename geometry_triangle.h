/*
 * Quadrature and discretization for triangle in 2D.
 *
 * copyright@ Yimin Zhong <yzhong@math.utexas.edu>. All Rights Reserved.
 *
 *
 *  1. volumetric discretization with triangle cells.
 *  2. high order quadrature rules. (Rokhlin or Gimbatus).
 *
 */

#ifndef GEOMETRY_TRIANGLE_H
#define GEOMETRY_TRIANGLE_H


#include <iostream>
#include <cassert>
#include <cstring>
#include <unordered_map>
#include "utils.h"

extern "C" {
#include "triangle.h"
}

class geometry_triangle {
public:
    geometry_triangle(const geometry_triangle &) = delete;

    geometry_triangle();

    ~geometry_triangle();

    void set_points(Array<double> &_points);

    void set_facets(Array<int> &_facets);

    void build(std::string, geometry_triangle &out);

    void refine(std::string, geometry_triangle &out);

    scalar_t getX(int index);

    scalar_t getY(int index);

    void getAllArea();

    scalar_t getDistance(int fromIndex, int toIndex);

    void getNearField(vector<vector<int>> &near);

    int getTriangle(int index);

    int numberofpoints;
    int numberoftriangles;
    vector<scalar_t> arealist;

    struct triangulateio _meta;
private:
    scalar_t getArea(int index);

    void setProperty() {
        numberofpoints = _meta.numberofpoints;
        numberoftriangles = _meta.numberoftriangles;
    }

};


#endif //GEOMETRY_TRIANGLE_H
