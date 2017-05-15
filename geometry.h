/*
 * Quadrature and discretization in 2D.
 *
 * copyright@ Yimin Zhong <yzhong@math.utexas.edu>. All Rights Reserved.
 *
 *
 *  1. volumetric discretization with triangle cells.
 *  2. high order quadrature rules. (Rokhlin or Gimbatus).
 *
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H


#include <iostream>
#include <cassert>
#include <cstring>
#include <unordered_map>
#include "utils.h"

extern "C" {
#include "triangle.h"
}

template<typename T>
class Array {
protected:
    T *data;
    size_t size;
    bool owner;
public:
    Array(size_t _n) {
        data = (T *) malloc(_n * sizeof(T));
        memset(data, 0, _n * sizeof(T));
        size = _n;
        owner = true;
    }

    Array(size_t _n, bool _o) {
        data = (T *) malloc(_n * sizeof(T));
        size = _n;
        owner = _o;
    }

    Array(std::initializer_list<T> l) {
        size_t index = 0;
        for (auto it = l.begin(); it != l.end(); it++) {
            index++;
        }

        data = (T *) malloc(index * sizeof(T));
        T *ptr = data;
        for (auto it = l.begin(); it != l.end(); it++) {
            *(ptr++) = *it;
        }
        size = index;
        owner = true;
    }

    Array(T *_data, size_t _n, bool _o) {
        if (_o) {
            data = (T *) malloc(_n * sizeof(T));
            memcpy(data, _data, _n * sizeof(T));
            size = _n;
            owner = _o;
        } else {
            data = _data;
            size = _n;
            owner = _o;
        }
    }

    ~Array() {
        if (owner) {
            free(data);
            data = nullptr;
        }
    }

    size_t get_size() {
        return size;
    }

    bool get_owner() {
        return owner;
    }

    T *get_data() {
        return data;
    }

    T &operator[](size_t i) {
        assert(i < size);
        return data[i];
    }

    const T &operator[](size_t i) const {
        assert(i < size);
        return data[i];
    }
};

class geometry {
public:
    geometry(const geometry &) = delete;

    geometry();

    ~geometry();

    void set_points(Array<double> &_points);

    void set_facets(Array<int> &_facets);

    void build(std::string, geometry &out);

    void refine(std::string, geometry &out);

    scalar_t getX(int index);

    scalar_t getY(int index);

    scalar_t getArea(int index);

    void getAllArea();

    scalar_t getDistance(int fromIndex, int toIndex);

    void getNearField(vector<vector<int>> &near);

    int getTriangle(int index);

    int numberofpoints;
    int numberoftriangles;
    vector<scalar_t> arealist;

    struct triangulateio _meta;
private:
    void setProperty() {
        numberofpoints = _meta.numberofpoints;
        numberoftriangles = _meta.numberoftriangles;
    }

};


#endif //GEOMETRY_H
