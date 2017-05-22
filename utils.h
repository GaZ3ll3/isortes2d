/*
 *  utility for general use.
 *
 *  copyright@ Yimin Zhong <yzhong@math.utexas.edu> All Rights Reserved.
 *
 */

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <utility>
#include <functional>
#include <chrono>
#include <iomanip>
#include "cblas.h"

#if !defined __extern_always_inline && defined __clang__
# if defined __GNUC_STDC_INLINE__ || defined __GNUC_GNU_INLINE__
#  define __extern_inline extern __inline __attribute__ ((__gnu_inline__))
#  define __extern_always_inline \
  extern __always_inline __attribute__ ((__gnu_inline__))
# else
#  define __extern_inline extern __inline
#  define __extern_always_inline extern __always_inline
# endif
#endif

#ifdef RUN_OMP

#include "omp.h"

#endif

#define EPS 1e-12
#define DIM 2 /* 2D FMM, do not change */

#define SQR(X) ((X)*(X))

#ifdef DISP
#define RUN(s, func){ \
std::chrono::steady_clock::time_point begin =std::chrono::steady_clock::now(); \
func;\
std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();\
std::cout << std::setw(15)<< s << " "  << std::setprecision(5) << std::setw(8) << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000.0 << " seconds"<<std::endl;\
}
#else
#define RUN(s, func){\
func;\
}
#endif


typedef int index_t;
typedef double scalar_t;
typedef bool bool_t;

using std::unordered_set;
using std::vector;
using std::queue;
using std::unordered_map;
using std::unordered_set;
using std::max;
using std::min;
using std::sqrt;
using std::acos;
using std::atan;

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


#endif //UTILS_H
