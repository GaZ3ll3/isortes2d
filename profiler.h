//
// Created by lurker on 5/24/17.
//

#ifndef ISORTES_PROFILER_H
#define ISORTES_PROFILER_H

#include <stack>
#include "utils.h"

class profiler {
public:
    profiler() {
        _timing_set.clear();
        _total_time = 0.;
        clocking = false;
    }

    ~profiler() {
        // display all time table.
        for (auto kv : _timing_set) {
            std::cout << std::fixed << std::setprecision(3) << std::setw(30) << kv.first << std::setw(15)
                      << kv.second / _total_time * 100 << "%" << std::setw(15) << kv.second << " seconds" << std::endl;
        }

        std::cout << std::setw(30) << "counted time" << std::setw(31) << _total_time << " seconds" << std::endl;
    }

    void tic(std::string s = "") {
        if (clocking) return;

        clocking = true;
        cur_task = s;

        if (_timing_set.find(cur_task) == _timing_set.end()) {
            _timing_set[cur_task] = 0.;
        }

        _begin = std::chrono::steady_clock::now();
    }

    void toc(bool count = true) {
        if (!clocking) return;
        _end = std::chrono::steady_clock::now();
        clocking = false;
        scalar_t local_time = std::chrono::duration_cast<std::chrono::microseconds>(_end - _begin).count() / 1000000.0;
        _timing_set[cur_task] += local_time;
        if (count) {
            _total_time += local_time;
        }
    }

private:
    std::chrono::steady_clock::time_point _begin;
    std::chrono::steady_clock::time_point _end;

    bool clocking;
    scalar_t _total_time;
    std::string cur_task;

    unordered_map<std::string, scalar_t> _timing_set;
    unordered_map<std::string, int> _timing_tag;
};


#endif //ISORTES_PROFILER_H