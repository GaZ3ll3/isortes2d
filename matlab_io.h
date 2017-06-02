/*
 * output csv subroutine.
 *
 * copyright@ Yimin Zhong. yzhong@math.utexas.edu.
 */

#ifndef MATLAB_IO_H
#define MATLAB_IO_H

#include "utils.h"
#include "bbfmm.h"
#include "linalg.h"

void write_to_csv(vector<scalar_t> data, std::string filename) {
    std::ofstream outfile;
    outfile.open(filename);
    for (auto it : data) {
        outfile << it << "\n";
    }
    outfile.close();
}


void write_to_csv(bbfmm::Vector data, std::string filename) {
    std::ofstream outfile;
    outfile.open(filename);
    for (int id = 0; id < data.row(); ++id) {
        outfile << data(id) << "\n";
    }
    outfile.close();
}

void write_to_csv(vector<point> data, std::string filename, std::string sep = " ") {
    std::ofstream outfile;
    outfile.open(filename);
    for (auto it : data) {
        outfile << it.x << sep << it.y << "\n";
    }
    outfile.close();
}


#endif //MATLAB_IO_H
