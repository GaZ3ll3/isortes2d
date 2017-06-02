/*
 * configuration input subroutine.
 *
 * copyright@ Yimin Zhong. yzhong@math.utexas.edu. All Rights Reserved.
 *
 */

#ifndef CONFIG_H
#define CONFIG_H

#include "utils.h"
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>


class Config {
public:
    Config();

    ~Config();

    std::map<std::string, std::string> options;

    void parse(std::istream &cfgFile);

    void print();
};

#endif //CONFIG_H
