//
// Created by lurker on 5/27/17.
//

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
