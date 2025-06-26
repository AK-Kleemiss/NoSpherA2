#pragma once

#include <filesystem>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <sstream>

class primitive
{

public:
    primitive(int c, int t, std::vector<double> e, std::vector<double> coef) : center(c), type(t), exp(e), coefficient(coef) {};
    primitive(int c, int t, double e, double coef) : center(c), type(t), exp({e}), coefficient({coef}) {};

    int center, type;
    std::vector<double> exp, coefficient;
};


template <class T>
T fromString(const std::string& s)
{
    std::istringstream stream(s);
    T t;
    stream >> t;
    return t;
}

template <class T>
std::vector<T> split_string(const std::string& input, const std::string delimiter)
{
    std::string input_copy = input + delimiter; // Need to add one delimiter in the end to return all elements
    std::vector<T> result;
    size_t pos = 0;
    while ((pos = input_copy.find(delimiter)) != std::string::npos)
    {
        result.emplace_back(fromString<T>(input_copy.substr(0, pos)));
        input_copy.erase(0, pos + delimiter.length());
    }
    return result;
};

const std::array<std::vector<primitive>, 118> read_basis_set(std::filesystem::path basis_path);
std::string write_basis_set(std::ofstream& file, const std::string basis_name, std::array<std::vector<primitive>, 118> basis_set);