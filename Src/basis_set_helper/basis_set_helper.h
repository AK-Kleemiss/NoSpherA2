#pragma once

#include <filesystem>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <unordered_map>

class primitive
{
private:
    int center, type;
    double exp, coefficient;
    double norm_const = -10;

public:
    primitive() : center(0), type(0), exp(0.0), coefficient(0.0) {};
    primitive(int c, int t, double e, double coef) : center(c), type(t), exp(e), coefficient(coef) {};

    int get_center() const
    {
        return center;
    };
    int get_type() const
    {
        return type;
    };
    double get_exp() const
    {
        return exp;
    };
    double get_coef() const
    {
        return coefficient;
    };
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
        result.push_back(fromString<T>(input_copy.substr(0, pos)));
        input_copy.erase(0, pos + delimiter.length());
    }
    return result;
};

const std::array<std::vector<primitive>, 118> read_basis_set(std::filesystem::path basis_path);
void write_basis_sets(const std::filesystem::path basis_dir, const std::unordered_map<std::string, std::array<std::vector<primitive>, 118>> basis_sets);