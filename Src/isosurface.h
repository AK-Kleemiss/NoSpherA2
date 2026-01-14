#pragma once
#include "convenience.h"

// Suppose we store each face color as an (R,G,B) triple
typedef std::array<int, 3> RGB;
class cube;

// A simple triangle struct composed of three vertices
class Triangle {
private:
    d3 v1, v2, v3;
    RGB colour;
    int colour_index;
public:
    Triangle(d3 v1, d3 v2, d3 v3) : v1(v1), v2(v2), v3(v3), colour({ 0, 0, 0 }), colour_index(0) {};
    Triangle(d3 v1, d3 v2, d3 v3, RGB colour) : v1(v1), v2(v2), v3(v3), colour(colour), colour_index(0) {};
    Triangle(d3 v1, d3 v2, d3 v3, RGB colour, int colour_index) : v1(v1), v2(v2), v3(v3), colour(colour), colour_index(colour_index) {};
    Triangle() = default;
    Triangle(const Triangle&) = default;
    Triangle& operator=(const Triangle&) = default;
    ~Triangle() = default;

    RGB get_colour() const { return colour; };
    void set_colour(const RGB& given) { colour = given; };
    void set_colour_index(const int& index) { colour_index = index; };
    void set_colour_index(const size_t& index) { colour_index = (int)index; };
    int get_colour_index() const { return colour_index; };
    d3 get_v(const int& nr) const {
        if (nr == 1)
            return v1;
        else if (nr == 2)
            return v2;
        else if (nr == 3)
            return v3;
        else
            return { 0, 0, 0 };
    };
    double calc_area() const {
        d3 a = { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
        d3 b = { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };
        return 0.5 * array_length(cross(a, b));
    };
    double calc_inner_volume() const {
        d3 a = { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
        d3 b = { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };
        return a_dot(cross(b, a), v1) / 6.;
    };
    d3 calc_center() const {
        return { (v1[0] + v2[0] + v3[0]) / 3., (v1[1] + v2[1] + v3[1]) / 3., (v1[2] + v2[2] + v3[2]) / 3. };
    };
};
std::vector<Triangle> marchingCubes(const cube& volumeData, const double isoVal);//, const int subdivisionLevel = 2);
bool writeObj(const std::filesystem::path& filename, const std::vector<Triangle>& triangles);
bool writeColourObj(const std::filesystem::path& filename, std::vector<Triangle>& triangles);
bool writeMTL(const std::string& mtlFilename, std::vector<Triangle>& triangles);
void get_colour(Triangle& t, const cube& volumeData, std::array<std::array<int, 3>, 3> Colourcode, double low_lim, double high_lim);
void get_colour(Triangle& t, double(*func)(const d3&, const WFN&), const WFN& wavy, std::array<std::array<int, 3>, 3> Colourcode, double low_lim, double high_lim);
double calc_d_i(const d3& p_t, const WFN& wavy);
