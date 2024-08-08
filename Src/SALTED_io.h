#pragma once

#include "convenience.h"
#include "npy.h"

#if has_RAS
#include "H5Cpp.h"
#include <H5public.h>
#include <H5Dpublic.h>
#include <H5Opublic.h>

vec2 readHDF5(H5::H5File file, std::string dataset_name);
#endif

template <class T>
std::vector<T> readVectorFromFile(const std::string &filename);

template <typename T>
std::unordered_map<int, std::vector<T>> read_fps(std::string filename, int lmax_max);

template <typename Scalar>
void read_npy(std::string filename, std::vector<Scalar> &data);

struct Config
{
    std::string predict_filename;
    bool average;
    bool field;
    bool sparsify;
    int ncut;
    std::vector<std::string> species;
    double rcut1;
    double rcut2;
    int nang1;
    int nang2;
    int nrad1;
    int nrad2;
    double sig1;
    double sig2;
    std::vector<std::string> neighspe1;
    std::vector<std::string> neighspe2;
    double zeta;
    int Menv;
    int Ntrain;
    float trainfrac;
    std::string dfbasis;

    // The are initialized to -1 to check if they are set
    // They are set in the populateFromFile function as static_cast<int>(neighspeX.size())
    int nspe1 = -1;
    int nspe2 = -1;

    void populateFromFile(const std::string &filename);
    std::vector<std::string> parseVector(const std::string &str);
};