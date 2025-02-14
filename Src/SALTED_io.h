#pragma once

#include "convenience.h"
#include "npy.h"
#include "H5Cpp.h"
#include <H5public.h>
#include <H5Dpublic.h>
#include <H5Opublic.h>

// vec2 readHDF5(H5::H5File file, std::string dataset_name);

template <typename T>
void readHDF5Data(H5::DataSet &dataset, std::vector<T> &data);

template <typename T>
std::vector<T> readHDF5(H5::H5File file, std::string dataset_name, std::vector<hsize_t> &dims_out);
//template <typename T>
//Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent>> readHDF5(H5::H5File file, std::string dataset_name, std::vector<hsize_t>& dims_out);

std::filesystem::path find_first_h5_file(const std::filesystem::path &directory_path);

template <class T>
std::vector<T> readVectorFromFile(const std::filesystem::path &filename);

template <typename T>
std::unordered_map<int, std::vector<T>> read_fps(std::filesystem::path &filename, int lmax_max);

template <typename Scalar>
void read_npy(std::filesystem::path &filename, std::vector<Scalar> &data);

struct Config
{
public:
    bool from_h5;
    std::filesystem::path h5_filename = "";
    std::filesystem::path predict_filename;
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
    double trainfrac;
    std::string dfbasis;

    // The are initialized to -1 to check if they are set
    // They are set in the populateFromFile function as static_cast<int>(neighspeX.size())
    int nspe1 = -1;
    int nspe2 = -1;

    void populateFromFile(const std::filesystem::path &filename);
    void populateFromFile(const H5::H5File file);

private:
    std::vector<std::string> parseVector(const std::string &str);
    void populate_config(const std::filesystem::path &dataset_name, const int &data);
    void populate_config(const std::filesystem::path &dataset_name, const float &data);
    void populate_config(const std::filesystem::path &dataset_name, const double &data);
    void populate_config(const std::filesystem::path &dataset_name, const std::string &data);
    void populate_config(const std::filesystem::path &dataset_name, const std::vector<std::string> &data);

    void handle_int_dataset(const std::filesystem::path &dataset_name, H5::DataSet &dataSet);
    void handle_float_dataset(const std::filesystem::path &dataset_name, H5::DataSet &dataSet);
    void handle_string_dataset(const std::filesystem::path &dataset_name, H5::DataSet &dataSet);

    static herr_t op_func(hid_t loc_id, const char *name, const H5L_info_t *info, void *operator_data);

    std::unordered_map<H5T_class_t, std::function<void(const std::string &, H5::DataSet &)>> handlers = {
        {H5T_INTEGER, std::bind(&Config::handle_int_dataset, this, std::placeholders::_1, std::placeholders::_2)},
        {H5T_FLOAT, std::bind(&Config::handle_float_dataset, this, std::placeholders::_1, std::placeholders::_2)},
        {H5T_STRING, std::bind(&Config::handle_string_dataset, this, std::placeholders::_1, std::placeholders::_2)},
        {H5T_ENUM, std::bind(&Config::handle_int_dataset, this, std::placeholders::_1, std::placeholders::_2)}};
};
