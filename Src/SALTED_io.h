#pragma once

#include "convenience.h"
#include "npy.h"

std::filesystem::path find_first_salted_file(const std::filesystem::path &directory_path);

template <class T>
std::vector<T> readVectorFromFile(const std::filesystem::path &filename);

template <typename T>
std::unordered_map<int, std::vector<T>> read_fps(std::filesystem::path &filename, int lmax_max);

template <typename Scalar>
void read_npy(std::filesystem::path &filename, std::vector<Scalar> &data);


struct Config
{
public:
    bool from_binary;
    std::filesystem::path salted_filename;
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
private:
    std::vector<std::string> parseVector(const std::string &str);
};


class SALTED_BINARY_FILE {
private:
    const std::string MAGIC_NUMBER = "SALTD";
    static const int HEADER_SIZE = 5;
    enum DataType { INT32 = 0, FLOAT64 = 1, STRING = 2 };

    std::filesystem::path filepath;
    std::ifstream file;

    int32_t version;
    int32_t numBlocks;
    std::map<std::string, size_t> table_of_contents;
	int conf_location = 0; //Initialized after gen_config was called
	int header_end = -1;

    void open_file();
    bool read_header();
    

    template <typename T>
    void read_dataset(std::vector<T>& data, ivec& dims);

    std::string read_string_remove_NULL(const int lengh);

    template <typename T>
    T read_generic_blocks(const std::string& key, std::function<void(T&, int)> process_block);
    std::unordered_map<std::string, vec2> read_lambda_based_data(const std::string& key);

public:
    SALTED_BINARY_FILE(const std::filesystem::path& fpath) : filepath(fpath) {
        open_file();
        err_checkf(read_header(), "Error reading header!", std::cout);
    };

    ~SALTED_BINARY_FILE() {
        if (file.is_open()) {
            file.close();
        }
    };

    void populate_config(Config &config_in);


    std::unordered_map<int, std::vector<int64_t>> read_fps();
    std::unordered_map<std::string, vec> read_averages();
    std::unordered_map<int, vec> read_wigners();
    vec read_weights();
    std::unordered_map<std::string, vec2> read_projectors();
    std::unordered_map<std::string, vec2> read_features();
};