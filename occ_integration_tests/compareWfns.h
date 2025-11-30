//
// Created by lucas on 11/30/25.
//

#ifndef NOSPHERA2_COMPAREWFNS_H
#define NOSPHERA2_COMPAREWFNS_H

struct VecSize
{
    int x{1};
    int y{1};
    bool is_empty() const { return (x==1) && (y==1); }
};

std::string get_shape(const occ::Mat& mat)
{
    return std::format("({}, {})", mat.rows(), mat.cols());
}

template <typename T>
bool compare_vectors(std::vector<T> vecNOS, std::vector<T> vecOCC, std::string label,
    VecSize size = VecSize(), double tol = 1e-12, bool compare_diff_size = false, bool print_on_success = true);

#endif //NOSPHERA2_COMPAREWFNS_H
