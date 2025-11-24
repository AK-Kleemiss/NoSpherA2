//
// Created by lucas on 11/21/25.
//
#include <algorithm>
#include <cstdio>
#include <occ/qm/wavefunction.h>
#include <Eigen/Eigen>
#include <iostream>
#include <ostream>
#include "../Src/convenience.h"
#include "../Src/wfn_class.h"
using namespace occ::qm;
using std::iostream;
std::string get_shape(const occ::Mat& mat)
{
    return std::format("({}, {})", mat.rows(), mat.cols());
}

template <typename T>
ivec argsort(const std::vector<T>& v) {
    ivec idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::stable_sort(idx.begin(), idx.end(),
       [&v](int i1, int i2) { return v[i1] < v[i2]; });
    return idx;
}
template <typename T>
std::vector<T> apply_permutation(const std::vector<T>& v, const ivec& indices) {
    std::vector<T> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = v[indices[i]];
    }
    return result;
}
template <typename T>
std::vector<T> sort_vec_by_other(const std::vector<T>& vec1, const std::vector<T>& other)
{
    assert (vec1.size() == other.size());
    return apply_permutation(vec1, argsort(other));
}


WFN wfn_from_nos(const std::string &filepath)
{
    return WFN(filepath);
}

WFN original_approach(Wavefunction& wfn)
{
    WFN nos_wfn(wfn);
    auto shells = wfn.basis.shells();
    auto shell2atom  = wfn.basis.shell_to_atom();
    vec exponentsnos; // push_back_exponent emulation
    vec con_coefs;
    ivec types; // push_back_type emulation
    ivec center; // push_back_center emulation
    int cumm{1};
    const occ::qm::MolecularOrbitals& mo = wfn.mo;
    int r_u_ro_switch{0};
    switch (mo.kind)
    {
        case occ::qm::General:
            r_u_ro_switch = 2;
            break;
        case occ::qm::Unrestricted:
            r_u_ro_switch = 1;
            break;
        case occ::qm::Restricted:
            r_u_ro_switch = 0;
            break;
    }
    for (int i = 0; i<shells.size(); i++)
    {
        auto shell = shells[i];
        occ::Vec confac;

        auto exponents = shell.exponents;
        std::cout << "Shell: " << i << "\nsize of Exponents:" << get_shape(exponents) << "\n";
        auto contraction = shell.contraction_coefficients;
        std::cout << "Contraction: " << "size: " << get_shape(contraction) << "\n";
        std::cout << "Shell2atom: " << shell2atom[i] << "\n";
        int l = shell.l;
        std::cout << "Shell type, l: " << l << "\n";
        std::cout << "calcs:" << "\n";
        confac = Eigen::pow(pow(2, (4*l+3))*Eigen::pow(exponents.array(), 2*l+3), 0.25);
        auto scaled_contraction = contraction*confac.transpose();
        std::cout << confac << "\n";
        int n_cart = (l+1)*(l+2)/2;
        double exp;
        for (int j = 0; j < exponents.size(); j++)
        {
            exp = exponents(j);
            for (int cart=0; cart < n_cart; cart++)
            {
                con_coefs.push_back(scaled_contraction(j));
                exponentsnos.push_back(exp);
                types.push_back(cumm + cart);
            }
        }
        cumm += n_cart;
        std::cout << contraction.transpose() * exponents.matrix() << "\n";
        std::cout << "----" << std::endl;
    }
    return nos_wfn;
}



int main()
{
    std::string filepath("/home/lucas/CLionProjects/NoSpherA2/occ_integration_tests/alanine.owf.fchk");
    auto wfn_nos = WFN(filepath);
    Wavefunction wfn = Wavefunction::load(filepath);
    auto wfn_from_occ = WFN(wfn);
    std::cout <<wfn_from_occ.sort_wfn(wfn_from_occ.check_order(true), true);
    std::cout << wfn_nos.sort_wfn(wfn_nos.check_order(true), true);
    auto centersvecocc = wfn_from_occ.get_centers();
    auto centersvecnos = wfn_nos.get_centers();

    auto typesvecocc = wfn_from_occ.get_types();
    auto typesvecnos = wfn_nos.get_types();

    auto expsvecocc = wfn_from_occ.get_exponents();
    auto expsvecnos = wfn_nos.get_exponents();
    bool differ = false;
    if (typesvecnos.size() == typesvecocc.size())
    {
        for (int i = 0; i<typesvecnos.size(); i++)
        {
            // std::cout << "idx: " << i << " typesvecnos: " << typesvecnos[i] << " typesvecocc: " << typesvecocc[i] << "\n";
            if (typesvecnos[i] != typesvecocc[i])
            {
                differ=true;
            }
        }
    }
    if (differ)
        std::cout << "[types] from occ constructor and from file differ.\n";
    else
        std::cout << "[types] from occ constructor and from file are equal.\n";


    if (expsvecocc == expsvecnos)
        std::cout << "[exponents] from occ constructor and from file are equal.\n";
    else {
        std::cout << "[exponents] from occ constructor and from file differ.\n";
        for (int i=0; i<expsvecnos.size(); i++)
        {
            if (expsvecnos[i] != expsvecocc[i])
                std::cout << "idx: " << i << " expnos[" << centersvecnos[i] <<"]: " << expsvecnos[i] << " expocc[" << centersvecocc[i] <<"]:  " <<expsvecocc[i] << "\n";
        }
    }
    return 0;
}
void occ_native_approach(Wavefunction& wfn)
{
    auto shell2atom  = wfn.basis.shell_to_atom();
    // original_approach(wfn);
    auto shells = wfn.basis.shells();
    std::vector<occ::Mat> shell_vec;
    int old_atom{-1};
    int new_atom{0};
    double norm;
    for (int i = 0; i<shells.size(); i++)
    {
        auto shell = shells[i];
        int l = shells[i].l;
        if (old_atom != new_atom)
        {
            std::cout << "---------\n" << "Atom: " << new_atom << "\n";
            old_atom = new_atom;
        }
        new_atom = shell2atom[i];
        std::cout << "shell: " << i << "\n" << shell.coeffs_normalized_for_libecpint() << "\n";
        std::cout << "coefficients gto_norm: " << "\n";

        auto exponents = shells[i].exponents;
        for (int j=0; j< shell.contraction_coefficients.size(); j++)
        {
            norm = occ::qm::gto_norm(l, exponents(j));
            std::cout << "gto_norm : " << norm*shell.contraction_coefficients(j) << "\n";
        }
    }
}
