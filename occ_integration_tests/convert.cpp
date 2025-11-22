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
using occ::qm::Wavefunction;
using std::iostream;
std::string get_shape(const occ::Mat& mat)
{
    return std::format("({}, {})", mat.rows(), mat.cols());
}

void original_approach(Wavefunction& wfn)
{
    auto shells = wfn.basis.shells();
    auto shell2atom  = wfn.basis.shell_to_atom();
    vec exponentsnos; // push_back_exponent emulation
    vec con_coefs;
    ivec types; // push_back_type emulation
    ivec center; // push_back_center emulation
    int cumm{1};
    for (int i = 0; i<shells.size(); i++)
    {
        occ::Vec confac;
        auto exponents = shells[i].exponents;
        std::cout << "Shell: " << i << "\nsize of Exponents:" << get_shape(exponents) << "\n";
        auto contraction = shells[i].contraction_coefficients;
        std::cout << "Contraction: " << "size: " << get_shape(contraction) << "\n";
        std::cout << "Shell2atom: " << shell2atom[i] << "\n";
        int l = shells[i].l;
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

int main()
{
    Wavefunction wfn = Wavefunction::load("/home/lucas/CLionProjects/NoSpherA2/occ_integration_tests/alanine.owf.fchk");

    return 0;
}
