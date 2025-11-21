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

std::string get_shape(const occ::Mat& mat)
{
    return std::format("({}, {})", mat.rows(), mat.cols());
}

int main()
{
    using std::iostream;
    Wavefunction wfn = Wavefunction::load("/home/lucas/CLionProjects/NoSpherA2/temp_tests/alanine.owf.fchk");
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
        auto scaled_contraction = contraction*confac;
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

    // int total_count = 0;
    // for (auto shell : shells)
    // {
    //     for (int i=0; i< shell.num_primitives(); i++)
    //     {
    //         std::cout << shell.coeff_normalized(0,i) << " ";
    //         total_count++;
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "\nTotal Count: " << total_count << std::endl;
    // std::cout << shells[0].contraction_coefficients << std::endl;
    // wfn.save("/home/lucas/CLionProjects/NoSpherA2/temp_tests/alanine_saved.owf.fchk");
    // occ::Mat S_AB = wfn.compute_overlap_matrix();
      // ABo.symmetric_orthonormalize_molecular_orbitals(S_AB);
    // std::cout << wfn.mo.D << std::endl;
    return 0;
}
