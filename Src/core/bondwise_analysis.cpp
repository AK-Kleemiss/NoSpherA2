#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "bondwise_analysis.h"
#include "properties.h"
#include "libCintMain.h"
#include "nos_math.h"
#include "integration_params.h"
#include "b2c.h"
#include <occ/qm/hf.h>

namespace {
    struct OhOperation {
        std::array<int, 3> permutation;
        std::array<int, 3> signs;
    };

    struct CartesianTransform {
        ivec destination;
        ivec signs;
    };

    int cartesian_shell_size(const int l) {
        return (l + 1) * (l + 2) / 2;
    }

    int atomic_shell_size(const int l, const bool cartesian) {
        return cartesian ? cartesian_shell_size(l) : 2 * l + 1;
    }

    int tonto_period_number(const int atomic_number) {
        if (atomic_number < 1)
            return 0;

        int period = 1;
        int noble = 0;
        while (true) {
            const int n = (period + 2) / 2;
            noble += 2 * n * n;
            if (atomic_number <= noble)
                return period;
            period++;
        }
    }

    int tonto_column_number(const int atomic_number) {
        if (atomic_number < 1)
            return 0;

        int period = 1;
        int noble = 0;
        while (true) {
            const int n = (period + 2) / 2;
            noble += 2 * n * n;
            if (atomic_number <= noble)
                return atomic_number - (noble - 2 * n * n);
            period++;
        }
    }

    int tonto_ground_state_multiplicity(const int atomic_number) {
        const int period = tonto_period_number(atomic_number);
        const int column = tonto_column_number(atomic_number);
        switch (period) {
        case 0:
            return 1;
        case 1:
        case 2:
        case 3:
            switch (column) {
            case 8: return 1;
            case 1:
            case 3:
            case 7: return 2;
            case 2:
            case 4:
            case 6: return 3;
            case 5: return 4;
            default: return 1;
            }
        case 4:
        case 5:
            switch (column) {
            case 2:
            case 12:
            case 18: return 1;
            case 1:
            case 3:
            case 11:
            case 13:
            case 17: return 2;
            case 4:
            case 10:
            case 14:
            case 16: return 3;
            case 5:
            case 9:
            case 15: return 4;
            case 6:
            case 8: return 5;
            case 7: return 6;
            default: return 1;
            }
        case 6:
        case 7:
            switch (column) {
            case 2:
            case 16:
            case 26:
            case 32: return 1;
            case 1:
            case 3:
            case 15:
            case 17:
            case 25:
            case 27:
            case 31: return 2;
            case 4:
            case 14:
            case 18:
            case 24:
            case 28:
            case 30: return 3;
            case 5:
            case 13:
            case 19:
            case 23:
            case 29: return 4;
            case 6:
            case 12:
            case 20:
            case 22: return 5;
            case 7:
            case 11:
            case 21: return 6;
            case 8:
            case 10: return 7;
            case 9: return 8;
            default: return 1;
            }
        default:
            return 1;
        }
    }

    int tonto_core_electron_count(const int atomic_number) {
        if (atomic_number < 1)
            return 0;

        int noble = 0;
        int period = 1;
        while (true) {
            const int n = (period + 2) / 2;
            const int next_noble = noble + 2 * n * n;
            if (atomic_number <= next_noble)
                return noble;
            noble = next_noble;
            period++;
        }
    }

    occ::gto::AOBasis build_occ_atomic_basis_from_wfn_atom(
        const atom &atm, const e_origin origin, const bool cartesian) {
        std::vector<occ::core::Atom> occ_atoms{ { atm.get_charge(), 0.0, 0.0, 0.0 } };
        std::vector<occ::gto::Shell> shells;
        const auto basis_set = atm.get_basis_set();
        int primitive_idx = 0;

        while (primitive_idx < static_cast<int>(basis_set.size())) {
            const int shell_id = static_cast<int>(basis_set[primitive_idx].get_shell());
            const int l = origin == e_origin::OCC
                ? static_cast<int>(basis_set[primitive_idx].get_type())
                : static_cast<int>(basis_set[primitive_idx].get_type()) - 1;
            err_checkf(l >= 0,
                "Encountered an invalid shell angular momentum while building an OCC atomic basis.",
                std::cout);
            std::vector<double> exponents;
            std::vector<double> coefficients;
            while (primitive_idx < static_cast<int>(basis_set.size()) &&
                static_cast<int>(basis_set[primitive_idx].get_shell()) == shell_id) {
                exponents.push_back(basis_set[primitive_idx].get_exponent());
                coefficients.push_back(basis_set[primitive_idx].get_coefficient());
                primitive_idx++;
            }

            occ::gto::Shell shell(l, exponents, { coefficients }, { 0.0, 0.0, 0.0 });
            shell.kind = cartesian ? occ::gto::Shell::Kind::Cartesian : occ::gto::Shell::Kind::Spherical;
            shell.incorporate_shell_norm();
            shells.push_back(shell);
        }

        occ::gto::AOBasis result(occ_atoms, shells, "wfn-atomic-basis");
        result.set_pure(!cartesian);
        result.set_ecp_electrons({ atm.get_ECP_electrons() });
        return result;
    }

    dMatrix2 eigen_matrix_to_dmatrix2(const occ::Mat &matrix) {
        dMatrix2 result(matrix.rows(), matrix.cols());
        for (int row = 0; row < matrix.rows(); ++row)
            for (int col = 0; col < matrix.cols(); ++col)
                result(row, col) = matrix(row, col);
        return result;
    }

    class ScopedOccLogLevel {
    public:
        explicit ScopedOccLogLevel(const spdlog::level::level_enum level) {
            auto logger = spdlog::default_logger();
            if (logger != nullptr) {
                had_logger = true;
                previous_level = logger->level();
                logger->set_level(level);
            }
        }

        ~ScopedOccLogLevel() {
            if (!had_logger)
                return;
            if (auto logger = spdlog::default_logger(); logger != nullptr)
                logger->set_level(previous_level);
        }

    private:
        bool had_logger = false;
        spdlog::level::level_enum previous_level = spdlog::level::info;
    };

    dMatrix2 compute_tonto_style_atomic_density(
        const atom &atm, const e_origin origin, const bool cartesian) {
        ScopedOccLogLevel quiet_occ_logs(spdlog::level::err);
        const occ::gto::AOBasis basis = build_occ_atomic_basis_from_wfn_atom(atm, origin, cartesian);
        const int effective_atomic_number = atm.get_charge() - atm.get_ECP_electrons();
        const int multiplicity = std::max(1, tonto_ground_state_multiplicity(effective_atomic_number));
        const bool restricted = multiplicity == 1 && (effective_atomic_number % 2 == 0);
        const auto spin_kind = restricted
            ? occ::qm::SpinorbitalKind::Restricted
            : occ::qm::SpinorbitalKind::Unrestricted;

        occ::qm::HartreeFock hf(basis);
        occ::qm::SCF<occ::qm::HartreeFock> scf(hf, spin_kind);
        scf.set_charge_multiplicity(0, multiplicity);
        scf.compute_scf_energy();

        occ::qm::MolecularOrbitals mo =
            occ::io::conversion::orb::to_gaussian_order(basis, scf.wavefunction().mo);
        mo.update_occupied_orbitals();
        mo.update_density_matrix();

        if (spin_kind == occ::qm::SpinorbitalKind::Restricted)
            return eigen_matrix_to_dmatrix2(mo.D);

        const occ::Mat spin_summed =
            occ::qm::block::a(mo.D) + occ::qm::block::b(mo.D);
        return eigen_matrix_to_dmatrix2(spin_summed);
    }

    const std::vector<OhOperation> &oh_operations() {
        static const std::vector<OhOperation> operations = [] {
            std::vector<OhOperation> result;
            result.reserve(48);
            std::array<int, 3> permutation{ 0, 1, 2 };
            do {
                for (int mask = 0; mask < 8; ++mask) {
                    result.push_back({
                        permutation,
                        { (mask & 1) ? -1 : 1,
                          (mask & 2) ? -1 : 1,
                          (mask & 4) ? -1 : 1 }
                        });
                }
            } while (std::next_permutation(permutation.begin(), permutation.end()));
            return result;
            }();
        return operations;
    }

    std::vector<std::array<int, 3>> cartesian_exponents(const int l) {
        const int size = cartesian_shell_size(l);
        // constants::type_vector currently contains Cartesian components through h.
        const int first_type = l * (l + 1) * (l + 2) / 6 + 1;
        std::vector<std::array<int, 3>> result(size);
        for (int component = 0; component < size; ++component) {
            int exponent[3];
            constants::type2vector(first_type + component, exponent);
            err_checkf(exponent[0] >= 0 && exponent[0] + exponent[1] + exponent[2] == l,
                "Cartesian O_h symmetrization does not support angular momentum l=" +
                std::to_string(l) + ".", std::cout);
            result[component] = { exponent[0], exponent[1], exponent[2] };
        }
        return result;
    }

    std::vector<std::array<int, 3>> libcint_cartesian_exponents(const int l) {
        std::vector<std::array<int, 3>> result;
        result.reserve(cartesian_shell_size(l));
        for (int lx = l; lx >= 0; --lx) {
            for (int ly = l - lx; ly >= 0; --ly)
                result.push_back({ lx, ly, l - lx - ly });
        }
        return result;
    }

    CartesianTransform cartesian_transform_from_exponents(
        const std::vector<std::array<int, 3>> &exponents,
        const OhOperation &operation) {
        CartesianTransform transform{ ivec(exponents.size()), ivec(exponents.size(), 1) };

        for (int source = 0; source < static_cast<int>(exponents.size()); ++source) {
            std::array<int, 3> transformed{ 0, 0, 0 };
            int sign = 1;
            for (int axis = 0; axis < 3; ++axis) {
                transformed[operation.permutation[axis]] = exponents[source][axis];
                if (operation.signs[axis] < 0 && (exponents[source][axis] & 1))
                    sign = -sign;
            }

            const auto destination = std::find(exponents.begin(), exponents.end(), transformed);
            err_checkf(destination != exponents.end(),
                "O_h operation produced an unknown Cartesian component.", std::cout);
            transform.destination[source] = static_cast<int>(destination - exponents.begin());
            transform.signs[source] = sign;
        }
        return transform;
    }

    dMatrix2 cartesian_transform_matrix(const int l, const OhOperation &operation,
        const bool libcint_order = false) {
        const auto exponents = libcint_order
            ? libcint_cartesian_exponents(l)
            : cartesian_exponents(l);
        const auto sparse = cartesian_transform_from_exponents(exponents, operation);
        const int size = cartesian_shell_size(l);
        dMatrix2 result(size, size);
        std::fill(result.container().begin(), result.container().end(), 0.0);
        for (int source = 0; source < size; ++source)
            result(sparse.destination[source], source) = sparse.signs[source];
        return result;
    }

    dMatrix2 spherical_transform_matrix(const int l, const OhOperation &operation) {
        const dMatrix2 cart_to_spherical = cart2sph(l, true);
        const dMatrix2 cart_transform = cartesian_transform_matrix(l, operation, true);
        const int n_cart = cartesian_shell_size(l);
        const int n_spherical = 2 * l + 1;

        // Columns of C contain the real spherical functions in the Cartesian basis.
        // Solve T_cart C = C T_sph using the left pseudoinverse of C.  This keeps
        // libcint's real-spherical ordering and phase convention authoritative.
        dMatrix2 gram(n_spherical, n_spherical);
        for (int i = 0; i < n_spherical; ++i)
            for (int j = 0; j < n_spherical; ++j)
                for (int cart = 0; cart < n_cart; ++cart)
                    gram(i, j) += cart_to_spherical(cart, i) * cart_to_spherical(cart, j);
        const dMatrix2 inverse_gram = LAPACKE_invert(gram, 1E-12);

        dMatrix2 transformed_cart(n_cart, n_spherical);
        for (int i = 0; i < n_cart; ++i)
            for (int j = 0; j < n_spherical; ++j)
                for (int cart = 0; cart < n_cart; ++cart)
                    transformed_cart(i, j) += cart_transform(i, cart) * cart_to_spherical(cart, j);

        dMatrix2 result(n_spherical, n_spherical);
        for (int i = 0; i < n_spherical; ++i) {
            for (int j = 0; j < n_spherical; ++j) {
                for (int k = 0; k < n_spherical; ++k) {
                    for (int cart = 0; cart < n_cart; ++cart) {
                        result(i, j) += inverse_gram(i, k) *
                            cart_to_spherical(cart, k) * transformed_cart(cart, j);
                    }
                }
                if (std::abs(result(i, j)) < 1E-12)
                    result(i, j) = 0.0;
            }
        }
        return result;
    }
} // namespace

void symmetrize_atomic_matrix_oh(dMatrix2 &matrix, const ivec &shell_angular_momenta,
    const bool spherical) {
    const bool square = matrix.extent(0) == matrix.extent(1);
    err_checkf(square, "Cannot O_h-symmetrize a non-square atomic matrix.", std::cout);
    if (!square)
        return;

    ivec shell_offsets(shell_angular_momenta.size() + 1, 0);
    for (int shell = 0; shell < static_cast<int>(shell_angular_momenta.size()); ++shell) {
        const int l = shell_angular_momenta[shell];
        const bool supported = l >= 0 && l <= 5;
        err_checkf(supported,
            "O_h symmetrization supports shells from s through h.", std::cout);
        if (!supported)
            return;
        shell_offsets[shell + 1] = shell_offsets[shell] +
            (spherical ? 2 * l + 1 : cartesian_shell_size(l));
    }

    const bool correct_size = shell_offsets.back() == static_cast<int>(matrix.extent(0));
    err_checkf(correct_size,
        "Atomic matrix size does not match its Cartesian shell description.", std::cout);
    if (!correct_size)
        return;

    dMatrix2 symmetrized(matrix.extent(0), matrix.extent(1));
    std::fill(symmetrized.container().begin(), symmetrized.container().end(), 0.0);

    const auto &operations = oh_operations();
    for (const auto &operation : operations) {
        if (spherical) {
            std::vector<dMatrix2> transforms;
            transforms.reserve(shell_angular_momenta.size());
            for (const int l : shell_angular_momenta)
                transforms.push_back(spherical_transform_matrix(l, operation));

            for (int shell_a = 0; shell_a < static_cast<int>(shell_angular_momenta.size()); ++shell_a) {
                const int first_a = shell_offsets[shell_a];
                const int size_a = 2 * shell_angular_momenta[shell_a] + 1;
                for (int shell_b = 0; shell_b < static_cast<int>(shell_angular_momenta.size()); ++shell_b) {
                    const int first_b = shell_offsets[shell_b];
                    const int size_b = 2 * shell_angular_momenta[shell_b] + 1;
                    for (int transformed_a = 0; transformed_a < size_a; ++transformed_a) {
                        for (int transformed_b = 0; transformed_b < size_b; ++transformed_b) {
                            for (int a = 0; a < size_a; ++a) {
                                for (int b = 0; b < size_b; ++b) {
                                    symmetrized(first_a + transformed_a, first_b + transformed_b) +=
                                        transforms[shell_a](transformed_a, a) *
                                        matrix(first_a + a, first_b + b) *
                                        transforms[shell_b](transformed_b, b);
                                }
                            }
                        }
                    }
                }
            }
            continue;
        }

        std::vector<CartesianTransform> transforms;
        transforms.reserve(shell_angular_momenta.size());
        for (const int l : shell_angular_momenta)
            transforms.push_back(cartesian_transform_from_exponents(
                cartesian_exponents(l), operation));

        for (int shell_a = 0; shell_a < static_cast<int>(shell_angular_momenta.size()); ++shell_a) {
            const int first_a = shell_offsets[shell_a];
            const int size_a = cartesian_shell_size(shell_angular_momenta[shell_a]);
            const auto &transform_a = transforms[shell_a];
            for (int shell_b = 0; shell_b < static_cast<int>(shell_angular_momenta.size()); ++shell_b) {
                const int first_b = shell_offsets[shell_b];
                const int size_b = cartesian_shell_size(shell_angular_momenta[shell_b]);
                const auto &transform_b = transforms[shell_b];
                for (int a = 0; a < size_a; ++a) {
                    const int transformed_a = first_a + transform_a.destination[a];
                    for (int b = 0; b < size_b; ++b) {
                        const int transformed_b = first_b + transform_b.destination[b];
                        symmetrized(transformed_a, transformed_b) +=
                            transform_a.signs[a] * transform_b.signs[b] *
                            matrix(first_a + a, first_b + b);
                    }
                }
            }
        }
    }

    const double inverse_order = 1.0 / static_cast<double>(operations.size());
    for (auto &value : symmetrized.container())
        value *= inverse_order;
    matrix = std::move(symmetrized);
}

#ifdef NSA2DEBUG
void print_dmatrix2(const dMatrix2 &EVC2, const std::string name) {
    std::cout << std::endl << name << ":\n";
    for (int i = 0; i < EVC2.extent(0); i++) {
        for (int j = 0; j < EVC2.extent(1); j++)
            std::cout << std::setw(14) << std::setprecision(8) << std::fixed << EVC2(i, j) << " ";
        std::cout << std::endl;
    }
}
#endif


int compute_dens(WFN &wavy, bool debug, int *np, double *origin, double *gvector, double *incr, std::string &outname, bool rho, bool rdg, bool eli, bool lap) {
    properties_options opts;
    opts.lap = lap;
    opts.eli = eli;
    opts.rdg = rdg;
    opts.rho = rho;
    opts.resolution = *incr;
    opts.NbSteps = { np[0], np[1], np[2] };

    std::vector<cube> cubes = { {opts.NbSteps, wavy.get_ncen(), opts.rho || opts.eli || opts.lap},
        {},
        {opts.NbSteps, wavy.get_ncen(), opts.rdg},
        {},
        {opts.NbSteps, wavy.get_ncen(), opts.eli},
        {opts.NbSteps, wavy.get_ncen(), opts.lap},
        {},
        {},
        {},
        {},
        {},
        {},
        {} };

    std::string Oname = outname;
    std::vector<int> ntyp;
    for (int i = 0; i < 3; i++) {
        if (debug) {
            std::cout << "gvector before: ";
            for (int j = 0; j < 3; j++) std::cout << gvector[j + 3 * i] << " ";
            std::cout << std::endl;
        }
        for (int j = 0; j < 3; j++)
            gvector[j + 3 * i] = incr[i] * gvector[j + 3 * i];
        if (debug) {
            for (int j = 0; j < 3; j++) std::cout << gvector[j + 3 * i] << " ";
            std::cout << std::endl;
        }
    }
    for (int i = 0; i < 3; i++) {
        cubes[cube_type::Rho].set_origin(i, origin[i]);
        cubes[cube_type::RDG].set_origin(i, origin[i]);
        cubes[cube_type::Eli].set_origin(i, origin[i]);
        cubes[cube_type::Lap].set_origin(i, origin[i]);
        for (int j = 0; j < 3; j++) {
            cubes[cube_type::Rho].set_vector(i, j, gvector[j + 3 * i]);
            cubes[cube_type::RDG].set_vector(i, j, gvector[j + 3 * i]);
            cubes[cube_type::Eli].set_vector(i, j, gvector[j + 3 * i]);
            cubes[cube_type::Lap].set_vector(i, j, gvector[j + 3 * i]);
        }
    }
    cubes[cube_type::Rho].set_path(Oname + "_rho.cube");
    cubes[cube_type::RDG].set_path(Oname + "_rdg.cube");
    cubes[cube_type::Eli].set_path(Oname + "_eli.cube");
    cubes[cube_type::Lap].set_path(Oname + "_lap.cube");

    //opt.NbAtoms[0]=wavy.get_ncen();
    std::cout << "\n   .      ___________________________________________________________      .\n";
    std::cout << "  *.                                                                      *.\n";
    std::cout << "  *.  Wavefunction              : " << std::setw(20) << wavy.get_path().filename() << " / " << std::setw(5) << wavy.get_ncen() << " atoms      *.\n";
    std::cout << "  *.  OutPut filename Prefix    : " << std::setw(40) << Oname << "*.\n";
    std::cout << "  *.                                                                      *.\n";
    std::cout << "  *.  gridBox Min               : " << std::setw(11) << std::setprecision(6) << origin[0] << std::setw(12) << origin[1] << origin[2] << "     *.\n";
    std::cout << "  *.  gridBox Max               : " << std::setw(11) << std::setprecision(6) << (origin[0] + gvector[0] * np[0] + gvector[3] * np[1] + gvector[6] * np[2]) << std::setw(12) << (origin[1] + gvector[1] * np[0] + gvector[4] * np[1] + gvector[7] * np[2]) << (origin[2] + gvector[2] * np[0] + gvector[5] * np[1] + gvector[8] * np[2]) << "     *.\n";
    std::cout << "  *.  Increments(bohr)          : " << std::setw(11) << std::setprecision(6) << array_length(d3{ gvector[0], gvector[1], gvector[2] }) << std::setw(12) << array_length(d3{ gvector[3], gvector[4], gvector[5] }) << std::setw(12) << array_length(d3{ gvector[6], gvector[7], gvector[8] }) << "     *.\n";
    std::cout << "  *.  NbSteps                   : " << std::setw(11) << opts.NbSteps[0] << std::setw(12) << opts.NbSteps[1] << std::setw(12) << opts.NbSteps[2] << "     *.\n";
    std::cout << "  *.                                                                      *.\n";
    std::cout << "  *.  Number of primitives      :     " << std::setw(5) << wavy.get_nex() << "                               *.\n";
    std::cout << "  *.  Number of MOs             :       " << std::setw(3) << wavy.get_nmo() << "                               *.\n";

    cube dummy({ 0, 0, 0 });
    Calc_Prop(
        cubes,
        wavy,
        20.0,
        std::cout,
        false,
        false
    );

    std::cout << "  *.                                                                      *.\n";
    std::cout << "  *.  Writing .cube files ...                                             *.\n";
    if (rho && !rdg) {
        cubes[cube_type::Rho].set_path(Oname + "_rho.cube");
        cubes[cube_type::Rho].write_file(true, true);
    }
    if (rdg) {
        cubes[cube_type::Rho].set_path(Oname + "_signed_rho.cube");
        cubes[cube_type::Rho].write_file(true);
        cubes[cube_type::Rho].set_path(Oname + "_rho.cube");
        cubes[cube_type::Rho].write_file(true, true);
        cubes[cube_type::RDG].write_file(true);
    }
    if (eli) cubes[cube_type::Eli].write_file(true);
    if (lap) cubes[cube_type::Lap].write_file(true);

    std::cout << "  *                                                                       *\n";
    std::cout << "          ___________________________________________________________\n";
    return 0;
};

bond do_bonds(WFN &wavy,
    int mode_sel, bool mode_leng, bool mode_res,
    double res[], bool cub, double boxsize[],
    int atom1, int atom2, int atom3,
    const bool &debug, const bool &bohr,
    int runnumber,
    bool rho, bool rdg, bool eli, bool lap) {
    bond results{ "","","","",false,false,false };
    std::string line("");
    std::vector <std::string> label;
    label.resize(3);
    double coords1[3], coords2[3], coords3[3];
    int na = 0;
    double z[3], x[3], y[3], help[3], size[3], s2[3];
    int np[3];
    double ang2bohr;
    if (bohr)ang2bohr = 0.52917720859;
    else ang2bohr = 1.0;
    na = wavy.get_ncen();
    if (atom1 <= 0 || atom2 <= 0 || atom3 <= 0 || mode_sel <= 0 || mode_sel > 4 || atom1 > na || atom2 > na || atom3 > na || atom1 == atom2 || atom2 == atom3 || atom1 == atom3)
    {
        std::cout << "Invalid selections of atoms or mode_sel, please try again!\n";
        return(results);
    }
    if (debug) std::cout << "No. of Atoms selected: " << atom1 << atom2 << atom2 << "\n";
    for (int i = 0; i < 3; i++)
        coords1[i] = wavy.get_atom_coordinate(atom1 - 1, i);
    for (int i = 0; i < 3; i++)
        coords2[i] = wavy.get_atom_coordinate(atom2 - 1, i);
    for (int i = 0; i < 3; i++)
        coords3[i] = wavy.get_atom_coordinate(atom3 - 1, i);
    label[0] = wavy.get_atom_label(atom1 - 1);
    label[1] = wavy.get_atom_label(atom2 - 1);
    label[2] = wavy.get_atom_label(atom3 - 1);
    if (debug)
    {
        std::cout << "The Atoms found corresponding to your selection are:\n";
        std::cout << label[0] << " " << coords1[0] << " " << coords1[1] << " " << coords1[2] << "\n";
        std::cout << label[1] << " " << coords2[0] << " " << coords2[1] << " " << coords2[2] << "\n";
        std::cout << label[2] << " " << coords3[0] << " " << coords3[1] << " " << coords3[2] << "\n";
    }
    double znorm = sqrt((coords2[0] - coords1[0]) * (coords2[0] - coords1[0]) + (coords2[1] - coords1[1]) * (coords2[1] - coords1[1]) + (coords2[2] - coords1[2]) * (coords2[2] - coords1[2]));
    z[0] = (coords2[0] - coords1[0]) / znorm;
    z[1] = (coords2[1] - coords1[1]) / znorm;
    z[2] = (coords2[2] - coords1[2]) / znorm;
    help[0] = coords3[0] - coords2[0];
    help[1] = coords3[1] - coords2[1];
    help[2] = coords3[2] - coords2[2];
    y[0] = z[1] * help[2] - z[2] * help[1];
    y[1] = z[2] * help[0] - z[0] * help[2];
    y[2] = z[0] * help[1] - z[1] * help[0];
    double ynorm = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
    y[0] = y[0] / ynorm;
    y[1] = y[1] / ynorm;
    y[2] = y[2] / ynorm;
    x[0] = (z[1] * y[2] - z[2] * y[1]);
    x[1] = (z[2] * y[0] - z[0] * y[2]);
    x[2] = (z[0] * y[1] - z[1] * y[0]);
    z[0] = -z[0];
    z[1] = -z[1];
    z[2] = -z[2];

    double hnorm = 0.0;
    double ringhelp[3] = {};
    if (mode_sel == 2 || mode_sel == 1) hnorm = znorm;
    else if (mode_sel == 3) hnorm = sqrt((coords3[0] - coords1[0]) * (coords3[0] - coords1[0]) + (coords3[1] - coords1[1]) * (coords3[1] - coords1[1]) + (coords3[2] - coords1[2]) * (coords3[2] - coords1[2]));
    else if (mode_sel == 4)
    {
        for (int r = 0; r < 3; r++) ringhelp[r] = (coords3[r] + coords1[r] + coords2[r]) / 3;
        hnorm = (sqrt((coords3[0] - ringhelp[0]) * (coords3[0] - ringhelp[0]) + (coords3[1] - ringhelp[1]) * (coords3[1] - ringhelp[1]) + (coords3[2] - ringhelp[2]) * (coords3[2] - ringhelp[2]))
            + sqrt((coords1[0] - ringhelp[0]) * (coords1[0] - ringhelp[0]) + (coords1[1] - ringhelp[1]) * (coords1[1] - ringhelp[1]) + (coords1[2] - ringhelp[2]) * (coords1[2] - ringhelp[2]))
            + sqrt((coords2[0] - ringhelp[0]) * (coords2[0] - ringhelp[0]) + (coords2[1] - ringhelp[1]) * (coords2[1] - ringhelp[1]) + (coords2[2] - ringhelp[2]) * (coords2[2] - ringhelp[2]))) / 3;
    }
    if (debug)
    {
        std::cout << "Your three vectors are:\n";
        std::cout << "X= " << x[0] << " " << x[1] << " " << x[2] << "\n";
        std::cout << "Y= " << y[0] << " " << y[1] << " " << y[2] << "\n";
        std::cout << "Z= " << z[0] << " " << z[1] << " " << z[2] << "\n";
        std::cout << "From here on mode_leng and cub are used\n";
    }
    for (int r = 0; r < 3; r++)
    {
        if (res[r] < 0)
        {
            std::cout << "Wrong input in res! Try again!\n";
            return(results);
        }
        if (boxsize[r] < 0)
        {
            std::cout << "Wrong input for box scaling! Try again!\n";
            return(results);
        }
        else if (boxsize[r] >= 50)
        {
            std::cout << "Come on, be realistic!\nTry Again!\n";
            return(results);
        }
    }
    if (!mode_res)
    {
        if (debug) std::cout << "Determining resolution\n";
        if (mode_leng)
        {
            if (debug) std::cout << "mres=true; mleng=true\n";
            switch (mode_sel)
            {
            case 1:
                size[0] = ceil(9 * znorm / ang2bohr) / 10;
                break;
            case 2:
                size[0] = ceil(15 * znorm / ang2bohr) / 10;
                break;
            case 3:
                size[0] = ceil(15 * hnorm / ang2bohr) / 10;
                break;
            case 4:
                size[0] = ceil(30 * hnorm / ang2bohr) / 10;
                break;
            }
            for (int r = 0; r < 3; r++)
            {
                if (boxsize[r] != 0) size[r] = boxsize[r] / ang2bohr;
                else if (r > 0 || cub == true) size[r] = size[0];
                s2[r] = size[r] / 2;
            }
        }
        else
        {
            if (debug) std::cout << "mres=true; mleng=false\n";
            for (int r = 0; r < 3; r++)
            {
                if (cub && r > 0)
                {
                    if (debug) std::cout << "I'm making things cubic!\n";
                    s2[r] = s2[0];
                    continue;
                }
                switch (mode_sel)
                {
                case 1:
                case 2:
                    size[r] = znorm / ang2bohr;
                    break;
                case 3:
                case 4:
                    size[r] = hnorm / ang2bohr;
                    break;
                }
                if (boxsize[r] != 0) size[r] = size[r] * boxsize[r];
                else {
                    switch (mode_sel)
                    {
                    case 1:
                        size[r] = ceil(9 * znorm / ang2bohr) / 10;
                        break;
                    case 2:
                        size[r] = ceil(15 * znorm / ang2bohr) / 10;
                        break;
                    case 3:
                        size[r] = ceil(15 * hnorm / ang2bohr) / 10;
                        break;
                    case 4:
                        size[r] = ceil(30 * hnorm / ang2bohr) / 10;
                        break;
                    }
                }
                s2[r] = size[r] / 2;
            }
        }
        for (int r = 0; r < 3; r++) np[r] = (int)round((2 * s2[r]) / res[r]) + 1;
    }
    else
    {
        if (debug) std::cout << "Boxsize is used\n";
        for (int r = 0; r < 3; r++)
        {
            if (cub) np[r] = (int)res[0];
            else np[r] = (int)res[r];
        }
        if (mode_leng)
        {
            if (debug) std::cout << "mode_leng=true; using boxsize\n";
            switch (mode_sel)
            {
            case 1:
                size[0] = ceil(9 * znorm / ang2bohr) / 10;
                break;
            case 2:
                size[0] = ceil(15 * znorm / ang2bohr) / 10;
                break;
            case 3:
                size[0] = ceil(15 * hnorm / ang2bohr) / 10;
                break;
            case 4:
                size[0] = ceil(30 * hnorm / ang2bohr) / 10;
                break;
            }
            for (int r = 0; r < 3; r++)
            {
                if (debug) std::cout << r + 1 << ". Axis:";
                if (boxsize[r] != 0) size[r] = boxsize[r];
                else if (r > 0 || cub == true) size[r] = size[0];
                s2[r] = size[r] / 2;
            }
        }
        else
        {
            if (debug) std::cout << "Enter multiplicator for length (z,y,x)\n";
            for (int r = 0; r < 3; r++)
            {
                if (cub == 1 && r > 0)
                {
                    if (debug) std::cout << "I'm making things cubic!\n";
                    s2[r] = s2[0];
                    size[r] = size[0];
                    continue;
                }
                if (debug) std::cout << r + 1 << ". Axis:";
                switch (mode_sel)
                {
                case 1:
                case 2:
                    size[r] = znorm / ang2bohr;
                    break;
                case 3:
                case 4:
                    size[r] = hnorm / ang2bohr;
                    break;
                }
                if (boxsize[r] != 0) size[r] = size[r] * boxsize[r];
                else {
                    switch (mode_sel)
                    {
                    case 1:
                        size[r] = ceil(9 * znorm / ang2bohr) / 10;
                        break;
                    case 2:
                        size[r] = ceil(15 * znorm / ang2bohr) / 10;
                        break;
                    case 3:
                        size[r] = ceil(15 * hnorm / ang2bohr) / 10;
                        break;
                    case 4:
                        size[r] = ceil(30 * hnorm / ang2bohr) / 10;
                        break;
                    }
                }
                s2[r] = size[r] / 2;
            }
        }
    }
    double o[3] = {};
    if (debug) std::cout << "This is the origin:\n";
    for (int d = 0; d < 3; d++)
    {
        switch (mode_sel)
        {
        case 1:
            o[d] = coords2[d] - (s2[0] * x[d] + s2[1] * y[d] + s2[2] * z[d]) / ang2bohr;
            break;
        case 2:
            o[d] = (coords1[d] - (coords1[d] - coords2[d]) / 2) - (s2[0] * x[d] + s2[1] * y[d] + s2[2] * z[d]) / ang2bohr;
            break;
        case 3:
            o[d] = (coords3[d] - (coords3[d] - coords1[d]) / 2) - (s2[0] * x[d] + s2[1] * y[d] + s2[2] * z[d]) / ang2bohr;
            break;
        case 4:
            o[d] = ringhelp[d] - (s2[0] * x[d] + s2[1] * y[d] + s2[2] * z[d]) / ang2bohr;
            break;
        }
        if (debug) std::cout << o[d] << "\n";
    }
    std::string outname = { "" };
    outname += wavy.get_path().generic_string();
    outname += "_";
    outname += label[0];
    outname += std::to_string(atom1);
    outname += "_";
    outname += label[1];
    outname += std::to_string(atom2);
    outname += "_";
    outname += label[2];
    outname += std::to_string(atom3);
    outname += "_";
    outname += std::to_string(runnumber);
    double v[9];
    double incr[3];
    for (int i = 0; i < 3; i++) {
        v[i] = x[i];
        v[i + 3] = y[i];
        v[i + 6] = z[i];
        incr[i] = size[i] / np[i];
        //incr[i]=res[i];
    }
    if (compute_dens(wavy, debug, np, o, v, incr, outname, rho, rdg, eli, lap) == 0) {
        results.success = true;
        results.filename = outname;
    }
    return(results);
}

int autobonds(bool debug, WFN &wavy, const std::filesystem::path &inputfile, const bool &bohr) {
    char inputFile[2048] = "";
    if (!exists(inputfile))
    {
        std::cout << "No input file specified! I am going to make an example for you in input.example!\n";
        std::ofstream example("input.example");
        example << "!COMMENT LINES START WITH ! AND CAN ONLY BE ABOVE THE FIRST SWITCHES AND NUMBERS!\n!First row contains the following switches: rho (dens), Reduced Density Gradient (RDG), ELI-d (eli) and Laplacian of rho (lap)\n!Following rows each contain a bond you want to investigate.\n!The key to read these numbers is:\n!orientation_selection(1-4)\n!      1=atom1 centered\n!      2=bond atom1 atom2\n!      3=bond atom1 and atom3\n!      4=ring centroid of the three atoms\n!length selection(0/1)\n!      0=box will contain multiplicator of bondlength between atom1 and atom2\n!      1=box will contain length in angstrom (?)\n!resolution selection(0/1)\n!      0=res will contain number of gridpoints\n!      1=res will contain distance between gridpoints\n!res1 res2 res3\n!      based on selection above resolution in x,y,z direction (either gridpoints or distance between points)\n!cube selection (0/1)\n!      0= all selections are considered\n!      1= all selections in x-direction will be applied in the y and z direction, as well, making it a cube\n!box1 box2 box3\n!      based on selection above size of the box/cube (either multiplicator of bondlength or absolute length)\n!atom1 atom2 atom3\n!      the atoms in your wavefunction file (counting from 1) to be used as references\n! BELOW THE INPUT SECTION STARTS\n!rho, rdg, eli, lap\n!mode_sel(int) mode_leng(bool) mode_res(bool) res1(double) res2(double) res3(double) cube(bool) boxsize1(float) boxsize2(float) boxsize3(float) atom1(int) atom2(int) atom3(int)\n 1    1    1    1\n            2             1               1              20.0         20.0         20.0         0          5.0             5.1             5.2             1          2          3\n";
        example.close();
        return 0;
    }
    std::ifstream input(inputfile.c_str());
    if (!input.good())
    {
        std::cout << inputFile << " does not exist or is not readable!\n";
        return 0;
    }
    input.seekg(0);
    std::string line("");
    getline(input, line);
    std::string comment("!");
    while (line.compare(0, 1, comment) == 0) getline(input, line);
    int rho, rdg, eli, lap;
    if (line.length() < 10) return 0;
    else
    {
        std::istringstream iss(line);
        if (!(iss >> rho >> rdg >> eli >> lap)) {
            std::cerr << "Error parsing line for rho, rdg, eli, and lap values." << std::endl;
            return 0; // Handle the error appropriately
        }
    }
    int errorcount = 0, runnumber = 0;
    do {
        getline(input, line);
        if (line.length() < 10) continue;
        runnumber++;
        int sel = 0, leng = 0, cube = 0, mres = 0, a1 = 0, a2 = 0, a3 = 0;
        double res[3], box[3];
        bool bleng = false, bres = false, bcube = false;
        if (line.length() > 1)
        {
            std::istringstream iss(line);
            if (!(iss >> sel >> leng >> mres >> res[0] >> res[1] >> res[2] >> cube >> box[0] >> box[1] >> box[2] >> a1 >> a2 >> a3)) {
                std::cerr << "Error parsing line for bond parameters." << std::endl;
                continue; // Skip this line and move to the next
            }
        }
        if (leng == 1) { bleng = true; if (debug) std::cout << "leng=true\n"; }
        else { bleng = false; if (debug) std::cout << "leng=false\n"; }
        if (cube == 1) { bcube = true; if (debug) std::cout << "cube=true\n"; }
        else { bcube = false; if (debug) std::cout << "cube=false\n"; }
        if (mres == 1) { bres = true; if (debug) std::cout << "mres=true\n"; }
        else { bres = false; if (debug) std::cout << "mres=false\n"; }
        if (debug) std::cout << "running calculations for line " << runnumber << " of the input file:\n\n";
        bond work{ "","","","",false,false,false };
        work = do_bonds(wavy, sel, bleng, bres, res, bcube, box, a1, a2, a3, debug, bohr, runnumber, rho == 1, rdg == 1, eli == 1, lap == 1);
        if (!work.success)
        {
            std::cout << "!!!!!!!!!!!!problem somewhere during the calculations, see messages above!!!!!\n";
            errorcount++;
        }
        else {
            if (rho == 1) wavy.push_back_cube(work.filename + "_rho.cube", false, false);
            if (rdg == 1) wavy.push_back_cube(work.filename + "_rdg.cube", false, false);
            if (eli == 1) wavy.push_back_cube(work.filename + "_eli.cube", false, false);
            if (lap == 1) wavy.push_back_cube(work.filename + "_lap.cube", false, false);
            if (rho == 1 && rdg == 1) wavy.push_back_cube(work.filename + "_signed_rho.cube", false, false);
        }
    } while (!input.eof());
    std::cout << "\n  *   Finished all calculations! " << runnumber - errorcount << " out of " << runnumber << " were successful!            *\n";
    return 1;
}

std::vector<std::pair<int, int>> get_bonded_atom_pairs(const WFN &wavy) {
    std::vector<std::pair<int, int>> bonds;
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        for (int j = i + 1; j < wavy.get_ncen(); j++)
        {
            const double distance = array_length(wavy.get_atom_pos(i), wavy.get_atom_pos(j));
            const double svdW = constants::ang2bohr(constants::covalent_radii[wavy.get_atom_charge(i)] + constants::covalent_radii[wavy.get_atom_charge(j)]);
            if (distance < 1.25 * svdW)
            {
#ifdef NSA2DEBUG
                std::cout << "Bond between " << i << " (" << wavy.get_atom_charge(i) << ") and " << j << " (" << wavy.get_atom_charge(j) << ") with distance " << distance << " and svdW " << svdW << std::endl;
#endif // NSA2DEBUG
                bonds.push_back(std::make_pair(i, j));
            }
        }
    }
    return bonds;
}


//Assuming square matrices
vec change_basis_sq(const vec &in, const vec &transformation, int size) {

    // new = transform^T * in * tranform
    vec result(size * size);
    vec temp(size * size);
    // first we do temp = t^T * i
    cblas_dgemm(CblasRowMajor,
        CblasTrans, CblasNoTrans,
        size, size, size,
        1.0,
        transformation.data(), size,
        in.data(), size,
        0.0,
        temp.data(), size);

    //Then we do res = temp * t
    cblas_dgemm(CblasRowMajor,
        CblasNoTrans, CblasNoTrans,
        size, size, size,
        1.0,
        temp.data(), size,
        transformation.data(), size,
        0.0,
        result.data(), size);
    return result;
}

//For non-square transformation matrices
//performs res = trans * in * trans^T in case of forward
//and res = trans^T * int * trans in case of !forward
dMatrix2 change_basis_general(const dMatrix2 &in, const dMatrix2 &transformation, bool forward = true) {
    //Checks should be handled by dot itself
    //int a, b, c, d;
    //c = transformation.extent(1); //cols of t
    //d = transformation.extent(0); //rows of t
    //a = in.extent(1); //cols of in
    //b = in.extent(0); //rows of in
    //err_checkf(a == b, "Input matrix must be square.", std::cout);
    //err_checkf(c == b, "Incompatible matrix dimensions for basis change.", std::cout);
    if (!forward) {
        dMatrix2 temp = dot<dMatrix2>(transformation, in, false, false); // temp = t * in
        dMatrix2 result = dot<dMatrix2>(temp, transformation, false, true); // res = temp * t^T
        return result;
    }
    else {
        dMatrix2 temp = dot<dMatrix2>(transformation, in, true, false); // temp = t^T * in
        dMatrix2 result = dot<dMatrix2>(temp, transformation, false, false); // res = temp * t
        return result;
    }
}

/**
 * Calculates Atomic Natural Orbitals for a specific atom.
 *
 * @param D_full      Pointer to the full Density Matrix (N_basis x N_basis)
 * @param S_full      Pointer to the full Overlap Matrix (N_basis x N_basis)
 * @param full_stride The leading dimension of the full matrices (usually N_basis)
 * @param atom_indices A vector containing the indices (0-based) of the basis functions for this atom
 * @return NAOResult containing sorted occupancies and coefficients
 */
Roby_information::NAOResult Roby_information::calculateAtomicNAO(const dMatrix2 &D_full,
    const dMatrix2 &S_full,
    const std::vector<int> &atom_indices,
    const ivec &shell_angular_momenta,
    const bool spherical,
    const double occupancy_cutoff,
    const int leading_orbitals_to_skip) {

    err_checkf(D_full.extent(0) == D_full.extent(1), "Density matrix D must be square.", std::cout);
    err_checkf(S_full.extent(0) == S_full.extent(1), "Overlap matrix S must be square.", std::cout);
    err_checkf(D_full.extent(0) == S_full.extent(0), "Density and Overlap matrices must be of the same size.", std::cout);

    const int n = static_cast<int>(atom_indices.size());

    // 1. Memory Allocation for Submatrices
    // Using flat std::vectors to ensure contiguous memory for MKL
    vec D_sub(n * n, 0.0); // called P in tonto
    vec S_sub(n * n, 0.0); // called S in tonto
    get_submatrices(D_full, S_full, D_sub, S_sub, atom_indices);

    // Tonto's atomic spherical averaging applies the local point group to the
    // atom-centred density before the ANOs are constructed.  O_h is used here
    // because it averages all Cartesian directions without mixing radial shells.
    if (!shell_angular_momenta.empty()) {
        dMatrix2 atomic_density = reshape<dMatrix2>(D_sub, Shape2D(n, n));
        symmetrize_atomic_matrix_oh(atomic_density, shell_angular_momenta, spherical);
        D_sub = atomic_density.container();
    }
    vec Rho(n * n);        // To store target density

    vec V = S_sub;
    vec W(n);
    // make V = Sqrt(S)
    const vec Temp = mat_sqrt(V, W);

#ifdef NSA2DEBUG
    print_dmatrix2(reshape<dMatrix2>(Temp, Shape2D(n, n)), "projection matrix V");
#endif

    const vec X = change_basis_sq(D_sub, Temp, n);
#ifdef NSA2DEBUG
    print_dmatrix2(reshape<dMatrix2>(X, Shape2D(n, n)), "projected density X");
#endif

    vec occu(n, 0);
    vec P = X;
#ifdef NSA2DEBUG
    err_checkf(isSymmetricViaEigenvalues<vec>(P, n), "Transformed matrix not symmetric!", std::cout);
#endif
    make_Eigenvalues(P, occu);
#ifdef NSA2DEBUG
    std::cout << "Eigenvalues of projected density P:\n";
    for (int i = 0; i < n; i++) {
        std::cout << std::setw(14) << std::setprecision(8) << std::fixed << W[i] << " ";
    }
    print_dmatrix2(reshape<dMatrix2>(P, Shape2D(n, n)), "Projected density P");
#endif

    for (int i = 0; i < n; i++)
        W[i] = abs(W[i]) < 1E-10 ? 0.0 : 1.0 / W[i];

    vec Temp2(n * n, 0.0);
    double *T;
    int in, jn;
    for (int i = 0; i < n; i++) {
        in = i * n;
        for (int j = 0; j < n; j++) {
            jn = j * n;
            T = &Temp2[in + j];
            for (int k = 0; k < n; k++)
                *T += V[in + k] * W[k] * V[jn + k];
        }
    }

    V.clear();
    W.clear();

#ifdef NSA2DEBUG
    print_dmatrix2(reshape<dMatrix2>(Temp2, Shape2D(n, n)), "Back projection S^-0.5");
#endif

    cblas_dgemm(CblasRowMajor,
        CblasNoTrans, CblasNoTrans,
        n, n, n,
        1.0,
        Temp2.data(), n,
        P.data(), n,
        0.0,
        Rho.data(), n);

#ifdef NSA2DEBUG
    print_dmatrix2(reshape<dMatrix2>(Rho, Shape2D(n, n)), "resulting NAO");
#endif

    NAOResult result;
    result.eigenvalues = occu;
    result.eigenvectors = Rho; // Rho now contains eigenvectors
    result.matrix_elements = atom_indices;
    result.sub_DM = D_sub;
    result.sub_OM = S_sub;

    // Create an index vector to sort descending
    ivec idx(n);
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(), [&](int i1, int i2) {
        return result.eigenvalues[i1] > result.eigenvalues[i2];
        });

    // Reorder based on sorted indices
    vec sorted_evals;
    vec sorted_evecs;
    sorted_evecs.reserve(n * n);

    const int skip_orbitals = std::clamp(leading_orbitals_to_skip, 0, n);
    for (int i = 0; i < n; i++) {
        if (i < skip_orbitals)
            continue;
        int original_idx = idx[i];
        if (occupancy_cutoff >= 0.0 && result.eigenvalues[original_idx] < occupancy_cutoff)
            continue;
        sorted_evals.emplace_back(result.eigenvalues[original_idx]);

        for (int row = 0; row < n; row++) {
            sorted_evecs.emplace_back(result.eigenvectors[row * n + original_idx]);
        }
    }

    result.eigenvalues = sorted_evals;
    result.eigenvectors = sorted_evecs;

    return result;
}

/**
 * Computes Final NAOs by symmetrically orthogonalizing the Pre-NAOs.
 *
 * @param C_PNAO    Global matrix of Pre-NAOs (N x N). Columns are orbitals.
 * @param S_AO      Original Overlap Matrix in AO basis (N x N).
 * @param n         Number of basis functions (N).
 * @return          Matrix of final NAOs (N x N), Columns are orbitals.
 */

 /* CURRENTLY NOT IN USE
 std::vector<double> orthogonalizePNAOs(const vec& C_PNAO,
     const vec& S_AO,
     int n) {

     vec S_PNAO(n * n);
     vec Temp(n * n); // Intermediate buffer
     vec C_NAO(n * n); // Result

     // 1. Compute Overlap in PNAO basis: S_PNAO = C_PNAO^T * S_AO * C_PNAO

     // Step A: Temp = S_AO * C_PNAO
     // S_AO is symmetric.
     cblas_dsymm(CblasRowMajor,
         CblasLeft, CblasUpper,
         n, n,
         1.0, S_AO.data(), n,
         C_PNAO.data(), n,
         0.0, Temp.data(), n);

     // Step B: S_PNAO = C_PNAO^T * Temp
     // C_PNAO is not symmetric, so we use dgemm.
     // Transpose the first matrix (C_PNAO^T).
     cblas_dgemm(CblasRowMajor,
         CblasTrans, CblasNoTrans,
         n, n, n,
         1.0, C_PNAO.data(), n,
         Temp.data(), n,
         0.0, S_PNAO.data(), n);

     // 2. Compute S_PNAO^(-1/2) using Eigendecomposition
     // S_PNAO is symmetric (and positive definite).

     vec W(n); // Eigenvalues
     // We can overwrite S_PNAO with eigenvectors to save memory,
     // but let's keep it clear. Copy S_PNAO to 'U' (Eigenvectors).
     vec U = S_PNAO;

     // LAPACKE_dsyevd: Computes all eigenvalues and eigenvectors
     err_checkf(LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, U.data(), n, W.data()) == 0, "Eigenvalue computation failed.", std::cout);

     // Construct S^-1/2 = U * Lambda^(-1/2) * U^T
     std::fill(S_PNAO.begin(), S_PNAO.end(), 0.0); // Reuse S_PNAO to store S^-1/2

     for (int k = 0; k < n; k++) {
         double scale = 1.0 / std::sqrt(W[k]);
         // Add contribution of k-th eigenvector: scale * (v_k * v_k^T)
         // v_k is the k-th ROW of U.
         const double* v_k = &U[k * n];

         // This is a rank-1 update (dger), but we can just sum manually or loop
         // Since we need the full matrix for the next multiplication
         for (int i = 0; i < n; i++) {
             for (int j = 0; j < n; j++) {
                 S_PNAO[i * n + j] += scale * v_k[i] * v_k[j];
             }
         }
     }

     // 3. Transform PNAOs to NAOs
     // C_NAO = C_PNAO * S^(-1/2)
     // C_PNAO (Columns are PNAOs) * S_inv_sqrt (transformation matrix)

     cblas_dgemm(CblasRowMajor,
         CblasNoTrans, CblasNoTrans,
         n, n, n,
         1.0, C_PNAO.data(), n,
         S_PNAO.data(), n, // This is now S^-1/2
         0.0,
         C_NAO.data(), n);

     return C_NAO;
 }
 */

double Roby_information::projection_matrix_and_expectation(const ivec &indices, const ivec &eigvals, const ivec &eigvecs, dMatrix2 *given_NAO, dMatrix2 *proj_out) {
    const int n = indices.size();
    //vec D_Sub(n * n, 0.0);
    vec S_Sub(n * n, 0.0);
    get_submatrix(overlap_matrix, S_Sub, indices);
    dMatrix2 S = reshape<dMatrix2>(S_Sub, Shape2D(n, n));
    int atom = -1;
    //dMatrix2 D = reshape<dMatrix2>(D_Sub, Shape2D(n, n));

    auto zero_projection = [&]() {
        dMatrix2 zero(n, n);
        for (int r = 0; r < n; ++r)
            for (int c = 0; c < n; ++c)
                zero(r, c) = 0.0;
        if (atom >= 0 && proj_out == nullptr) {
            projection_matrices.push_back(zero);
            overlap_matrices.push_back(S);
        }
        if (proj_out != nullptr)
            *proj_out = zero;
        return 0.0;
        };

    dMatrix2 NAOs;
    // When an explicit NAO matrix is provided together with row/col selectors,
    // use it directly — this takes priority over the full-system shortcut below.
    if (given_NAO != nullptr && !eigvals.empty() && !eigvecs.empty()) {
        const int n1 = eigvals.size();
        const int n2 = eigvecs.size();
        vec NAO_sub(n1 * n2);
        get_submatrix(*given_NAO, NAO_sub, eigvals, eigvecs);
        NAOs = reshape<dMatrix2>(NAO_sub, Shape2D(n1, n2));
    }
    else if (given_NAO == nullptr && eigvals.empty() && eigvecs.empty()
        && n == static_cast<int>(overlap_matrix.extent(0))) {
        err_checkf(static_cast<int>(total_NAOs.extent(1)) == n,
            "RGBI total NAO matrix column count (" + std::to_string(total_NAOs.extent(1)) +
            ") does not match the basis-function count (" + std::to_string(n) + ").",
            std::cout);
        NAOs = total_NAOs;
    }
    else {
        //TODO: assign subspace NAOs from NAOResults for a given atom
        for (auto NAO : this->NAOs) {
            // if matrix_elemnts of this NAO are identical to indices, resahpe NAO.eigenvectors to the correct shape
            if (NAO.matrix_elements == indices) {
                NAOs = reshape<dMatrix2>(NAO.eigenvectors, Shape2D(NAO.eigenvalues.size(), n));
                atom = NAO.atom_index;
                break;
            }
        }
    }
    //For atom groups we can use submatrices of total_NAOs
    if (NAOs.size() == 0) {
        const int n1 = eigvals.size();
        const int n2 = eigvecs.size();
        if (n1 == 0 || n2 == 0)
            return zero_projection();
        vec NAO_sub(n1 * n2);
        if (given_NAO == nullptr)
            given_NAO = &total_NAOs;
        get_submatrix(*given_NAO, NAO_sub, eigvals, eigvecs);
        NAOs = reshape<dMatrix2>(NAO_sub, Shape2D(n1, n2));
    }
#ifdef NSA2DEBUG
    print_dmatrix2(NAOs, "W in projection matrix making");
#endif

    auto X = change_basis_general(S, transpose(NAOs));
#ifdef NSA2DEBUG
    print_dmatrix2(X, "new Basis");
#endif

    auto Y = LAPACKE_invert(X);
#ifdef NSA2DEBUG
    print_dmatrix2(Y, "Pseudo inverse of Y");
#endif

    //making the projection matrix
    X = change_basis_general(Y, transpose(NAOs), false);
#ifdef NSA2DEBUG
    print_dmatrix2(X, "Backtransformed X");
#endif

    if (atom >= 0 && proj_out == nullptr) {
        projection_matrices.push_back(X);
        overlap_matrices.push_back(S);
    }
    if (proj_out != nullptr)
        *proj_out = X;

    dMatrix2 S_rect = get_rectangle(overlap_matrix, indices);
#ifdef NSA2DEBUG
    print_dmatrix2(S_rect, "S_rect:");
#endif

    //overlap projection
    auto W = change_basis_general(X, S_rect);
#ifdef NSA2DEBUG
    print_dmatrix2(W, "Overlap transformed W");
#endif

    const double expect = trace_product<double>(W, density_matrix);
    return expect;
}

double Roby_information::Roby_population_analysis(const ivec atoms) {
    ivec bf_indices;
    if (atoms.size() == 0) {
        const int n_basis_functions = static_cast<int>(overlap_matrix.extent(0));
        bf_indices.reserve(n_basis_functions);
        for (int index = 0; index < n_basis_functions; ++index)
            bf_indices.push_back(index);
    }
    else {
        bf_indices = atoms;
    }
    double P = projection_matrix_and_expectation(bf_indices);
    return P;
}

void Roby_information::computeAllAtomicNAOs(WFN &wavy, const bool symmetrize, const bool use_ano_basis) {
    const int N_atoms = wavy.get_ncen();
    const std::vector<atom> ats = wavy.get_atoms();
    NAOs.reserve(N_atoms);
    ano_fallback_atoms.clear();

    density_matrix = wavy.get_dm();

    Int_Params basis(wavy);
    vec S_full;
    if (wavy.get_d_f_switch())
        compute2C<Overlap2C_CRT>(basis, S_full);
    else
        compute2C<Overlap2C_SPH>(basis, S_full);

    overlap_matrix = reshape<dMatrix2>(S_full, Shape2D(density_matrix.extent(0), density_matrix.extent(1)));

#ifdef NSA2DEBUG
    print_dmatrix2(overlap_matrix, "Overlap matrix");
#endif

    //err_checkf()

    const double occupancy_cutoff = 0.17;

    int last_index = 0;
    ivec2 indices(wavy.get_ncen());
    for (auto &a : ats) {
        indices[a.get_nr() - 1].reserve(density_matrix.extent(0) / N_atoms); // Rough estimate
        ivec shell_angular_momenta;
        int current_shell = -1;
        int nr_indices = 0;
        std::vector<basis_set_entry> basis_set = a.get_basis_set();
        for (auto &bf : basis_set) {
            if (bf.get_shell() != current_shell) {
                current_shell++;
                const int l = wavy.get_origin() == e_origin::OCC
                    ? static_cast<int>(bf.get_type())
                    : static_cast<int>(bf.get_type()) - 1;
                err_checkf(l >= 0,
                    "Encountered an invalid shell angular momentum while building RGBI atomic NAOs.",
                    std::cout);
                shell_angular_momenta.push_back(l);
                nr_indices = atomic_shell_size(l, wavy.get_d_f_switch());
                if (wavy.get_origin() == e_origin::tonto && bf.get_type() == 3) {
                    // 2 - 4; 3 - 6; 5 - 6
                    swap_rows_cols_symm(overlap_matrix, last_index + 1, last_index + 3);
                    swap_rows_cols_symm(overlap_matrix, last_index + 2, last_index + 5);
                    swap_rows_cols_symm(overlap_matrix, last_index + 4, last_index + 5);
                }

                //if (wavy.get_origin() != e_origin::tonto && wavy.get_origin() != e_origin::OCC)
                //    bf.get_type() == 1 ? nr_indices = 1 : (bf.get_type() == 2 ? nr_indices = 3 : (bf.get_type() == 3 ? nr_indices = 5 : nr_indices = 7));
                //else {
                //    
                //    bf.get_type() == 1 ? nr_indices = 1 : (bf.get_type() == 2 ? nr_indices = 3 : (bf.get_type() == 3 ? nr_indices = 6 : nr_indices = 10));
                //    if (bf.get_type() == 3 && wavy.get_origin() == e_origin::tonto) {
                //        //2 - 4; 3 - 6; 5 - 6
                //        swap_rows_cols_symm(overlap_matrix, last_index + 1, last_index + 3);
                //        swap_rows_cols_symm(overlap_matrix, last_index + 2, last_index + 5);
                //        swap_rows_cols_symm(overlap_matrix, last_index + 4, last_index + 5);
                //    }
                //}
                for (int i = 0; i < nr_indices; i++) {
                    indices[a.get_nr() - 1].push_back(last_index);
                    last_index++;
                }
            }
        }

        // The GBW reader converts ORCA components to PySCF/libcint order and
        // groups an atom's shells by increasing angular momentum.  Mirror that
        // layout here so each symmetry block describes the corresponding DM rows.
        if (wavy.get_origin() == e_origin::gbw)
            std::stable_sort(shell_angular_momenta.begin(), shell_angular_momenta.end());

        const bool spherical = !wavy.get_d_f_switch();

        auto make_molecular_fallback = [&]() {
            auto fallback = calculateAtomicNAO(density_matrix, overlap_matrix,
                indices[a.get_nr() - 1],
                symmetrize ? shell_angular_momenta : ivec{},
                spherical,
                occupancy_cutoff);
            fallback.atom_index = a.get_nr() - 1;
            return fallback;
        };
        if (use_ano_basis) {
            try {
                const dMatrix2 atomic_density =
                    compute_tonto_style_atomic_density(a, wavy.get_origin(), wavy.get_d_f_switch());
                const int n_local = static_cast<int>(indices[a.get_nr() - 1].size());
                if (static_cast<int>(atomic_density.extent(0)) != n_local ||
                    static_cast<int>(atomic_density.extent(1)) != n_local) {
                    throw std::runtime_error(
                        "Atomic OCC density size does not match the loaded WFN basis layout.");
                }
                const int effective_atomic_number = a.get_charge() - a.get_ECP_electrons();
                const int core_orbitals_to_skip =
                    std::max(0, tonto_core_electron_count(effective_atomic_number) / 2);
                ivec local_indices(indices[a.get_nr() - 1].size());
                std::iota(local_indices.begin(), local_indices.end(), 0);
                vec S_sub(n_local * n_local, 0.0);
                get_submatrix(overlap_matrix, S_sub, indices[a.get_nr() - 1]);
                const dMatrix2 atomic_overlap = reshape<dMatrix2>(S_sub, Shape2D(n_local, n_local));
                auto ano = calculateAtomicNAO(atomic_density, atomic_overlap,
                    local_indices,
                    symmetrize ? shell_angular_momenta : ivec{},
                    spherical,
                    occupancy_cutoff,
                    core_orbitals_to_skip);
                ano.sub_OM = S_sub;
                ano.sub_DM = atomic_density.container();
                ano.matrix_elements = indices[a.get_nr() - 1];
                ano.atom_index = a.get_nr() - 1;
                NAOs.emplace_back(std::move(ano));
            }
            catch (const std::exception &e) {
                std::cout << "Warning: RGBI ANO state build failed for atom "
                    << a.get_nr() - 1 << " (" << a.get_label()
                    << "). Falling back to molecular local orbitals. Reason: "
                    << e.what() << '\n';
                ano_fallback_atoms.push_back(a.get_nr() - 1);
                NAOs.emplace_back(make_molecular_fallback());
            }
            catch (...) {
                std::cout << "Warning: RGBI ANO state build failed for atom "
                    << a.get_nr() - 1 << " (" << a.get_label()
                    << "). Falling back to molecular local orbitals.\n";
                ano_fallback_atoms.push_back(a.get_nr() - 1);
                NAOs.emplace_back(make_molecular_fallback());
            }
        }
        else {
            NAOs.emplace_back(make_molecular_fallback());
        }
        NAOs.back().atom_index = a.get_nr() - 1;
    }
#ifdef NSA2DEBUG
    print_dmatrix2(overlap_matrix, "Overlap matrix repaired");
#endif
}

ivec Roby_information::find_eigenvalue_pairs(const vec &eigvals, const double tolerance) {
    const int n = eigvals.size();
    ivec pairs(n, -1);
    for (int i = 0; i < n; i++) {
        if (pairs[i] >= 0)
            continue; // already paired
        for (int j = 0; j < n; j++) {
            if (pairs[j] >= 0)
                continue; // already paired
            if (abs(abs(eigvals[j]) - 1.0) < tolerance)
                continue; // skip values close ot +/- one, as they reside only on one atom
            if (abs(eigvals[i] + eigvals[j]) < tolerance) {
                pairs[i] = j;
                pairs[j] = i;
                break;
            }
        }
        if (pairs[i] == -1)
            pairs[i] = i;
    }
    return pairs;
}

void Roby_information::transform_Ionic_eigenvectors_to_Ionic_orbitals(
    dMatrix2 &EVC,
    const vec &eigvals,
    const ivec &pairs,
    const int index_a,
    const int index_b,
    const ivec &pair_matrix_indices)
{
    double fp, fm, fa, fb, s, c, s2;
    const int n_ab = EVC.extent(0);
    const int n_eigvals = eigvals.size();
    const int n_a = projection_matrices[index_a].extent(0);
    const int n_b = projection_matrices[index_b].extent(0);
    err_checkf(n_ab == n_a + n_b, "Inconsitent size in projection matrices?!", std::cout);
    err_checkf(n_eigvals == EVC.extent(1), "Inconsitency between EVC and eigvals", std::cout);

    dMatrix1 A(n_a), B(n_b);
    dMatrix1 EVC_column(n_ab);
    dMatrix2 PAS(n_a, n_ab), PBS(n_b, n_ab);

    const ivec *indices_a = nullptr;
    const ivec *indices_b = nullptr;
    for (const auto &NAO : NAOs) {
        if (NAO.atom_index == index_a)
            indices_a = &NAO.matrix_elements;
        if (NAO.atom_index == index_b)
            indices_b = &NAO.matrix_elements;
    }
    err_checkf(indices_a != nullptr, "No NAO data found for atom " + std::to_string(index_a), std::cout);
    err_checkf(indices_b != nullptr, "No NAO data found for atom " + std::to_string(index_b), std::cout);

    vec Sub_overlap(n_a * n_ab);
    get_submatrix(overlap_matrix, Sub_overlap, *indices_a, pair_matrix_indices);
    dMatrix2 Sa = reshape<dMatrix2>(Sub_overlap, Shape2D(n_a, n_ab));
    PAS = dot<dMatrix2>(projection_matrices[index_a], Sa, false, false);
    Sub_overlap.clear(); Sub_overlap.resize(n_b * n_ab);
    get_submatrix(overlap_matrix, Sub_overlap, *indices_b, pair_matrix_indices);
    dMatrix2 Sb = reshape<dMatrix2>(Sub_overlap, Shape2D(n_b, n_ab));
    PBS = dot<dMatrix2>(projection_matrices[index_b], Sb, false, false);
#ifdef NSA2DEBUG
    print_dmatrix2(PAS, "PAS");
    print_dmatrix2(PBS, "PBS");
#endif // NSA2DEBUG


    for (int i = 0; i < n_eigvals; i++) {
        if (pairs[i] < 0) continue;
        if (pairs[i] == i) continue;
        if (eigvals[i] < eigvals[pairs[i]]) continue;
#ifdef NSA2DEBUG
        std::cout << "Doing i=" << i << std::endl;
#endif
        for (int a = 0; a < n_ab; a++)
            EVC_column(a) = EVC(a, i);

#ifdef NSA2DEBUG
        print_dmatrix2(reshape<dMatrix2>(EVC_column, Shape2D(n_ab, 1)), "slice used");
#endif

        s = eigvals[i];
        s2 = s * s;
        if (abs(s2 - 1.0) < 1E-8) c = 0.0;
        else c = sqrt(1.0 - s2);

        fm = sqrt(1 - c) / s;
        fp = sqrt(1 + c) / s;
        fa = 0.5 * ((fm + fp) + c * (fm - fp));
        fb = 0.5 * (c * (fm + fp) + (fm - fp));

        if (abs(fa - 1.0) > 1E-8) {
            A = dot_BLAS<dMatrix1, dMatrix2>(PAS, EVC_column, false);
            for (int a = 0; a < n_a; a++)
                A(a) /= fa;
        }
        if (abs(fa - 1.0) > 1E-8) {
            B = dot_BLAS<dMatrix1, dMatrix2>(PBS, EVC_column, false);
            for (int b = 0; b < n_b; b++)
                B(b) /= fb;
        }

#ifdef NSA2DEBUG
        std::cout << "fa: " << fa << std::endl << "fb: " << fb << std::endl;
        print_dmatrix2(reshape<dMatrix2>(A, Shape2D(n_a, 1)), "A");
        print_dmatrix2(reshape<dMatrix2>(B, Shape2D(n_b, 1)), "B");
#endif

        //build antibonding state in pairs[i]
        fa = 0.5 * (fm - fp);
        fb = 0.5 * (fm + fp);

        for (int a = 0; a < n_a; a++) {
            EVC(a, pairs[i]) = fa * A(a);
        }
        for (int b = n_a; b < n_ab; b++) {
            EVC(b, pairs[i]) = fb * B(b - n_a);
        }
    }
}

std::map<char, dMatrix2> Roby_information::make_covalent_from_ionic(
    const dMatrix2 &theta_I,
    const vec &eigvals,
    const ivec &pairs) {
    std::map<char, dMatrix2> res; // A = angle, V = eigen_value, T = Theta_vector
    const int size = eigvals.size();
    res.emplace('A', dMatrix2(size, 1));
    res.emplace('V', dMatrix2(size, 1));
    res.emplace('T', dMatrix2(theta_I.extent(0), theta_I.extent(1)));

    for (int i = 0; i < size; i++) {
        res['A'](i, 0) = 90.0;
        res['V'](i, 0) = 0.0;
    }

    for (int val = 0; val < size; val++) {
        if (pairs[val] < 0) continue;
        if (pairs[val] == val) continue;
        const double s = eigvals[val];
        if (s < eigvals[pairs[val]]) continue;

        for (int i = 0; i < theta_I.extent(0); i++) {
            res['T'](i, val) = constants::INV_SQRT2 * (theta_I(i, val) + theta_I(i, pairs[val]));
            res['T'](i, pairs[val]) = constants::INV_SQRT2 * (theta_I(i, val) - theta_I(i, pairs[val]));
        }

        const double c = sqrt(1.0 - s * s);
        res['V'](val, 0) = c;
        res['V'](pairs[val], 0) = -c;

        res['A'](val, 0) = atan2(s, c) * constants::INV_PI_180;
    }

    return res;
}

// Assembles the block-projected PAS matrix for a group: each atom's projection matrix
// is multiplied by the rectangular overlap block (atom bf rows × bond bf cols) and
// the results are stacked vertically.  Result shape: n_bf_G × n_ab.
dMatrix2 Roby_information::build_group_PAS(const ivec &group_atom_indices, const ivec &bond_bf_indices, const dMatrix2 &P_G) {
    const int n_ab = static_cast<int>(bond_bf_indices.size());
    // Collect group BF indices in the same sorted-atom order as P_G was built
    ivec group_bf;
    for (int ai : group_atom_indices)
        for (const auto &NAO : NAOs)
            if (NAO.atom_index == ai)
                for (int idx : NAO.matrix_elements)
                    group_bf.push_back(idx);
    const int n_G = static_cast<int>(group_bf.size());
    err_checkf(n_G > 0, "RGBI group PAS has no basis functions for the requested atom group.", std::cout);
    err_checkf(n_ab > 0, "RGBI group PAS has no pair basis functions.", std::cout);
    vec sub_S(n_G * n_ab, 0.0);
    get_submatrix(overlap_matrix, sub_S, group_bf, bond_bf_indices);
    dMatrix2 S_G = reshape<dMatrix2>(sub_S, Shape2D(n_G, n_ab));
    // PAS = P_G (n_G × n_G) × S_G (n_G × n_ab) → (n_G × n_ab)
    return dot<dMatrix2>(P_G, S_G, false, false);
}

// Group-aware version of transform_Ionic_eigenvectors_to_Ionic_orbitals.
// Builds PAS/PBS from block-assembled group projection matrices, then applies
// the same bonding/antibonding rotation as the single-atom version.
void Roby_information::transform_group_Ionic_orbitals(
    dMatrix2 &EVC,
    const vec &eigvals,
    const ivec &pairs,
    const ivec &group_a_atoms,
    const ivec &group_b_atoms,
    const ivec &bond_bf_indices,
    const dMatrix2 &P_GA,
    const dMatrix2 &P_GB)
{
    const int n_ab = EVC.extent(0);
    const int n_eigvals = static_cast<int>(eigvals.size());

    dMatrix2 PAS = build_group_PAS(group_a_atoms, bond_bf_indices, P_GA);
    dMatrix2 PBS = build_group_PAS(group_b_atoms, bond_bf_indices, P_GB);

    const int n_a = static_cast<int>(PAS.extent(0));
    const int n_b = static_cast<int>(PBS.extent(0));
    err_checkf(n_ab == n_a + n_b, "Inconsistent group sizes in transform_group_Ionic_orbitals", std::cout);

    dMatrix1 A(n_a), B(n_b);
    dMatrix1 EVC_column(n_ab);

    for (int i = 0; i < n_eigvals; i++) {
        if (pairs[i] < 0) continue;
        if (pairs[i] == i) continue;
        if (eigvals[i] < eigvals[pairs[i]]) continue;

        for (int a = 0; a < n_ab; a++)
            EVC_column(a) = EVC(a, i);

        const double s = eigvals[i];
        const double s2 = s * s;
        const double c = (abs(s2 - 1.0) < 1E-8) ? 0.0 : sqrt(1.0 - s2);

        const double fm = sqrt(1.0 - c) / s;
        const double fp = sqrt(1.0 + c) / s;
        double fa = 0.5 * ((fm + fp) + c * (fm - fp));
        double fb = 0.5 * (c * (fm + fp) + (fm - fp));

        if (abs(fa - 1.0) > 1E-8) {
            A = dot_BLAS<dMatrix1, dMatrix2>(PAS, EVC_column, false);
            for (int a = 0; a < n_a; a++)
                A(a) /= fa;
        }
        if (abs(fb - 1.0) > 1E-8) {
            B = dot_BLAS<dMatrix1, dMatrix2>(PBS, EVC_column, false);
            for (int b = 0; b < n_b; b++)
                B(b) /= fb;
        }

        // Build antibonding partner in column pairs[i]
        fa = 0.5 * (fm - fp);
        fb = 0.5 * (fm + fp);

        for (int a = 0; a < n_a; a++)
            EVC(a, pairs[i]) = fa * A(a);
        for (int b = n_a; b < n_ab; b++)
            EVC(b, pairs[i]) = fb * B(b - n_a);
    }
}

void Roby_information::computeGroupAnalysis(const ivec2 &group_defs, const vec &atom_pops, const ivec &atom_charges) {
    RGBI_groups.clear();
    const int n_groups = static_cast<int>(group_defs.size());

    // Helper: build NAO row/col index sets for a list of sorted atom indices
    auto collect_nao_indices = [&](const ivec &atoms, ivec &eigenvals_out, ivec &eigenvecs_out) {
        eigenvals_out.clear();
        eigenvecs_out.clear();
        int start_val = 0;
        for (const auto &NAO : NAOs) {
            bool wanted = std::find(atoms.begin(), atoms.end(), NAO.atom_index) != atoms.end();
            if (wanted) {
                for (int i = 0; i < static_cast<int>(NAO.eigenvalues.size()); i++)
                    eigenvals_out.push_back(start_val + i);
                for (int idx : NAO.matrix_elements)
                    eigenvecs_out.push_back(idx);
            }
            start_val += static_cast<int>(NAO.eigenvalues.size());
        }
        };

    // Helper: collect BF indices for a list of atom indices in the given order.
    // Caller must pass atoms in sorted order so BF ordering matches the ionic-operator layout.
    // Scans NAOs by atom_index rather than using direct array indexing, so the result is
    // correct regardless of the order atoms were processed in computeAllAtomicNAOs.
    auto collect_bf_indices = [&](const ivec &atoms, ivec &bf_out) {
        bf_out.clear();
        for (int ai : atoms)
            for (const auto &NAO : NAOs)
                if (NAO.atom_index == ai)
                    for (int idx : NAO.matrix_elements)
                        bf_out.push_back(idx);
        };

    auto atom_list_string = [](const ivec &atoms) {
        std::string result;
        for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
            if (i > 0)
                result += ",";
            result += std::to_string(atoms[i]);
        }
        return result;
        };

    // Pre-compute per-group data: sorted atoms, BF indices, NAO indices, projection, population
    struct GroupData {
        ivec sorted_atoms;
        ivec bf;
        ivec nao_evals;
        ivec nao_evecs;
        dMatrix2 P;
        double pop = 0.0;
        std::string elem_list; // e.g. "N,H,H,H"
    };
    std::vector<GroupData> gdata(n_groups);
    for (int gi = 0; gi < n_groups; gi++) {
        gdata[gi].sorted_atoms = group_defs[gi];
        std::sort(gdata[gi].sorted_atoms.begin(), gdata[gi].sorted_atoms.end());
        collect_bf_indices(gdata[gi].sorted_atoms, gdata[gi].bf);
        collect_nao_indices(gdata[gi].sorted_atoms, gdata[gi].nao_evals, gdata[gi].nao_evecs);
        err_checkf(!gdata[gi].bf.empty(),
            "RGBI group G" + std::to_string(gi) + " has no basis functions for atoms " +
            atom_list_string(gdata[gi].sorted_atoms) + ".", std::cout);
        gdata[gi].pop = projection_matrix_and_expectation(
            gdata[gi].bf, gdata[gi].nao_evals, gdata[gi].nao_evecs, nullptr, &gdata[gi].P);
        {
            std::map<int, int> charge_count;
            for (int ai : gdata[gi].sorted_atoms)
                charge_count[atom_charges[ai]]++;
            std::vector<std::pair<int, int>> by_charge(charge_count.begin(), charge_count.end());
            std::sort(by_charge.begin(), by_charge.end(), [](const auto &a, const auto &b) {
                return a.first > b.first; // heaviest element first
                });
            for (const auto &[charge, count] : by_charge) {
                gdata[gi].elem_list += constants::atnr2letter(charge);
                if (count > 1) gdata[gi].elem_list += std::to_string(count);
            }
        }
    }

    // ---- Header ----
    std::cout << "\n\nRoby-Gould Bond Indices (RGBI) - Group Analysis\n";
    std::cout << "----------------------------------------------\n";
    std::cout << "Groups defined:\n";
    for (int g = 0; g < n_groups; g++) {
        std::cout << "  G" << g << ": atoms ";
        for (int k = 0; k < static_cast<int>(gdata[g].sorted_atoms.size()); k++) {
            if (k > 0) std::cout << ", ";
            const int ai = gdata[g].sorted_atoms[k];
            std::cout << ai << "(" << constants::atnr2letter(atom_charges[ai]) << ")";
        }
        std::cout << "\n";
    }
    std::cout << "----------------------------------------------\n";

    // ---- Group populations ----
    std::cout << "\nGroup populations:\n";
    for (int g = 0; g < n_groups; g++) {
        const std::string label = "G" + std::to_string(g) + " (" + gdata[g].elem_list + ")";
        std::cout << "  Population of " << std::left << std::setw(24) << label
            << std::right << std::fixed << std::setprecision(5) << gdata[g].pop << "\n";
    }
    std::cout << "----------------------------------------------\n";

    // ---- Pair analysis ----
    for (int gi = 0; gi < n_groups; gi++) {
        for (int gj = gi + 1; gj < n_groups; gj++) {
            const ivec &ga_sorted = gdata[gi].sorted_atoms;
            const ivec &gb_sorted = gdata[gj].sorted_atoms;
            const ivec &ga_bf = gdata[gi].bf;
            const ivec &gb_bf = gdata[gj].bf;
            const dMatrix2 &P_GA = gdata[gi].P;
            const dMatrix2 &P_GB = gdata[gj].P;
            const double pop_ga = gdata[gi].pop;
            const double pop_gb = gdata[gj].pop;

            // bond_bf: ga's BFs first, gb's second — order MUST match ionic operator layout
            ivec bond_bf;
            bond_bf.insert(bond_bf.end(), ga_bf.begin(), ga_bf.end());
            bond_bf.insert(bond_bf.end(), gb_bf.begin(), gb_bf.end());

            // Pair population
            ivec bond_evals = gdata[gi].nao_evals, bond_evecs = gdata[gi].nao_evecs;
            bond_evals.insert(bond_evals.end(), gdata[gj].nao_evals.begin(), gdata[gj].nao_evals.end());
            bond_evecs.insert(bond_evecs.end(), gdata[gj].nao_evecs.begin(), gdata[gj].nao_evecs.end());
            err_checkf(!bond_bf.empty(),
                "RGBI group pair G" + std::to_string(gi) + "-G" + std::to_string(gj) +
                " has no basis functions.", std::cout);
            const bool pair_spans_full_basis =
                static_cast<int>(bond_bf.size()) == static_cast<int>(overlap_matrix.extent(0));
            // Full-basis group pairs use the same projection path as the total RGBI population.
            // This preserves the requested group BF ordering used by the established references.
            const double pair_pop = pair_spans_full_basis
                ? projection_matrix_and_expectation(bond_bf)
                : projection_matrix_and_expectation(bond_bf, bond_evals, bond_evecs);

            // Build ionic operator [P_GA | 0; 0 | -P_GB] using dense group projections
            const int n_a = static_cast<int>(ga_bf.size());
            const int n_b = static_cast<int>(gb_bf.size());
            const int n_ab = n_a + n_b;
            dMatrix2 Ionic_Operator(n_ab, n_ab);
            for (int r = 0; r < n_ab; r++)
                for (int c = 0; c < n_ab; c++)
                    Ionic_Operator(r, c) = 0.0;
            for (int r = 0; r < n_a; r++)
                for (int c = 0; c < n_a; c++)
                    Ionic_Operator(r, c) = P_GA(r, c);
            for (int r = 0; r < n_b; r++)
                for (int c = 0; c < n_b; c++)
                    Ionic_Operator(n_a + r, n_a + c) = -P_GB(r, c);

            // Sub-overlap → sqrt → solve eigenproblem
            vec S_Sub(n_ab * n_ab, 0.0);
            get_submatrix(overlap_matrix, S_Sub, bond_bf);
            vec V = S_Sub;
            vec W(n_ab);
            const vec Temp = mat_sqrt(V, W);
            dMatrix2 A = reshape<dMatrix2>(Temp, Shape2D(n_ab, n_ab));
            dMatrix2 SI = LAPACKE_invert(A);

            auto X = change_basis_general(Ionic_Operator, transpose(A), true);
            vec ionic_eigenvals(X.extent(0));
            make_Eigenvalues(X.container(), ionic_eigenvals);

            auto EVC = dot<dMatrix2>(SI, X);

            // Prune near-zero eigenvalues
            ivec non_zero;
            for (int i = 0; i < static_cast<int>(ionic_eigenvals.size()); i++)
                if (abs(ionic_eigenvals[i]) > 1E-5)
                    non_zero.push_back(i);

            vec pruned_eigvals;
            for (int idx : non_zero)
                pruned_eigvals.push_back(ionic_eigenvals[idx]);

            EVC = transpose(EVC);
            auto EVC2 = transpose(get_rectangle(EVC, non_zero));
            EVC.container().clear();

            const int n0 = static_cast<int>(pruned_eigvals.size());
            auto pairs = find_eigenvalue_pairs(pruned_eigvals);

            transform_group_Ionic_orbitals(EVC2, pruned_eigvals, pairs, ga_sorted, gb_sorted, bond_bf, P_GA, P_GB);

            auto covalent_info = make_covalent_from_ionic(EVC2, pruned_eigvals, pairs);

            // Per-orbital populations
            vec covalent_popul(n0), ionic_popul(n0);
            ivec vals;
            for (int i = 0; i < static_cast<int>(EVC2.extent(0)); i++)
                vals.emplace_back(i);
            for (int i = 0; i < n0; i++) {
                auto temp = transpose(covalent_info['T']);
                covalent_popul[i] = projection_matrix_and_expectation(bond_bf, { i }, vals, &temp);
                temp = transpose(EVC2);
                ionic_popul[i] = projection_matrix_and_expectation(bond_bf, { i }, vals, &temp);
            }

            const double zero_angle_cutoff = 1E-2 * constants::INV_PI_180;
            group_bond_index_result result;
            result.group_index_first = gi;
            result.group_index_second = gj;
            result.atoms_first = ga_sorted;
            result.atoms_second = gb_sorted;
            result.pair_population = pair_pop;
            result.population_first = pop_ga;
            result.population_second = pop_gb;
            result.covalent = 0.0;
            result.ionic = 0.0;

            // Lone-pair orbitals localized almost entirely within one group have
            // eigvals near ±1.  In the group basis the inter-group contamination
            // can push these to ~0.996 rather than exactly 1, so they evade the
            // angle cutoff (84.6° < 89.99°) and corrupt the ionic index with large
            // contributions of the wrong sign.  Exclude any pair whose positive
            // eigval exceeds this threshold.
            constexpr double lone_pair_eigval_threshold = 0.99;
            for (int i = 0; i < n0; i++) {
                if (covalent_info['A'](i, 0) < zero_angle_cutoff || covalent_info['A'](i, 0) > 90.0 - zero_angle_cutoff)
                    continue;
                if (pruned_eigvals[i] < pruned_eigvals[pairs[i]])
                    continue;
                if (pruned_eigvals[i] > lone_pair_eigval_threshold)
                    continue;
                if (pairs[i] != i) {
                    result.covalent += 0.5 * (covalent_popul[i] - covalent_popul[pairs[i]]);
                    result.ionic += 0.5 * (ionic_popul[i] - ionic_popul[pairs[i]]);
                }
            }

            const double b2 = result.covalent * result.covalent + result.ionic * result.ionic;
            result.percent_covalent_Pyth = (b2 > 1E-12) ? 100.0 * (result.covalent * result.covalent / b2) : 0.0;
            const double b = sqrt(b2);
            result.total = b;
            result.percent_covalent_Arakai = (b > 1E-12) ? 200.0 * abs(asin(result.covalent / b)) / constants::PI : 0.0;

            RGBI_groups.push_back(result);
        }
    }

    // ---- Results table ----
    const std::string sep(90, '-');
    std::cout << "\n" << std::left << std::setw(14) << "Group" << std::right
        << std::setw(8) << "n_G1"
        << std::setw(8) << "n_G2"
        << std::setw(8) << "n_G1G2"
        << std::setw(8) << "s_G1G2"
        << std::setw(8) << "Cov."
        << std::setw(8) << "Ion."
        << std::setw(8) << "Tot."
        << std::setw(8) << "Pyth."
        << std::setw(8) << "Arak."
        << "\n" << sep << "\n";
    for (const auto &res : RGBI_groups) {
        const std::string label = "G" + std::to_string(res.group_index_first)
            + " - G" + std::to_string(res.group_index_second);
        std::cout << std::left << std::setw(14) << label << std::right
            << std::fixed << std::setprecision(3)
            << std::setw(8) << res.population_first
            << std::setw(8) << res.population_second
            << std::setw(8) << res.pair_population
            << std::setw(8) << (res.population_first + res.population_second - res.pair_population)
            << std::setw(8) << res.covalent
            << std::setw(8) << res.ionic
            << std::setw(8) << res.total
            << std::setw(8) << res.percent_covalent_Pyth
            << std::setw(8) << res.percent_covalent_Arakai
            << "\n";
    }
    std::cout << sep << "\n";
}

Roby_information::Roby_information(WFN &wavy, const ivec3 &group_sets, const bool symmetrize, const bool use_ano_basis) {
    auto bonds = get_bonded_atom_pairs(wavy);
    const char *orbital_label = use_ano_basis ? "ANOs" : "NAOs";
    std::cout << "Calculating " << orbital_label << " for all atoms...                 " << std::flush;
    computeAllAtomicNAOs(wavy, symmetrize, use_ano_basis);
    std::cout << " ...done!" << std::endl;
    if (use_ano_basis) {
        if (ano_fallback_atoms.empty()) {
            std::cout << "ANO fallback summary: no atom-level fallbacks were needed." << std::endl;
        }
        else {
            std::cout << "ANO fallback summary: used molecular local-orbital fallback for "
                << ano_fallback_atoms.size() << "/" << wavy.get_ncen() << " atoms";
            std::cout << " (";
            for (int i = 0; i < static_cast<int>(ano_fallback_atoms.size()); ++i) {
                if (i > 0)
                    std::cout << ", ";
                std::cout << ano_fallback_atoms[i];
            }
            std::cout << ")." << std::endl;
        }
    }
    Shape2D NAOs_size;
    NAOs_size.cols = static_cast<int>(overlap_matrix.extent(0));
    for (size_t atom_idx = 0; atom_idx < NAOs.size(); atom_idx++) {
#ifdef NSA2DEBUG
        std::cout << "Atom " << atom_idx + 1 << " NAO Occupancies:\n";
#endif
        for (size_t i = 0; i < NAOs[atom_idx].eigenvalues.size(); i++) {
#ifdef NSA2DEBUG
            std::cout << "  NAO " << i + 1 << ": " << std::setprecision(6) << NAOs[atom_idx].eigenvalues[i] << "\n";
#endif
            NAOs_size.rows++;
        }
#ifdef NSA2DEBUG
        std::cout << std::endl;
#endif
    }
    //Create a matrix of size n_NAOs x n_bf and fill with subblocks of the evecs
    vec NAO_matrix(NAOs_size.cols * NAOs_size.rows, 0.0);
    total_NAOs = reshape<dMatrix2>(NAO_matrix, Shape2D(NAOs_size.rows, NAOs_size.cols));
    Shape2D temp = { 0,0 };
    for (auto NAO : NAOs) {
        const int ME = NAO.matrix_elements.size();
        const int NEV = NAO.eigenvalues.size();
        for (int i = 0; i < NEV; i++) {
            const int row = temp.rows + i;
            for (int j = 0; j < ME; j++) {
                const int index_eigenvector = i * ME + j;
                const int col = NAO.matrix_elements[j];
                total_NAOs(row, col) = NAO.eigenvectors[index_eigenvector];
            }
        }
        temp.rows += NAO.eigenvalues.size();
    }
#ifdef NSA2DEBUG
    print_dmatrix2(total_NAOs, "Global NAO Matrix");
#endif
    // total_NAOs has n_NAOs rows (each row is an individual NAO in the basis set)
    // and n_bfs columns, where each value corresponds to a basis function from here on
    std::cout << "Calculating Population for all atoms...           " << std::flush;
    const double all_atom_population = Roby_population_analysis({});
    std::cout << " ...done." << std::endl;
    std::cout << std::endl << "Total Population: " << all_atom_population << "\n\n";
    vec atom_pops(NAOs.size(), 0.0);
    projection_matrices.clear();
    projection_matrices.resize(NAOs.size());
    overlap_matrices.clear();
    overlap_matrices.resize(NAOs.size());
#ifndef NSA2DEBUG
    ProgressBar *pb = new ProgressBar(NAOs.size(), 40, "-", " ", "Calculating Atomic Populations");
#endif
    for (auto NAO : NAOs) {
        dMatrix2 P;
        atom_pops[NAO.atom_index] = projection_matrix_and_expectation(NAO.matrix_elements, {}, {}, nullptr, &P);
        projection_matrices[NAO.atom_index] = P;
#ifndef NSA2DEBUG
        pb->update();
#endif
    }
#ifndef NSA2DEBUG
    delete pb;
#endif
    for (int i = 0; i < atom_pops.size(); i++) {
        std::cout << "Population of atom " << i << ": " << atom_pops[i] << std::endl;
    }

#ifndef NSA2DEBUG
    pb = new ProgressBar(bonds.size(), 40, "-", " ", "Calculating Bond Populations");
#endif
    //now perform bond analysis for all bonded atoms
    for (auto bond : bonds) {
        //std::cout << std::endl << "---------------------------- Atom Pair: " << bond.first << " " << bond.second << " ----------------------\n";
        ivec bond_indices, bond_eigenvecs, bond_eigenvals;
        //gather basis function indices for both atoms
        for (const auto &NAO : NAOs) {
            if (NAO.atom_index == bond.first || NAO.atom_index == bond.second)
                for (auto idx : NAO.matrix_elements)
                    bond_indices.push_back(idx);
        }
        //just in case: sort bond indices, easy since each atom's indices are already assumed sorted and each basis function belongs to only one atom once
        std::sort(bond_indices.begin(), bond_indices.end());

        //now determine the bond atom NAO indices
        int start_val = 0;
        for (auto NAO : NAOs) {
            if (NAO.atom_index == bond.first || NAO.atom_index == bond.second) {
                const int n1 = NAO.eigenvalues.size();
                for (int i = 0; i < n1; i++) {
                    bond_eigenvals.push_back(start_val + i);
                }
                for (int idx : NAO.matrix_elements)
                    bond_eigenvecs.push_back(idx);
            }
            start_val += NAO.eigenvalues.size();
        }

        //calcualte population using data from both atoms
        const double bond_population = projection_matrix_and_expectation(bond_indices, bond_eigenvals, bond_eigenvecs);
        //atom_pair_populations(bond.first, bond.second) = bond_population;
        //atom_pair_populations(bond.second, bond.first) = bond_population;
        std::cout << "Bond population between atom " << bond.first + 1 << " and atom " << bond.second + 1 << ": " << bond_population << std::endl;
        const int size_ion1 = projection_matrices[bond.first].extent(0) + projection_matrices[bond.second].extent(0),
            size_ion2 = projection_matrices[bond.first].extent(1) + projection_matrices[bond.second].extent(1),
            size1_1 = projection_matrices[bond.first].extent(0),
            size1_2 = projection_matrices[bond.first].extent(1),
            size2_1 = projection_matrices[bond.second].extent(0),
            size2_2 = projection_matrices[bond.second].extent(1);
        dMatrix2 Ionic_Operator(size_ion1, size_ion2);
        for (int i = 0; i < size1_1; i++) {
            for (int j = 0; j < size1_2; j++) {
                Ionic_Operator(i, j) = projection_matrices[bond.first](i, j);
            }
        }
        for (int i = 0; i < size2_1; i++) {
            for (int j = 0; j < size2_2; j++) {
                Ionic_Operator(i + size1_1, j + size1_2) = -projection_matrices[bond.second](i, j);
            }
        }
#ifdef NSA2DEBUG
        print_dmatrix2(Ionic_Operator, "Ionic Operator");
#endif
        // get matching suboverlap matrix
        const int n = bond_indices.size();
        vec S_Sub(n * n, 0.0);
        get_submatrix(overlap_matrix, S_Sub, bond_indices);
        //dMatrix2 S = reshape<dMatrix2>(S_Sub, Shape2D(n, n));

        vec V = S_Sub;
        vec W(n);
        // make V = Sqrt(S)
        const vec Temp = mat_sqrt(V, W);

        dMatrix2 A = reshape<dMatrix2>(Temp, Shape2D(n, n));
#ifdef NSA2DEBUG
        print_dmatrix2(A, "Overlap Sqrt SH");
#endif

        dMatrix2 SI = LAPACKE_invert(A);

#ifdef NSA2DEBUG
        print_dmatrix2(SI, "Overlap Pseudo Inverse");
#endif

        auto X = change_basis_general(Ionic_Operator, transpose(A), true);
#ifdef NSA2DEBUG
        print_dmatrix2(X, "Overlap Eigenproblem");
#endif

        // solve symmetric eigenproblem of X
        vec ionic_eigenvals(X.extent(0));
        make_Eigenvalues(X.container(), ionic_eigenvals);

#ifdef NSA2DEBUG
        std::cout << "Ionic eigenvalues between atom " << bond.first + 1 << " and atom " << bond.second + 1 << ":\n";
        for (size_t i = 0; i < ionic_eigenvals.size(); i++) {
            std::cout << std::setw(5) << i + 1 << ": " << std::setw(10) << std::setprecision(6) << ionic_eigenvals[i] << "\n";
        }
        std::cout << std::endl;
#endif

        auto EVC = dot<dMatrix2>(SI, X);
#ifdef NSA2DEBUG
        print_dmatrix2(EVC, "theta_I");
#endif

        ivec non_zero_indices;
        for (int i = 0; i < ionic_eigenvals.size(); i++) {
            if (abs(ionic_eigenvals[i]) > 1E-5)
                non_zero_indices.push_back(i);
        }
        const int n0 = non_zero_indices.size();
#ifdef NSA2DEBUG
        std::cout << "Non zero eigenvalues:\n";
        for (size_t i = 0; i < n0; i++) {
            std::cout << std::setw(5) << i + 1 << ": " << std::setw(10) << std::setprecision(6) << non_zero_indices[i] << "\n";
        }
        std::cout << std::endl;
#endif

        vec pruned_eigvals;
        for (int nzv = 0; nzv < n0; nzv++)
            pruned_eigvals.push_back(ionic_eigenvals[non_zero_indices[nzv]]);
        EVC = transpose(EVC);
        auto EVC2 = transpose(get_rectangle(EVC, non_zero_indices));
#ifdef NSA2DEBUG
        print_dmatrix2(EVC2, "theta_I after pruning:");
#endif
        EVC.container().clear();

        auto pairs = find_eigenvalue_pairs(pruned_eigvals);
#ifdef NSA2DEBUG
        std::cout << "Pairs:\n";
        for (size_t i = 0; i < n0; i++) {
            std::cout << std::setw(3) << i << ": " << std::setw(5) << std::setprecision(6) << pairs[i] << "\n";
        }
        std::cout << std::endl;
#endif

        transform_Ionic_eigenvectors_to_Ionic_orbitals(EVC2, pruned_eigvals, pairs, bond.first, bond.second, bond_indices);
#ifdef NSA2DEBUG
        print_dmatrix2(EVC2, "theta_I after Ionic");
#endif
        //EVC2 is theta_Ionic, now make theta_Covalent
        auto covalent_info = make_covalent_from_ionic(EVC2, pruned_eigvals, pairs);
#ifdef NSA2DEBUG
        print_dmatrix2(covalent_info['A'], "theta_Angles");
        print_dmatrix2(covalent_info['V'], "eigen_C");
        print_dmatrix2(covalent_info['T'], "theta_C");
#endif
        vec covalent_popul(n0), ionic_popul(n0);
        ivec vals;
        for (int i = 0; i < EVC2.extent(0); i++)
            vals.emplace_back(i);
        for (int i = 0; i < n0; i++) {
            //make covalent populations
            auto temp = transpose(covalent_info['T']);
            covalent_popul[i] = projection_matrix_and_expectation(bond_indices, { i }, vals, &(temp));
            temp = transpose(EVC2);
            ionic_popul[i] = projection_matrix_and_expectation(bond_indices, { i }, vals, &(temp));
        }
#ifdef NSA2DEBUG
        std::cout << "Covalent populations:\n";
        for (int i = 0; i < n0; i++)
            std::cout << "\t" << i << ":\t" << covalent_popul[i] << std::endl;
        std::cout << "Ionic populations:\n";
        for (int i = 0; i < n0; i++)
            std::cout << "\t" << i << ":\t" << ionic_popul[i] << std::endl;
#endif

        const double zero_angle_cutoff = 1E-2 * constants::INV_PI_180;

        vec cov_index(n0, 0.0), ion_index(n0, 0.0);
        double b_ab = 0;
        bond_index_result results{};
        for (int i = 0; i < n0; i++) {
            if (covalent_info['A'](i, 0) < zero_angle_cutoff || covalent_info['A'](i, 0) > 90 - zero_angle_cutoff)
                continue; // skip lone pairs

            if (pruned_eigvals[i] < pruned_eigvals[pairs[i]])
                continue; // skip antibonding

            if (pairs[i] != i) {
                cov_index[i] = 0.5 * (covalent_popul[i] - covalent_popul[pairs[i]]);
                ion_index[i] = 0.5 * (ionic_popul[i] - ionic_popul[pairs[i]]);
                results.covalent += cov_index[i];
                results.ionic += ion_index[i];
            }
            else if (pruned_eigvals[i] > 0.0)
                ion_index[i] = 0.5 * ionic_popul[i];
            else if (pruned_eigvals[i] < 0.0)
                ion_index[i] = -0.5 * ionic_popul[i];
        }

        b_ab = results.covalent * results.covalent + results.ionic * results.ionic;
        results.percent_covalent_Pyth = 100 * (results.covalent * results.covalent / b_ab);
        b_ab = sqrt(b_ab);
        results.percent_covalent_Arakai = 200 * abs(asin(results.covalent / b_ab)) / constants::PI;

        int el_a = wavy.get_atom_charge(bond.first);
        int el_b = wavy.get_atom_charge(bond.second);
        if (el_a > el_b) {
            results.atom_indices = bond;
            results.atom_element_nr = std::make_pair(el_a, el_b);
            results.total = b_ab;
            results.pair_population = bond_population;
            results.population_first = atom_pops[bond.first];
            results.population_second = atom_pops[bond.second];
        }
        else {
            results.atom_indices = std::make_pair(bond.second, bond.first);
            results.atom_element_nr = std::make_pair(el_b, el_a);
            results.total = b_ab;
            results.pair_population = bond_population;
            results.population_first = atom_pops[bond.second];
            results.population_second = atom_pops[bond.first];
            results.ionic = -results.ionic;
        }
        RGBI.push_back(results);
#ifndef NSA2DEBUG
        pb->update();
#endif
    }
#ifndef NSA2DEBUG
    delete pb;
#endif
    const double number_of_electrons = wavy.get_nr_electrons();

    std::cout << "\nRoby-Gould Bond Indices (RGBI) Analysis\n----------------------------------------------\n";
    std::cout << "Number of electrons in system:         " << number_of_electrons << "\n";
    std::cout << "Number of electrons in Roby Analysis:  " << all_atom_population << "\n";
    std::cout << "Percentage of electrons accounted for: " << std::setprecision(4) << (100.0 * all_atom_population / number_of_electrons) << " %\n";
    std::cout << "----------------------------------------------\n";

    std::cout << "\n\nAtom Nr        Els  "
        << std::setw(8) << "n_A"
        << std::setw(8) << "n_B"
        << std::setw(8) << "n_AB"
        << std::setw(8) << "s_AB"
        << std::setw(8) << "Cov."
        << std::setw(8) << "Ion."
        << std::setw(8) << "Tot."
        << std::setw(8) << "Pyth."
        << std::setw(8) << "Arak."
        << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------\n";

    //Sort results by Element of the heavier atom of the bonds and then within each group by the heavier of the second atom
    std::sort(RGBI.begin(), RGBI.end(), [](const bond_index_result &a, const bond_index_result &b) {
        if (a.atom_element_nr.first != b.atom_element_nr.first) {
            return a.atom_element_nr.first > b.atom_element_nr.first;
        }
        return a.atom_element_nr.second > b.atom_element_nr.second;
        });

    for (auto res : RGBI) {
        std::cout << std::setw(4) << res.atom_indices.first << " -" << std::setw(4) << res.atom_indices.second << "  "
            << std::setw(3) << constants::atnr2letter(res.atom_element_nr.first) << " -" << std::setw(3) << constants::atnr2letter(res.atom_element_nr.second)
            << std::fixed << std::setprecision(3) << std::setw(8) << res.population_first
            << std::fixed << std::setprecision(3) << std::setw(8) << res.population_second
            << std::fixed << std::setprecision(3) << std::setw(8) << res.pair_population
            << std::fixed << std::setprecision(3) << std::setw(8) << res.population_first + res.population_second - res.pair_population
            << std::fixed << std::setprecision(3) << std::setw(8) << res.covalent
            << std::fixed << std::setprecision(3) << std::setw(8) << res.ionic
            << std::fixed << std::setprecision(3) << std::setw(8) << res.total
            << std::fixed << std::setprecision(3) << std::setw(8) << res.percent_covalent_Pyth
            << std::fixed << std::setprecision(3) << std::setw(8) << res.percent_covalent_Arakai << std::endl;
    }
    std::cout << "--------------------------------------------------------------------------------------------\n";

    if (!group_sets.empty()) {
        const int N_atoms = wavy.get_ncen();
        ivec atom_charges(N_atoms);
        for (int i = 0; i < N_atoms; i++)
            atom_charges[i] = wavy.get_atom_charge(i);
        for (const auto &group_defs : group_sets)
            computeGroupAnalysis(group_defs, atom_pops, atom_charges);
    }
}

void bondwise_laplacian_plots(std::filesystem::path &wfn_name)
{
    char cwd[1024];
#ifdef _WIN32
    if (_getcwd(cwd, sizeof(cwd)) != NULL)
#else
    if (getcwd(cwd, sizeof(cwd)) != NULL)
#endif
        std::cout << "Current working dir: " << cwd << std::endl;
    WFN wavy(wfn_name);
    wavy.delete_unoccupied_MOs();
    std::cout << NoSpherA2_message(false) << std::endl << build_date << std::endl;

    err_checkf(wavy.get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);

    int points = 1001;

    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        for (int j = i + 1; j < wavy.get_ncen(); j++)
        {
            std::filesystem::path path = cwd;
            double distance = sqrt(pow(wavy.get_atom_coordinate(i, 0) - wavy.get_atom_coordinate(j, 0), 2) + pow(wavy.get_atom_coordinate(i, 1) - wavy.get_atom_coordinate(j, 1), 2) + pow(wavy.get_atom_coordinate(i, 2) - wavy.get_atom_coordinate(j, 2), 2));
            double svdW = constants::ang2bohr(constants::covalent_radii[wavy.get_atom_charge(i)] + constants::covalent_radii[wavy.get_atom_charge(j)]);
            if (distance < 1.35 * svdW)
            {
                std::cout << "Bond between " << i << " (" << wavy.get_atom_charge(i) << ") and " << j << " (" << wavy.get_atom_charge(j) << ") with distance " << distance << " and svdW " << svdW << std::endl;
                const vec bond_vec = { (wavy.get_atom_coordinate(j, 0) - wavy.get_atom_coordinate(i, 0)) / points, (wavy.get_atom_coordinate(j, 1) - wavy.get_atom_coordinate(i, 1)) / points, (wavy.get_atom_coordinate(j, 2) - wavy.get_atom_coordinate(i, 2)) / points };
                const double dr = distance / points;
                vec lapl(points, 0.0);
                const vec pos = { wavy.get_atom_coordinate(i, 0), wavy.get_atom_coordinate(i, 1), wavy.get_atom_coordinate(i, 2) };
#pragma omp parallel for schedule(dynamic)
                for (int k = 0; k < points; k++)
                {
                    d3 t_pos = { pos[0], pos[1], pos[2] };
                    t_pos[0] += k * bond_vec[0];
                    t_pos[1] += k * bond_vec[1];
                    t_pos[2] += k * bond_vec[2];
                    lapl[k] = wavy.computeLap(t_pos);
                }
                std::filesystem::path outname(wfn_name.string() + "_bondwise_laplacian_" + std::to_string(i) + "_" + std::to_string(j) + ".dat");
                path = path / std::filesystem::path(outname);
                std::ofstream result(path, std::ios::out);
                for (int k = 0; k < points; k++)
                {
                    result << std::setw(10) << std::scientific << std::setprecision(6) << dr * k << " " << std::setw(10) << std::scientific << std::setprecision(6) << lapl[k] << "\n";
                }
                result.flush();
                result.close();
            }
            else
            {
                std::cout << "No bond between " << i << " and " << j << " with distance " << distance << " and svdW " << svdW << std::endl;
            }
        }
    }
}

void ELI_analysis(const WFN &wavy, const options &opt) {
    err_checkf(wavy.get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);
    std::cout << "Analysing ELI basins in the wavefunction..." << std::endl;

    const double radius = opt.properties.radius;
    const double grid_spacing = opt.properties.resolution;
    properties_options prop_opt = opt.properties;
    WFN l_w = wavy;
    l_w.delete_unoccupied_MOs();
    readxyzMinMax_fromWFN(wavy, prop_opt, true);

    cube rho(prop_opt.NbSteps, l_w.get_ncen(), true);
    cube eli_cube(prop_opt.NbSteps, l_w.get_ncen(), true);
    rho.give_parent_wfn(l_w);
    eli_cube.give_parent_wfn(l_w);

    std::cout << "Calcualting grid of size " << prop_opt.NbSteps[0] << " x " << prop_opt.NbSteps[1] << " x " << prop_opt.NbSteps[2] << "..." << std::endl;
    std::cout << "Number of points: " << prop_opt.NbSteps[0] * prop_opt.NbSteps[1] * prop_opt.NbSteps[2] << std::endl;

    //Print a summary of the grid parameters
    std::cout << "Grid parameters:\n";
    std::cout << "  Min: (" << prop_opt.MinMax[0] << ", " << prop_opt.MinMax[1] << ", " << prop_opt.MinMax[2] << ")\n";
    std::cout << "  Max: (" << prop_opt.MinMax[3] << ", " << prop_opt.MinMax[4] << ", " << prop_opt.MinMax[5] << ")\n";

    const std::vector<atom> atoms = wavy.get_atoms();

    // print table of atom positions
    std::cout << "Atom positions:\n";
    for (int a = 0; a < atoms.size(); a++) {
        std::cout << "  Atom " << a << ": (" << atoms[a].get_coordinate(0) << ", " << atoms[a].get_coordinate(1) << ", " << atoms[a].get_coordinate(2) << ")\n";
    }

    vec stepsizes{ (prop_opt.MinMax[3] - prop_opt.MinMax[0]) / prop_opt.NbSteps[0],
                  (prop_opt.MinMax[4] - prop_opt.MinMax[1]) / prop_opt.NbSteps[1],
                  (prop_opt.MinMax[5] - prop_opt.MinMax[2]) / prop_opt.NbSteps[2] };

    for (int i = 0; i < 3; i++)
    {
        rho.set_origin(i, prop_opt.MinMax[i]);
        rho.set_vector(i, i, stepsizes[i]);
        eli_cube.set_origin(i, prop_opt.MinMax[i]);
        eli_cube.set_vector(i, i, stepsizes[i]);
    }

    rho.calc_dv();
    eli_cube.calc_dv();

    Calc_RhoEli(rho, eli_cube, l_w, radius);
    rho.set_path("rho.cube");
    eli_cube.set_path("eli.cube");
    //rho.write_file(true);
    //eli_cube.write_file(true);

    const double density_floor = std::max(1e-8, rho.max_value() * 1e-6);
    const std::vector<critical_point> density_critical_points = analyze_cube_critical_points(&rho, l_w, opt.debug, density_floor);
    std::cout << "Density Critical Points";
    if (!density_critical_points.empty())
        std::cout << " (" << density_critical_points.size() << " found)";
    std::cout << ":\n";
    if (density_critical_points.empty()) {
        std::cout << "  No critical points refined from the current cube grid.\n";
    }
    else {
        constexpr int nw = 15; // numeric field width — wide enough for large core densities
        for (size_t i = 0; i < density_critical_points.size(); i++) {
            const critical_point &cp = density_critical_points[i];

            std::string ownership_info;
            if (cp.type == "attractor" || cp.type == "bond") {
                double min_dist1 = std::numeric_limits<double>::max();
                double min_dist2 = std::numeric_limits<double>::max();
                int atom_index1 = -1;
                int atom_index2 = -1;
                for (int a = 0; a < l_w.get_ncen(); a++) {
                    const d3 apos = l_w.get_atom_pos(a);
                    const double dx = cp.position[0] - apos[0];
                    const double dy = cp.position[1] - apos[1];
                    const double dz = cp.position[2] - apos[2];
                    const double dist2 = dx * dx + dy * dy + dz * dz;
                    if (dist2 < min_dist1) {
                        min_dist2 = min_dist1;
                        atom_index2 = atom_index1;
                        min_dist1 = dist2;
                        atom_index1 = a;
                    }
                    else if (dist2 < min_dist2) {
                        min_dist2 = dist2;
                        atom_index2 = a;
                    }
                }

                if (cp.type == "attractor" && atom_index1 >= 0) {
                    ownership_info = "   " + l_w.get_atom_label(atom_index1) + std::to_string(atom_index1);
                }
                else if (cp.type == "bond" && atom_index1 >= 0 && atom_index2 >= 0) {
                    ownership_info = "   "
                        + l_w.get_atom_label(atom_index1) + std::to_string(atom_index1)
                        + "-"
                        + l_w.get_atom_label(atom_index2) + std::to_string(atom_index2)
                        + " bond";
                }
            }

            std::cout << "\n  CP " << std::right << std::setw(3) << i + 1
                << "  [" << std::left << std::setw(12) << cp.type << "]  "
                << (cp.converged ? "converged" : "NOT converged")
                << ownership_info << "\n";
            std::cout << std::fixed << std::setprecision(4) << std::right;
            std::cout << "    Position  :"
                << std::setw(nw) << cp.position[0]
                << std::setw(nw) << cp.position[1]
                << std::setw(nw) << cp.position[2] << "\n";
            std::cout << "    Rho       :" << std::setw(nw) << cp.density << "\n";
            std::cout << "    GradRho   :"
                << std::setw(nw) << cp.gradient[0]
                << std::setw(nw) << cp.gradient[1]
                << std::setw(nw) << cp.gradient[2]
                << "   |GradRho| :" << std::setw(nw) << cp.gradient_norm << "\n";
            std::cout << "    HessRho_EigVals:"
                << std::scientific << std::setprecision(4)
                << std::setw(nw) << cp.hessian_eigenvalues[0]
                << std::setw(nw) << cp.hessian_eigenvalues[1]
                << std::setw(nw) << cp.hessian_eigenvalues[2]
                << std::fixed << std::setprecision(4) << "\n";
            std::cout << "    HessRho_EigVecs v1:"
                << std::setw(nw) << cp.hessian_eigenvectors[0][0]
                << std::setw(nw) << cp.hessian_eigenvectors[0][1]
                << std::setw(nw) << cp.hessian_eigenvectors[0][2] << "\n";
            std::cout << "                    v2:"
                << std::setw(nw) << cp.hessian_eigenvectors[1][0]
                << std::setw(nw) << cp.hessian_eigenvectors[1][1]
                << std::setw(nw) << cp.hessian_eigenvectors[1][2] << "\n";
            std::cout << "                    v3:"
                << std::setw(nw) << cp.hessian_eigenvectors[2][0]
                << std::setw(nw) << cp.hessian_eigenvectors[2][1]
                << std::setw(nw) << cp.hessian_eigenvectors[2][2] << "\n";
            std::cout << "    DelSqRho  :" << std::setw(nw) << cp.laplacian << "\n";
            if (std::isfinite(cp.ellipticity))
                std::cout << "    Bond Ellipticity:" << std::setw(nw) << cp.ellipticity << "\n";
            if (std::isfinite(cp.virial_field) || std::isfinite(cp.kinetic_lagrangian) || std::isfinite(cp.kinetic_hamiltonian) || std::isfinite(cp.lagrangian_density)) {
                std::cout << "   V         :" << std::scientific << std::setprecision(4) << std::setw(nw) << cp.virial_field
                    << "   G         :" << std::fixed << std::setprecision(4) << std::setw(nw) << cp.kinetic_lagrangian
                    << "   K         :" << std::scientific << std::setprecision(4) << std::setw(nw) << cp.kinetic_hamiltonian
                    << "   L         :" << std::scientific << std::setprecision(4) << std::setw(nw) << cp.lagrangian_density
                    << std::fixed << std::setprecision(4) << "\n";
            }
        }
        std::cout << "\n";
    }

    std::pair<cubei, std::vector<d4>> qtaim_results = topological_cube_analysis(&rho, atoms, opt.debug, true, 0.0, 1e-10, radius);
    svec labels = assign_labels_to_basins(qtaim_results.second, atoms, opt.debug);

    std::pair<cubei, std::vector<d4>> eli_results = topological_cube_analysis(&eli_cube, atoms, opt.debug, false, 0.0, 1e-10, radius);
    svec eli_labels = assign_labels_to_basins(eli_results.second, atoms, opt.debug, 1);

    std::cout << "QTAIM Analysis:\n";
    integrate_values_in_basins(&rho, &(qtaim_results.first), labels, opt.debug);
    std::cout << "\n\nELI Analysis:\n";
    integrate_values_in_basins(&rho, &(eli_results.first), eli_labels, opt.debug);

}

// ---------------------------------------------------------------------------
// QTAIM_ELI_mask
// ---------------------------------------------------------------------------

void QTAIM_ELI_mask(
    cube& rho,
    cube& eli,
    WFN& parent_wfn,
    const std::vector<atom>& atoms,
    const std::vector<int>& selected_indices,
    double background_value,
    const std::filesystem::path& output_path,
    bool debug,
    std::ostream& log
) {
    err_checkf(!atoms.empty(), "QTAIM_ELI_mask: no atoms available.", log);
    err_checkf(rho.get_size(0) == eli.get_size(0) &&
               rho.get_size(1) == eli.get_size(1) &&
               rho.get_size(2) == eli.get_size(2),
               "QTAIM_ELI_mask: rho and eli grids have different dimensions.", log);

    // Validate atom indices
    for (int idx : selected_indices) {
        err_checkf(idx >= 0 && idx < (int)atoms.size(),
            "QTAIM_ELI_mask: atom index " + std::to_string(idx) +
            " is out of range (0.." + std::to_string((int)atoms.size() - 1) + ").", log);
    }

    log << "Running QTAIM topological analysis on density grid..." << std::endl;
    const double density_floor = std::max(1e-8, rho.max_value() * 1e-6);
    auto [basin_cube, maxima] = topological_cube_analysis(&rho, atoms, debug, false, density_floor);
    svec labels = assign_labels_to_basins(maxima, atoms, debug, 0);

    // Map selected atom indices → set of 1-based basin IDs
    // Label format: element + atom_index, e.g. "C0", "H3"
    std::set<int> selected_basins;
    for (int b = 0; b < (int)labels.size(); b++) {
        const std::string& lbl = labels[b];
        // Extract trailing digit sequence (the atom index)
        auto it = std::find_if(lbl.rbegin(), lbl.rend(),
                               [](char c) { return !std::isdigit(static_cast<unsigned char>(c)); });
        if (it == lbl.rend()) continue; // entire string is digits — skip
        const std::string num_str(it.base(), lbl.end());
        if (num_str.empty()) continue;
        int atom_idx = std::stoi(num_str);
        if (std::find(selected_indices.begin(), selected_indices.end(), atom_idx)
                != selected_indices.end()) {
            selected_basins.insert(b + 1); // basin IDs in cubei are 1-indexed
            if (debug)
                log << "  Basin " << (b + 1) << " (\"" << lbl << "\") selected.\n";
        }
    }

    if (selected_basins.empty()) {
        log << "Warning: no QTAIM basins matched the requested atom indices. "
               "Output will contain only the background value.\n";
    }

    const int nx = eli.get_size(0);
    const int ny = eli.get_size(1);
    const int nz = eli.get_size(2);

    // Pass 1: find tight bounding box of selected voxels
    int xmin = nx, xmax = -1, ymin = ny, ymax = -1, zmin = nz, zmax = -1;
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
            for (int z = 0; z < nz; z++) {
                if (selected_basins.count(basin_cube.get_value(x, y, z))) {
                    xmin = std::min(xmin, x); xmax = std::max(xmax, x);
                    ymin = std::min(ymin, y); ymax = std::max(ymax, y);
                    zmin = std::min(zmin, z); zmax = std::max(zmax, z);
                }
            }

    if (xmax < 0) {
        // No selected voxels — write a 1×1×1 cube with background value
        xmin = xmax = 0; ymin = ymax = 0; zmin = zmax = 0;
    }

    log << "Bounding box of selected basins: ["
        << xmin << ".." << xmax << "] x ["
        << ymin << ".." << ymax << "] x ["
        << zmin << ".." << zmax << "]\n";

    // Pass 2: build shrunk cube
    const std::array<int, 3> new_size = {xmax - xmin + 1, ymax - ymin + 1, zmax - zmin + 1};
    cube shrunk(new_size, (int)atoms.size(), true);
    shrunk.give_parent_wfn(parent_wfn);

    const auto new_origin = eli.get_pos(xmin, ymin, zmin);
    for (int i = 0; i < 3; i++) {
        shrunk.set_origin(i, new_origin[i]);
        for (int j = 0; j < 3; j++)
            shrunk.set_vector(i, j, eli.get_vector(i, j));
    }
    shrunk.calc_dv();

    for (int x = xmin; x <= xmax; x++)
        for (int y = ymin; y <= ymax; y++)
            for (int z = zmin; z <= zmax; z++) {
                double val = selected_basins.count(basin_cube.get_value(x, y, z))
                             ? eli.get_value(x, y, z)
                             : background_value;
                shrunk.set_value(x - xmin, y - ymin, z - zmin, val);
            }

    shrunk.set_comment1("QTAIM-masked ELI");
    shrunk.set_comment2("Selected atoms: " + [&]() {
        std::string s;
        for (int i = 0; i < (int)selected_indices.size(); i++) {
            if (i) s += ',';
            s += std::to_string(selected_indices[i]);
        }
        return s;
    }());
    shrunk.set_path(output_path);

    log << "Writing masked ELI cube to " << output_path.string() << " ..." << std::endl;
    err_checkf(shrunk.write_file(true), "QTAIM_ELI_mask: failed to write output cube.", log);
    log << "Done.\n";
}

// ---------------------------------------------------------------------------
// run_QTAIM_ELI_mask  (dispatch: cube-files mode vs WFN mode)
// ---------------------------------------------------------------------------

void run_QTAIM_ELI_mask(
    const std::filesystem::path& rho_or_wfn,
    const std::filesystem::path& eli_path,
    const std::vector<int>& selected_indices,
    double background_value,
    const options& opt,
    std::ostream& log
) {
    if (!eli_path.empty()) {
        // ---- Cube-files mode ----
        log << "Reading density cube: " << rho_or_wfn.string() << std::endl;
        WFN dummy;
        cube rho(rho_or_wfn, true, dummy, log);
        log << "Reading ELI cube:     " << eli_path.string() << std::endl;
        cube eli(eli_path, true, dummy, log);
        rho.give_parent_wfn(dummy);
        eli.give_parent_wfn(dummy);

        const std::vector<atom> atoms = rho.get_parent_wfn_atoms();
        const std::filesystem::path out =
            eli_path.parent_path() / (eli_path.stem().string() + "_qtaim_masked.cube");

        QTAIM_ELI_mask(rho, eli, dummy, atoms, selected_indices, background_value, out, opt.debug, log);
    } else {
        // ---- WFN mode: compute rho and eli cubes first ----
        log << "Loading wavefunction: " << rho_or_wfn.string() << std::endl;
        WFN wavy(rho_or_wfn, opt.debug);
        wavy.delete_unoccupied_MOs();

        properties_options prop_opt = opt.properties;
        readxyzMinMax_fromWFN(wavy, prop_opt, true);

        log << "Calculating density and ELI grid ("
            << prop_opt.NbSteps[0] << " x " << prop_opt.NbSteps[1] << " x " << prop_opt.NbSteps[2]
            << ") ..." << std::endl;

        cube rho(prop_opt.NbSteps, wavy.get_ncen(), true);
        cube eli(prop_opt.NbSteps, wavy.get_ncen(), true);
        rho.give_parent_wfn(wavy);
        eli.give_parent_wfn(wavy);

        const vec stepsizes{
            (prop_opt.MinMax[3] - prop_opt.MinMax[0]) / prop_opt.NbSteps[0],
            (prop_opt.MinMax[4] - prop_opt.MinMax[1]) / prop_opt.NbSteps[1],
            (prop_opt.MinMax[5] - prop_opt.MinMax[2]) / prop_opt.NbSteps[2]
        };
        for (int i = 0; i < 3; i++) {
            rho.set_origin(i, prop_opt.MinMax[i]);
            rho.set_vector(i, i, stepsizes[i]);
            eli.set_origin(i, prop_opt.MinMax[i]);
            eli.set_vector(i, i, stepsizes[i]);
        }
        rho.calc_dv();
        eli.calc_dv();

        Calc_RhoEli(rho, eli, wavy, prop_opt.radius);

        const std::vector<atom> atoms = wavy.get_atoms();
        const std::filesystem::path out =
            rho_or_wfn.parent_path() / "eli_qtaim_masked.cube";

        QTAIM_ELI_mask(rho, eli, wavy, atoms, selected_indices, background_value, out, opt.debug, log);
    }
}
