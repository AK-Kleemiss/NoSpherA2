#include "pch.h"
#include "wfn_class.h"
#include "properties.h"
#include "convenience.h"
#include "spherical_density.h"
#include "cube.h"
#include "constants.h"

void print_time(_time_point &start, _time_point &end, std::ostream &file) {
    if (get_sec(start, end) < 60)
        file << "Time to calculate Values: " << std::fixed << std::setprecision(0) << get_sec(start, end) << " s" << std::endl;
    else if (get_sec(start, end) < 3600)
        file << "Time to calculate Values: " << std::fixed << std::setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << std::endl;
    else
        file << "Time to calculate Values: " << std::fixed << std::setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << std::endl;
}

namespace {

bool is_within_radius(const d3 &pos, const std::vector<atom> &atoms, double radius_bohr)
{
    for (const atom &entry : atoms)
        if (array_length(pos, entry.get_pos()) < radius_bohr)
            return true;
    return false;
}

double sanitize_finite(double value)
{
    if (std::isnan(value) || std::isinf(value))
        return 0.0;
    return value;
}

template <typename EvalFn>
void evaluate_cube_in_radius(
    cube &target,
    bool wrap,
    const std::vector<atom> &atoms,
    double radius_bohr,
    EvalFn &&evaluate_inside)
{
    target.evaluate_on_grid(
        [&](const d3 &pos) {
            if (!is_within_radius(pos, atoms, radius_bohr))
                return 0.0;
            return evaluate_inside(pos);
        },
        wrap);
}

template <typename EvalFn>
void evaluate_cube_in_radius_mapped(
    cube &target,
    bool wrap,
    const std::vector<atom> &atoms,
    double radius_bohr,
    EvalFn &&evaluate_inside_mapped)
{
    target.evaluate_on_grid(
        [&](const d3 &pos, const i3 &raw_idx, const i3 &mapped_idx) {
            if (!is_within_radius(pos, atoms, radius_bohr))
                return 0.0;
            return evaluate_inside_mapped(pos, raw_idx, mapped_idx);
        },
        wrap);
}

template <typename EvalFn>
void evaluate_cube_near_atom_mapped(
    cube &target,
    bool wrap,
    const d3 &atom_pos,
    double radius_bohr,
    EvalFn &&evaluate_inside_mapped)
{
    target.evaluate_on_grid(
        [&](const d3 &pos, const i3 &raw_idx, const i3 &mapped_idx) {
            if (array_length(pos, atom_pos) >= radius_bohr)
                return 0.0;
            return evaluate_inside_mapped(pos, raw_idx, mapped_idx);
        },
        wrap);
}

struct PropValues {
    double rho = 0.0;
    double grad = 0.0;
    double elf = 0.0;
    double eli = 0.0;
    double lap = 0.0;
    double hess[9]{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
};

PropValues compute_prop_values(const std::vector<cube> &cubes, const WFN &wavy, const d3 &pos)
{
    PropValues values;

    const bool rdg_loaded = cubes[cube_type::RDG].get_loaded();
    const bool lap_loaded = cubes[cube_type::Lap].get_loaded();
    const bool elf_loaded = cubes[cube_type::Elf].get_loaded();
    const bool eli_loaded = cubes[cube_type::Eli].get_loaded();

    if (cubes[cube_type::ESP].get_loaded() && !rdg_loaded)
        values.rho = wavy.compute_dens(pos);

    if (rdg_loaded && lap_loaded && (elf_loaded || eli_loaded))
        wavy.computeValues(pos, values.rho, values.grad, values.hess, values.elf, values.eli, values.lap);
    else if (elf_loaded && eli_loaded && !rdg_loaded && !lap_loaded)
        wavy.computeELIELF(pos, values.elf, values.eli);
    else if (elf_loaded && !eli_loaded && !rdg_loaded && !lap_loaded)
        values.elf = wavy.computeELF(pos);
    else if (!elf_loaded && eli_loaded && !rdg_loaded && !lap_loaded)
        values.eli = wavy.computeELI(pos);
    else if (elf_loaded && eli_loaded && lap_loaded && !rdg_loaded)
        wavy.computeLapELIELF(pos, values.elf, values.eli, values.lap);
    else if (!elf_loaded && eli_loaded && lap_loaded && !rdg_loaded)
        wavy.computeLapELI(pos, values.eli, values.lap);
    else if (!elf_loaded && !eli_loaded && lap_loaded && !rdg_loaded)
        values.lap = wavy.computeLap(pos);
    else
        wavy.computeValues(pos, values.rho, values.grad, values.hess, values.elf, values.eli, values.lap);

    if (rdg_loaded)
        values.rho = get_lambda_1(values.hess) < 0 ? -values.rho : values.rho;

    return values;
}

void accumulate_prop_values(std::vector<cube> &cubes, const i3 &mapped_idx, const PropValues &values)
{
    const int x = mapped_idx[0];
    const int y = mapped_idx[1];
    const int z = mapped_idx[2];

    if (cubes[cube_type::RDG].get_loaded())
        cubes[cube_type::RDG].set_value(x, y, z, cubes[cube_type::RDG].get_value(x, y, z) + sanitize_finite(values.grad));
    if (cubes[cube_type::Lap].get_loaded())
        cubes[cube_type::Lap].set_value(x, y, z, cubes[cube_type::Lap].get_value(x, y, z) + sanitize_finite(values.lap));
    if (cubes[cube_type::Elf].get_loaded())
        cubes[cube_type::Elf].set_value(x, y, z, cubes[cube_type::Elf].get_value(x, y, z) + sanitize_finite(values.elf));
    if (cubes[cube_type::Eli].get_loaded())
        cubes[cube_type::Eli].set_value(x, y, z, cubes[cube_type::Eli].get_value(x, y, z) + sanitize_finite(values.eli));
}

} // namespace

void Calc_Spherical_Dens(
    cube &CubeSpher,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();

    vector<Thakkar> atom_models;
    for (int a = 0; a < 92; a++) {
        atom_models.emplace_back(a);
        atom_models[a].make_interpolator(1.005 * 1.005 * 1.005, 1E-7);
    }
    const double radius_bohr = constants::ang2bohr(radius);
    const vector<atom> wavy_atoms = wavy.get_atoms();

    evaluate_cube_in_radius(
        CubeSpher,
        wrap,
        wavy_atoms,
        radius_bohr,
        [&](const d3 &pos) {
            vector<double> dists(wavy.get_ncen(), 0.0);
            for (int a = 0; a < wavy.get_ncen(); a++)
                dists[a] = array_length(pos, wavy.get_atom_pos(a));

            double dens_all = 0.0;
            for (int a = 0; a < wavy.get_ncen(); a++)
                dens_all += atom_models[wavy.get_atom_charge(a) - 1].get_interpolated_density(dists[a]);
            return dens_all;
        });

    _time_point end = get_time();
    print_time(start, end, file);
};

void Calc_Static_Def(
    cube &CubeDEF,
    cube &CubeRho,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    vector<Thakkar> atoms;
    atoms.reserve(wavy.get_ncen());
    for (int a = 0; a < wavy.get_ncen(); a++)
        atoms.emplace_back(wavy.get_atom_charge(a));
    const double radius_bohr = constants::ang2bohr(radius);
    const vector<atom> wavy_atoms = wavy.get_atoms();

    evaluate_cube_in_radius_mapped(
        CubeDEF,
        wrap,
        wavy_atoms,
        radius_bohr,
        [&](const d3 &pos, const i3 &, const i3 &mapped_idx) {
            vector<double> dists(wavy.get_ncen(), 0.0);
            for (int a = 0; a < wavy.get_ncen(); a++)
                dists[a] = array_length(pos, wavy.get_atom_pos(a));

            double dens_all = 0.0;
            for (int a = 0; a < wavy.get_ncen(); a++)
                dens_all += atoms[a].get_radial_density(dists[a]);

            dens_all -= CubeRho.get_value(mapped_idx[0], mapped_idx[1], mapped_idx[2]);
            return -dens_all;
        });

    _time_point end = get_time();
    print_time(start, end, file);
};

void Calc_Static_Def(
    std::vector<cube> &Cubes,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    const double radius_bohr = constants::ang2bohr(radius);
    const vector<atom> wavy_atoms = wavy.get_atoms();

    evaluate_cube_in_radius_mapped(
        Cubes[cube_type::DEF],
        wrap,
        wavy_atoms,
        radius_bohr,
        [&](const d3 &, const i3 &, const i3 &mapped_idx) {
            const int x = mapped_idx[0];
            const int y = mapped_idx[1];
            const int z = mapped_idx[2];
            const double rho = Cubes[cube_type::Rho].get_value(x, y, z);
            const double spher = Cubes[cube_type::spherical_density].get_value(x, y, z);
            return rho - spher;
        });

    _time_point end = get_time();
    print_time(start, end, file);
};

void Calc_Hirshfeld(
    std::vector<cube> &Cubes,
    const WFN &wavy,
    double radius,
    int ignore_atom,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();

    vector<Thakkar> atoms;
    atoms.reserve(wavy.get_ncen());
    for (int a = 0; a < wavy.get_ncen(); a++)
        atoms.emplace_back(wavy.get_atom_charge(a));
    const double radius_bohr = constants::ang2bohr(radius);
    const vector<atom> focus_atom{ wavy.get_atoms()[ignore_atom] };

    evaluate_cube_in_radius_mapped(
        Cubes[cube_type::HDEF],
        wrap,
        focus_atom,
        radius_bohr,
        [&](const d3 &pos, const i3 &, const i3 &mapped_idx) {
            double dens_choice = 0.0;
            double dens_all = 0.0;
            for (int a = 0; a < wavy.get_ncen(); a++)
            {
                const double dist = array_length(pos, wavy.get_atom_pos(a));
                const double temp = atoms[a].get_radial_density(dist);
                if (ignore_atom == a)
                    dens_choice = temp;
                dens_all += temp;
            }

            if (std::abs(dens_all) < 1E-20)
                return 0.0;

            const double rho = Cubes[cube_type::Rho].get_value(mapped_idx[0], mapped_idx[1], mapped_idx[2]);
            return dens_choice / dens_all * rho - dens_choice;
        });

    _time_point end = get_time();
    print_time(start, end, file);
};

void Calc_Hirshfeld(
    cube &CubeHDEF,
    cube &CubeRho,
    cube &CubeSpherical,
    const WFN &wavy,
    double radius,
    int ignore_atom,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    Thakkar atom_model(wavy.get_atom_charge(ignore_atom));
    const double radius_bohr = constants::ang2bohr(radius);

    evaluate_cube_near_atom_mapped(
        CubeHDEF,
        wrap,
        wavy.get_atom_pos(ignore_atom),
        radius_bohr,
        [&](const d3 &pos, const i3 &, const i3 &mapped_idx) {
            const double dist = array_length(pos, wavy.get_atom_pos(ignore_atom));
            const int x = mapped_idx[0];
            const int y = mapped_idx[1];
            const int z = mapped_idx[2];
            const double spherical = CubeSpherical.get_value(x, y, z);
            if (std::abs(spherical) < 1E-20)
                return 0.0;

            const double dens_choice = atom_model.get_radial_density(dist);
            return (dens_choice / spherical * CubeRho.get_value(x, y, z)) - dens_choice;
        });

    _time_point end = get_time();
    print_time(start, end, file);
};

void Calc_Hirshfeld_atom(
    std::vector<cube> &Cubes,
    const WFN &wavy,
    double radius,
    int ignore_atom,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    Thakkar atom_model(wavy.get_atom_charge(ignore_atom));
    const double radius_bohr = constants::ang2bohr(radius);

    evaluate_cube_near_atom_mapped(
        Cubes[cube_type::Hirsh],
        wrap,
        wavy.get_atom_pos(ignore_atom),
        radius_bohr,
        [&](const d3 &pos, const i3 &, const i3 &mapped_idx) {
            const double dist = array_length(pos, wavy.get_atom_pos(ignore_atom));
            const double dens_choice = atom_model.get_radial_density(dist);
            const int x = mapped_idx[0];
            const int y = mapped_idx[1];
            const int z = mapped_idx[2];
            const double spherical = Cubes[cube_type::spherical_density].get_value(x, y, z);
            if (std::abs(spherical) < 1E-20)
                return 0.0;

            return dens_choice / spherical * Cubes[cube_type::Rho].get_value(x, y, z);
        });

    _time_point end = get_time();
    print_time(start, end, file);
};

void Calc_Rho(
    cube &CubeRho,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    const double radius_bohr = constants::ang2bohr(radius);
    const vector<atom> atoms = wavy.get_atoms();

    evaluate_cube_in_radius(
        CubeRho,
        wrap,
        atoms,
        radius_bohr,
        [&](const d3 &pos) {
            return wavy.compute_dens(pos);
        });

    _time_point end = get_time();
    print_time(start, end, file);
};

void Calc_Eli(
    cube &CubeEli,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    const double radius_bohr = constants::ang2bohr(radius);
    const vector<atom> atoms = wavy.get_atoms();

    evaluate_cube_in_radius(
        CubeEli,
        wrap,
        atoms,
        radius_bohr,
        [&](const d3 &pos) {
            return wavy.computeELI(pos);
        });

    _time_point end = get_time();
    print_time(start, end, file);
};

void Calc_RhoEli(
    cube &CubeRho,
    cube &CubeEli,
    const WFN &wavy,
    double radius)
{
    using namespace std;
    _time_point start = get_time();
    err_checkf(CubeRho.get_size(0) == CubeEli.get_size(0) && CubeRho.get_size(1) == CubeEli.get_size(1) && CubeRho.get_size(2) == CubeEli.get_size(2), "Cube sizes do not match", std::cout);

    const double radius_bohr = constants::ang2bohr(radius);
    const vector<atom> atoms = wavy.get_atoms();
    CubeEli.set_zero();

    evaluate_cube_in_radius_mapped(
        CubeRho,
        false,
        atoms,
        radius_bohr,
        [&](const d3 &pos, const i3 &, const i3 &mapped_idx) {
            double rho = 0.0;
            double eli = 0.0;
            wavy.computeRhoELI(pos, rho, eli);
            CubeEli.set_value(mapped_idx[0], mapped_idx[1], mapped_idx[2], eli);
            return rho;
        });

    _time_point end = get_time();
    print_time(start, end, std::cout);
};

void Calc_Rho_spherical_harmonics(
    cube &CubeRho,
    WFN &wavy,
    std::ostream &file)
{
    using namespace std;
    _time_point start = get_time();

    ProgressBar *progress = new ProgressBar(CubeRho.get_size(0), 50, "=", " ", "Calculating Rho");

#pragma omp parallel shared(CubeRho)
    {
        vec2 d(wavy.get_ncen());
        for (int i = 0; i < wavy.get_ncen(); i++)
            d[i].resize(16, 0.0);
        const int n = wavy.get_nmo(true);
        vec phi(n, 0.0);
        // #pragma omp for schedule(dynamic)
        for (int i = 0; i < CubeRho.get_size(0); i++)
        {
            for (int j = 0; j < CubeRho.get_size(1); j++)
                for (int k = 0; k < CubeRho.get_size(2); k++)
                    CubeRho.set_value(i, j, k, wavy.compute_dens(CubeRho.get_pos(i, j, k), d, phi));
            progress->update();
        }
    }
    delete (progress);

    _time_point end = get_time();
    print_time(start, end, file);
};

void Calc_MO_spherical_harmonics(
    cube &CubeMO,
    WFN &wavy,
    int MO,
    std::ostream &file,
    bool nodate)
{
    using namespace std;
    _time_point start = get_time();
    ProgressBar *progress = NULL;
    if (!nodate)
        progress = new ProgressBar(CubeMO.get_size(0), 50, "=", " ", "Calculating Values");

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < CubeMO.get_size(0); i++)
    {
        for (int j = 0; j < CubeMO.get_size(1); j++)
            for (int k = 0; k < CubeMO.get_size(2); k++)
                CubeMO.set_value(i, j, k, wavy.compute_MO_spherical(CubeMO.get_pos(i, j, k), MO));
        if (!nodate)
            progress->update();
    }
    if (!nodate)
    {
        delete (progress);

        _time_point end = get_time();
        print_time(start, end, file);
    }
};

void Calc_S_Rho(
    cube &Cube_S_Rho,
    WFN &wavy,
    std::ostream &file,
    bool &nodate)
{
    using namespace std;
    _time_point start = get_time();
    ProgressBar *progress = NULL;
    if (!nodate)
        progress = new ProgressBar(Cube_S_Rho.get_size(0), 50, "=", " ", "Calculating Values");

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < Cube_S_Rho.get_size(0); i++)
    {
        vec2 d(16);
        vec phi(wavy.get_nmo(), 0.0);
        for (int p = 0; p < 16; p++)
            d[p].resize(wavy.get_ncen());
        for (int j = 0; j < Cube_S_Rho.get_size(1); j++)
            for (int k = 0; k < Cube_S_Rho.get_size(2); k++)
                Cube_S_Rho.set_value(i, j, k, wavy.compute_spin_dens(Cube_S_Rho.get_pos(i, j, k), d, phi));
        if (!nodate)
            progress->update();
    }
    if (!nodate)
    {
        delete (progress);

        _time_point end = get_time();
        print_time(start, end, file);
    }
};

void Calc_Prop(
    std::vector<cube> &Cubes,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool test,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    const double radius_bohr = constants::ang2bohr(radius);
    const vector<atom> atoms = wavy.get_atoms();
    cube rho_contrib(Cubes[cube_type::Rho]);

    rho_contrib.evaluate_on_grid(
        [&](const d3 &pos_grid, const i3 &, const i3 &mapped_idx) {
            if (!is_within_radius(pos_grid, atoms, radius_bohr))
                return 0.0;

            const PropValues values = compute_prop_values(Cubes, wavy, pos_grid);
            accumulate_prop_values(Cubes, mapped_idx, values);
            return values.rho;
        },
        wrap);

#pragma omp parallel for schedule(dynamic)
    for (int x = 0; x < Cubes[cube_type::Rho].get_size(0); x++)
        for (int y = 0; y < Cubes[cube_type::Rho].get_size(1); y++)
            for (int z = 0; z < Cubes[cube_type::Rho].get_size(2); z++)
                Cubes[cube_type::Rho].set_value(
                    x,
                    y,
                    z,
                    Cubes[cube_type::Rho].get_value(x, y, z) + rho_contrib.get_value(x, y, z));

    if (!test)
    {
        _time_point end = get_time();
        print_time(start, end, file);
    }
};

void Calc_ESP(
    cube &CubeESP,
    const WFN &wavy,
    double radius,
    bool no_date,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();

    vec2 d2;
    d2.resize(wavy.get_ncen());
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        d2[i].resize(wavy.get_ncen(), 0.0);
        for (int j = 0; j < wavy.get_ncen(); j++)
        {
            if (i == j)
            {
                d2[i][j] = 0;
                continue;
            }
            d2[i][j] = pow(wavy.get_atom_coordinate(i, 0) - wavy.get_atom_coordinate(j, 0), 2) + pow(wavy.get_atom_coordinate(i, 1) - wavy.get_atom_coordinate(j, 1), 2) + pow(wavy.get_atom_coordinate(i, 2) - wavy.get_atom_coordinate(j, 2), 2);
        }
    }
    const double radius_bohr = constants::ang2bohr(radius);
    const vector<atom> atoms = wavy.get_atoms();

    evaluate_cube_in_radius(
        CubeESP,
        wrap,
        atoms,
        radius_bohr,
        [&](const d3 &pos) {
            return wavy.computeESP(pos, d2);
        });

    if (!no_date)
    {
        _time_point end = get_time();
        print_time(start, end, file);
    }
};

void Calc_MO(
    cube &CubeMO,
    int mo,
    const WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    err_checkf(mo <= wavy.get_nmo(), to_string(mo) + " bigger MO selected than " + to_string(wavy.get_nmo()) + " contained in the wavefunctions!", file);
    _time_point start = get_time();
    const double radius_bohr = constants::ang2bohr(radius);
    const vector<atom> atoms = wavy.get_atoms();

    evaluate_cube_in_radius(
        CubeMO,
        wrap,
        atoms,
        radius_bohr,
        [&](const d3 &pos) {
            return wavy.computeMO(pos, mo);
        });

    _time_point end = get_time();
    print_time(start, end, file);
};

namespace {

struct PromolecularFragmentDensities {
    double rho1 = 0.0;
    double rho2 = 0.0;

    double total() const
    {
        return rho1 + rho2;
    }
};

struct PromolecularAtom {
    d3 pos;
    int charge = 1;
    int fragment = 0;
};

std::vector<Thakkar> make_thakkar_interpolators()
{
    std::vector<Thakkar> atom_models;
    atom_models.reserve(92);
    for (int a = 0; a < 92; a++)
    {
        atom_models.emplace_back(a);
        atom_models[a].make_interpolator(1.005 * 1.005 * 1.005, 1E-7);
    }
    return atom_models;
}

double cube_value_clamped(const cube &source, int x, int y, int z)
{
    x = std::clamp(x, 0, source.get_size(0) - 1);
    y = std::clamp(y, 0, source.get_size(1) - 1);
    z = std::clamp(z, 0, source.get_size(2) - 1);
    return source.get_value(x, y, z);
}

double cube_axis_step(const cube &source, int axis)
{
    double sum = 0.0;
    for (int i = 0; i < 3; i++)
        sum += source.get_vector(i, axis) * source.get_vector(i, axis);
    return std::sqrt(sum);
}

PromolecularFragmentDensities promolecular_fragment_densities_at(
    const d3 &pos,
    const std::vector<PromolecularAtom> &atoms,
    const std::vector<Thakkar> &atom_models)
{
    PromolecularFragmentDensities result;
    for (const PromolecularAtom &atom : atoms)
    {
        const double contribution = atom_models[atom.charge - 1].get_interpolated_density(array_length(pos, atom.pos));
        if (atom.fragment == 1)
            result.rho1 += contribution;
        else
            result.rho2 += contribution;
    }
    return result;
}

bool is_promolecular_nci_point(
    const PromolecularFragmentDensities &densities,
    double total_density,
    double dominant_density_cutoff,
    double fragment_sum_cutoff)
{
    total_density = std::abs(total_density);
    if (total_density <= 1E-20)
        return false;

    const double fragment_density = densities.total();
    if (fragment_density <= 1E-20)
        return false;

    if (std::max(densities.rho1, densities.rho2) >= fragment_density * dominant_density_cutoff)
        return false;

    return fragment_density >= total_density * fragment_sum_cutoff;
}

double reduced_density_gradient_at(const cube &rho_cube, int x, int y, int z)
{
    const double rho = std::abs(rho_cube.get_value(x, y, z));
    if (rho <= 1E-20)
        return 0.0;

    const double hx = cube_axis_step(rho_cube, 0);
    const double hy = cube_axis_step(rho_cube, 1);
    const double hz = cube_axis_step(rho_cube, 2);

    const double gx = (cube_value_clamped(rho_cube, x + 1, y, z) - cube_value_clamped(rho_cube, x - 1, y, z)) / (2.0 * hx);
    const double gy = (cube_value_clamped(rho_cube, x, y + 1, z) - cube_value_clamped(rho_cube, x, y - 1, z)) / (2.0 * hy);
    const double gz = (cube_value_clamped(rho_cube, x, y, z + 1) - cube_value_clamped(rho_cube, x, y, z - 1)) / (2.0 * hz);

    const double grad_norm = std::sqrt(gx * gx + gy * gy + gz * gz);
    const double rdg_factor = 2.0 * std::pow(3.0 * constants::PI * constants::PI, 1.0 / 3.0);
    return grad_norm / (rdg_factor * std::pow(rho, 4.0 / 3.0));
}

double lambda2_at(const cube &rho_cube, int x, int y, int z)
{
    const double hx = cube_axis_step(rho_cube, 0);
    const double hy = cube_axis_step(rho_cube, 1);
    const double hz = cube_axis_step(rho_cube, 2);
    const double center = cube_value_clamped(rho_cube, x, y, z);

    double hessian[9]{};
    hessian[0] = (cube_value_clamped(rho_cube, x + 1, y, z) - 2.0 * center + cube_value_clamped(rho_cube, x - 1, y, z)) / (hx * hx);
    hessian[4] = (cube_value_clamped(rho_cube, x, y + 1, z) - 2.0 * center + cube_value_clamped(rho_cube, x, y - 1, z)) / (hy * hy);
    hessian[8] = (cube_value_clamped(rho_cube, x, y, z + 1) - 2.0 * center + cube_value_clamped(rho_cube, x, y, z - 1)) / (hz * hz);

    hessian[1] = hessian[3] =
        (cube_value_clamped(rho_cube, x + 1, y + 1, z) -
            cube_value_clamped(rho_cube, x + 1, y - 1, z) -
            cube_value_clamped(rho_cube, x - 1, y + 1, z) +
            cube_value_clamped(rho_cube, x - 1, y - 1, z)) / (4.0 * hx * hy);
    hessian[2] = hessian[6] =
        (cube_value_clamped(rho_cube, x + 1, y, z + 1) -
            cube_value_clamped(rho_cube, x + 1, y, z - 1) -
            cube_value_clamped(rho_cube, x - 1, y, z + 1) +
            cube_value_clamped(rho_cube, x - 1, y, z - 1)) / (4.0 * hx * hz);
    hessian[5] = hessian[7] =
        (cube_value_clamped(rho_cube, x, y + 1, z + 1) -
            cube_value_clamped(rho_cube, x, y + 1, z - 1) -
            cube_value_clamped(rho_cube, x, y - 1, z + 1) +
            cube_value_clamped(rho_cube, x, y - 1, z - 1)) / (4.0 * hy * hz);

    return get_lambda_1(hessian);
}

void add_promolecular_atoms(
    const WFN &fragment,
    int fragment_id,
    std::vector<PromolecularAtom> &atoms)
{
    for (int i = 0; i < fragment.get_ncen(); i++)
    {
        PromolecularAtom atom_entry;
        atom_entry.pos = fragment.get_atom_pos(i);
        atom_entry.charge = fragment.get_atom_charge(i);
        atom_entry.fragment = fragment_id;
        atoms.push_back(atom_entry);
    }
}

void add_atoms_to_combined_wfn(const WFN &fragment, WFN &combined)
{
    for (int i = 0; i < fragment.get_ncen(); i++)
        combined.push_back_atom(
            fragment.get_atom_label(i),
            fragment.get_atom_coordinate(i, 0),
            fragment.get_atom_coordinate(i, 1),
            fragment.get_atom_coordinate(i, 2),
            fragment.get_atom_charge(i));
}

std::string tcl_quote_path(const std::filesystem::path &path)
{
    const std::string input = path.generic_string();
    std::string output = "\"";
    for (const char c : input)
    {
        if (c == '\\' || c == '"' || c == '$' || c == '[' || c == ']')
            output += '\\';
        output += c;
    }
    output += '"';
    return output;
}

void write_promolecular_nci_vmd(
    const std::filesystem::path &xyz1,
    const std::filesystem::path &xyz2,
    const std::filesystem::path &output_base,
    std::ostream &log)
{
    const std::filesystem::path signed_rho_path = output_base.string() + "_signed_rho.cube";
    const std::filesystem::path rdg_path = output_base.string() + "_rdg.cube";
    const std::filesystem::path vmd_path = output_base.string() + "_nci.vmd";

    std::ofstream vmd_file(vmd_path, std::ios::out);
    err_checkf(vmd_file.good(), "Could not open " + vmd_path.string() + " for writing.", log);

    vmd_file
        << "# VMD visualization for promolecular NCI/RDG analysis\n"
        << "# Load with: vmd -e " << vmd_path.filename().string() << "\n"
        << "display projection Orthographic\n"
        << "display depthcue off\n"
        << "axes location Off\n"
        << "color Display Background white\n"
        << "color scale method BGR\n"
        << "\n"
        << "mol new " << tcl_quote_path(std::filesystem::absolute(xyz1)) << " type xyz waitfor all\n"
        << "set frag1 [molinfo top]\n"
        << "mol delrep 0 $frag1\n"
        << "mol representation CPK\n"
        << "mol color Name\n"
        << "mol selection all\n"
        << "mol material Opaque\n"
        << "mol addrep $frag1\n"
        << "\n"
        << "mol new " << tcl_quote_path(std::filesystem::absolute(xyz2)) << " type xyz waitfor all\n"
        << "set frag2 [molinfo top]\n"
        << "mol delrep 0 $frag2\n"
        << "mol representation CPK\n"
        << "mol color Name\n"
        << "mol selection all\n"
        << "mol material Opaque\n"
        << "mol addrep $frag2\n"
        << "\n"
        << "mol new " << tcl_quote_path(std::filesystem::absolute(signed_rho_path)) << " type cube first 0 last -1 step 1 waitfor 1 volsets {0 }\n"
        << "set nci [molinfo top]\n"
        << "mol delrep 0 $nci\n"
        << "mol addfile " << tcl_quote_path(std::filesystem::absolute(rdg_path)) << " type cube first 0 last -1 step 1 waitfor 1 volsets {0 } $nci\n"
        << "while {[molinfo $nci get numreps] > 0} {\n"
        << "    mol delrep 0 $nci\n"
        << "}\n"
        << "mol addrep $nci\n"
        << "mol modselect 0 $nci all\n"
        << "mol modstyle 0 $nci Isosurface 0.5 1 0 0 1 1\n"
        << "mol modcolor 0 $nci Volume 0\n"
        << "mol modmaterial 0 $nci Opaque\n"
        << "mol colupdate 0 $nci on\n"
        << "mol scaleminmax $nci 0 -2.0 2.0\n"
        << "color scale method BGR\n"
        << "color scale midpoint 0.5\n"
        << "\n"
        << "display resetview\n";

    log << "Wrote " << vmd_path << std::endl;
}

void write_promolecular_nci_plot_script(
    const std::filesystem::path &output_base,
    const properties_options &opts,
    std::ostream &log)
{
    const std::filesystem::path values_path = output_base.string() + "_values.dat";
    const std::filesystem::path plot_path = output_base.string() + "_plot.py";

    std::ofstream plot_file(plot_path, std::ios::out);
    err_checkf(plot_file.good(), "Could not open " + plot_path.string() + " for writing.", log);

    plot_file
        << "from pathlib import Path\n"
        << "import sys\n"
        << "\n"
        << "import matplotlib.pyplot as plt\n"
        << "import numpy as np\n"
        << "from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm\n"
        << "\n"
        << "def read_float_arg(name, default):\n"
        << "    if name not in sys.argv[1:]:\n"
        << "        return default\n"
        << "    idx = sys.argv.index(name)\n"
        << "    if idx + 1 >= len(sys.argv):\n"
        << "        raise SystemExit(f\"{name} requires a numeric value\")\n"
        << "    return float(sys.argv[idx + 1])\n"
        << "\n"
        << "dat_path = Path(__file__).with_name(" << std::quoted(values_path.filename().string()) << ")\n"
        << "data = np.loadtxt(dat_path, comments=\"#\")\n"
        << "if data.ndim == 1:\n"
        << "    data = data.reshape(1, -1)\n"
        << "\n"
        << "signed_rho = data[:, 0]\n"
        << "rdg = data[:, 1]\n"
        << "xmax = " << std::setprecision(16) << opts.promol_nci_rho_abs_max << "\n"
        << "ymax = " << std::setprecision(16) << opts.promol_nci_rdg_max << "\n"
        << "xmax = read_float_arg(\"-xmax\", None if xmax < 0.0 else xmax)\n"
        << "ymax = read_float_arg(\"-ymax\", None if ymax < 0.0 else ymax)\n"
        << "mask = np.ones_like(signed_rho, dtype=bool)\n"
        << "if xmax is not None:\n"
        << "    mask &= np.abs(signed_rho) <= xmax\n"
        << "if ymax is not None:\n"
        << "    mask &= rdg <= ymax\n"
        << "signed_rho = signed_rho[mask]\n"
        << "rdg = rdg[mask]\n"
        << "if signed_rho.size == 0:\n"
        << "    raise SystemExit(\"No points remain after applying plot limits\")\n"
        << "rho_min = float(np.min(signed_rho))\n"
        << "rho_max = float(np.max(signed_rho))\n"
        << "\n"
        << "cmap = LinearSegmentedColormap.from_list(\"nci_bgr\", [\"blue\", \"green\", \"red\"])\n"
        << "if rho_min < 0.0 < rho_max:\n"
        << "    norm = TwoSlopeNorm(vmin=rho_min, vcenter=0.0, vmax=rho_max)\n"
        << "else:\n"
        << "    norm = None\n"
        << "\n"
        << "fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)\n"
        << "scatter = ax.scatter(signed_rho, rdg, c=signed_rho, s=4, cmap=cmap, norm=norm, linewidths=0)\n"
        << "ax.axvline(0.0, color=\"0.65\", linewidth=0.8)\n"
        << "ax.set_xlabel(\"signed rho\")\n"
        << "ax.set_ylabel(\"RDG\")\n"
        << "if xmax is not None:\n"
        << "    ax.set_xlim(-xmax, xmax)\n"
        << "if ymax is not None:\n"
        << "    ax.set_ylim(0.0, ymax)\n"
        << "ax.set_title(dat_path.stem.replace(\"_values\", \"\"))\n"
        << "cbar = fig.colorbar(scatter, ax=ax)\n"
        << "cbar.set_label(\"signed rho\")\n"
        << "png_path = dat_path.with_name(dat_path.stem.replace(\"_values\", \"_plot\") + \".png\")\n"
        << "fig.savefig(png_path, dpi=300)\n"
        << "print(f\"Wrote {png_path}\")\n"
        << "if \"-show\" in sys.argv[1:]:\n"
        << "    plt.show()\n"
        << "else:\n"
        << "    plt.close(fig)\n";

    log << "Wrote " << plot_path << std::endl;
}

} // namespace

void promolecular_nci_analysis(
    const std::filesystem::path &xyz1,
    const std::filesystem::path &xyz2,
    const properties_options &opts,
    std::ostream &log)
{
    using namespace std;

    err_checkf(std::filesystem::exists(xyz1), "First XYZ file does not exist: " + xyz1.string(), log);
    err_checkf(std::filesystem::exists(xyz2), "Second XYZ file does not exist: " + xyz2.string(), log);

    WFN fragment1(e_origin::xyz);
    WFN fragment2(e_origin::xyz);
    fragment1.read_xyz(xyz1, log, false);
    fragment2.read_xyz(xyz2, log, false);

    WFN combined(e_origin::xyz);
    combined.set_path(xyz1.parent_path() / (xyz1.stem().string() + "_" + xyz2.stem().string() + ".xyz"));
    add_atoms_to_combined_wfn(fragment1, combined);
    add_atoms_to_combined_wfn(fragment2, combined);

    properties_options local_opts = opts;
    readxyzMinMax_fromWFN(combined, local_opts, true);
    err_checkf(local_opts.NbSteps[0] > 1 && local_opts.NbSteps[1] > 1 && local_opts.NbSteps[2] > 1,
        "Promolecular NCI grid is too small; decrease -resolution or increase -radius.", log);

    const std::filesystem::path output_base =
        xyz1.parent_path() / (xyz1.stem().string() + "_" + xyz2.stem().string());

    vector<PromolecularAtom> atoms;
    atoms.reserve(fragment1.get_ncen() + fragment2.get_ncen());
    add_promolecular_atoms(fragment1, 1, atoms);
    add_promolecular_atoms(fragment2, 2, atoms);
    const vector<Thakkar> atom_models = make_thakkar_interpolators();

    cube rho_cube(local_opts.NbSteps, combined.get_ncen(), true);
    cube signed_rho_cube(local_opts.NbSteps, combined.get_ncen(), true);
    cube rdg_cube(local_opts.NbSteps, combined.get_ncen(), true);
    rho_cube.give_parent_wfn(combined);
    signed_rho_cube.give_parent_wfn(combined);
    rdg_cube.give_parent_wfn(combined);

    for (int i = 0; i < 3; i++)
    {
        const double step = (local_opts.MinMax[i + 3] - local_opts.MinMax[i]) / local_opts.NbSteps[i];
        rho_cube.set_origin(i, local_opts.MinMax[i]);
        signed_rho_cube.set_origin(i, local_opts.MinMax[i]);
        rdg_cube.set_origin(i, local_opts.MinMax[i]);
        rho_cube.set_vector(i, i, step);
        signed_rho_cube.set_vector(i, i, step);
        rdg_cube.set_vector(i, i, step);
    }
    rho_cube.calc_dv();
    signed_rho_cube.calc_dv();
    rdg_cube.calc_dv();

    rho_cube.set_comment1("Promolecular density using Thakkar spherical atoms");
    rho_cube.set_comment2("from " + xyz1.string() + " and " + xyz2.string());
    signed_rho_cube.set_comment1("Promolecular signed density using Thakkar spherical atoms");
    signed_rho_cube.set_comment2("from " + xyz1.string() + " and " + xyz2.string());
    rdg_cube.set_comment1("Promolecular reduced density gradient using Thakkar spherical atoms");
    rdg_cube.set_comment2("from " + xyz1.string() + " and " + xyz2.string());

    log << "Promolecular NCI analysis for " << xyz1 << " and " << xyz2 << endl;
    log << "Grid points: " << local_opts.n_grid_points() << endl;
    log << "Dominant-fragment density discard cutoff: " << opts.promol_nci_rcut1 << endl;
    log << "Fragment-sum density keep cutoff: " << opts.promol_nci_rcut2 << endl;

    ProgressBar density_progress(rho_cube.get_size(0), 50, "=", " ", "Calculating promolecular rho");
#pragma omp parallel for schedule(dynamic)
    for (int x = 0; x < rho_cube.get_size(0); x++)
    {
        for (int y = 0; y < rho_cube.get_size(1); y++)
            for (int z = 0; z < rho_cube.get_size(2); z++)
            {
                const d3 pos = rho_cube.get_pos(x, y, z);
                const PromolecularFragmentDensities densities = promolecular_fragment_densities_at(pos, atoms, atom_models);
                const double rho = densities.total();
                rho_cube.set_value(x, y, z, rho);
                rdg_cube.set_value(
                    x,
                    y,
                    z,
                    is_promolecular_nci_point(densities, rho, opts.promol_nci_rcut1, opts.promol_nci_rcut2) ? 1.0 : 101.0);
            }
        density_progress.update();
    }
    std::cout << std::endl;

    i3 shrink_lower;
    i3 shrink_upper;
    if (rdg_cube.find_value_bounds(shrink_lower, shrink_upper, 101.0))
    {
        rho_cube.shrink_to_bounds(shrink_lower, shrink_upper);
        signed_rho_cube.shrink_to_bounds(shrink_lower, shrink_upper);
        rdg_cube.shrink_to_bounds(shrink_lower, shrink_upper);
        log << "Shrunk promolecular NCI work grid to "
            << rdg_cube.get_size(0) << " x "
            << rdg_cube.get_size(1) << " x "
            << rdg_cube.get_size(2)
            << " grid points, ignoring RDG mask value 101." << endl;
    }

    ofstream values_file(output_base.string() + "_values.dat", ios::out);
    err_checkf(values_file.good(), "Could not open " + output_base.string() + "_values.dat for writing.", log);
    values_file << "# signed_rho rdg\n";
    values_file << "# rcut1 " << opts.promol_nci_rcut1 << " rcut2 " << opts.promol_nci_rcut2 << "\n";

    unsigned long long kept_points = 0;
    std::vector<std::ostringstream> values_by_thread(1);
#ifdef _OPENMP
    values_by_thread.resize(omp_get_max_threads());
#endif
    ProgressBar rdg_progress(rho_cube.get_size(0), 50, "=", " ", "Calculating promolecular RDG");
#pragma omp parallel reduction(+ : kept_points)
    {
        int thread_id = 0;
#ifdef _OPENMP
        thread_id = omp_get_thread_num();
#endif
        std::ostringstream &local_values = values_by_thread[thread_id];
#pragma omp for schedule(dynamic)
    for (int x = 0; x < rho_cube.get_size(0); x++)
    {
        for (int y = 0; y < rho_cube.get_size(1); y++)
        {
            for (int z = 0; z < rho_cube.get_size(2); z++)
            {
                const double rho = rho_cube.get_value(x, y, z);
                if (std::abs(rdg_cube.get_value(x, y, z) - 101.0) <= 1E-12)
                {
                    signed_rho_cube.set_value(x, y, z, 0.0);
                    rdg_cube.set_value(x, y, z, 101.0);
                    continue;
                }

                const double lambda2 = lambda2_at(rho_cube, x, y, z);
                const double signed_rho = lambda2 < 0.0 ? -rho : rho;
                const double rdg = sanitize_finite(reduced_density_gradient_at(rho_cube, x, y, z));

                signed_rho_cube.set_value(x, y, z, signed_rho);
                rdg_cube.set_value(x, y, z, rdg);

                if (opts.promol_nci_rho_abs_max >= 0.0 && std::abs(signed_rho) > opts.promol_nci_rho_abs_max)
                    continue;
                if (opts.promol_nci_rdg_max >= 0.0 && rdg > opts.promol_nci_rdg_max)
                    continue;

                local_values << scientific << setprecision(10) << signed_rho << " " << rdg << "\n";
                kept_points++;
            }
        }
        rdg_progress.update();
    }
    }
    for (const std::ostringstream &local_values : values_by_thread)
        values_file << local_values.str();

    signed_rho_cube.set_path(output_base.string() + "_signed_rho.cube");
    rdg_cube.set_path(output_base.string() + "_rdg.cube");
    signed_rho_cube.write_file(true);
    rdg_cube.write_file(true);
    write_promolecular_nci_vmd(xyz1, xyz2, output_base, log);
    write_promolecular_nci_plot_script(output_base, opts, log);

    log << "Wrote " << signed_rho_cube.get_path() << endl;
    log << "Wrote " << rdg_cube.get_path() << endl;
    log << "Wrote " << output_base.string() + "_values.dat" << " with " << kept_points << " points" << endl;
}

void properties_calculation(options &opt)
{
    using namespace std;
    ofstream log2("NoSpherA2_cube.log", ios::out);
    auto _coutbuf = std::cout.rdbuf(log2.rdbuf()); // save and redirect
    log2 << NoSpherA2_message(opt.no_date);
    if (!opt.no_date)
    {
        log2 << build_date;
    }
    log2.flush();

    err_checkf(opt.wfn != "", "Error, no wfn file specified!", log2);
    WFN wavy(opt.wfn);
    if (opt.debug)
        log2 << "Starting calculation of properties" << endl;
    if (opt.properties.all_mos)
        for (int mo = 0; mo < wavy.get_nmo(); mo++)
            opt.properties.MO_numbers.push_back(mo);
    if (opt.debug)
        log2 << "Size of MOs: " << opt.properties.MO_numbers.size() << endl;

    vec2 cell_matrix;
    cell_matrix.resize(3);
    for (int i = 0; i < 3; i++)
        cell_matrix[i].resize(3, 0.0);
    if (opt.debug)
    {
        log2 << "cif|resolution|res(bohr)|radius|rad(bohr): " << opt.cif << "|" << opt.properties.resolution << "|" << constants::ang2bohr(opt.properties.resolution) << "|" << opt.properties.radius << "|" << constants::ang2bohr(opt.properties.radius) << endl;
        for (int a = 0; a < wavy.get_ncen(); a++)
            log2 << "Atom " << a << " at " << wavy.get_atom_coordinate(a, 0) << " " << wavy.get_atom_coordinate(a, 1) << " " << wavy.get_atom_coordinate(a, 2) << endl;
    }
    if (opt.cif != "")
        readxyzMinMax_fromCIF(opt.cif, opt.properties, cell_matrix);
    else
    {
        readxyzMinMax_fromWFN(wavy, opt.properties, true);
        for (int i = 0; i < 3; i++)
            cell_matrix[i][i] = constants::ang2bohr(opt.properties.resolution);
    }
    if (opt.debug)
    {
        log2 << "MinMax: ";
        for (int i = 0; i < 6; i++)
            log2 << setw(14) << scientific << opt.properties.MinMax[i];
        log2 << endl;
        log2 << "Steps: ";
        for (int i = 0; i < 3; i++)
            log2 << setw(14) << scientific << opt.properties.NbSteps[i];
        log2 << endl;
        log2 << "Cell Matrix:" << endl;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
                log2 << setw(14) << scientific << cell_matrix[i][j];
            log2 << endl;
        }
    }
    std::vector<cube> cubes;
    cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), true);
    cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.rdg);
    cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.elf);
    cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.eli);
    cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.lap);
    cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.esp);
    cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), true);
    cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.hdef);
    cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.def);
    cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.hirsh);
    cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.hirsh);
    cubes.emplace_back(opt.properties.NbSteps, wavy.get_ncen(), opt.properties.hdef || opt.properties.hirsh);

    for (cube &cube : cubes)
        cube.give_parent_wfn(wavy);


    for (int i = 0; i < 3; i++)
    {
        for (cube &cube : cubes)
            cube.set_origin(i, opt.properties.MinMax[i]);

        for (int j = 0; j < 3; j++)
        {
            for (cube &cube : cubes)
                cube.set_vector(i, j, cell_matrix[i][j]);
        }
    }
    for (cube &cube : cubes)
        cube.calc_dv();
    if (opt.debug)
        log2 << "Origins etc. are set up" << endl;
    cubes[cube_type::Rho].set_comment1("Calculated density using NoSpherA2");
    cubes[cube_type::RDG].set_comment1("Calculated reduced density gradient using NoSpherA2");
    cubes[cube_type::Elf].set_comment1("Calculated electron localization function using NoSpherA2");
    cubes[cube_type::Eli].set_comment1("Calculated same-spin electron localizability indicator using NoSpherA2");
    cubes[cube_type::Lap].set_comment1("Calculated laplacian of electron density using NoSpherA2");
    cubes[cube_type::ESP].set_comment1("Calculated electrostatic potential using NoSpherA2");
    cubes[cube_type::MO_val].set_comment1("Calcualted MO values using NoSpherA2");
    cubes[cube_type::HDEF].set_comment1("Calculated Atomic Hirshfeld deformation density values using NoSpherA2");
    cubes[cube_type::DEF].set_comment1("Calculated static deformation density values using NoSpherA2");
    cubes[cube_type::Hirsh].set_comment1("Calculated Hirshfeld atom density values using NoSpherA2");
    cubes[cube_type::Spin_Density].set_comment1("Calculated spin density using NoSpherA2");
    for (auto cube : cubes)
        cube.set_comment2("from " + wavy.get_path().string());

    cubes[cube_type::Rho].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rho.cube");
    cubes[cube_type::RDG].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rdg.cube");
    cubes[cube_type::Elf].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_elf.cube");
    cubes[cube_type::Eli].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_eli.cube");
    cubes[cube_type::Lap].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_lap.cube");
    cubes[cube_type::ESP].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_esp.cube");
    cubes[cube_type::DEF].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_def.cube");
    cubes[cube_type::Hirsh].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_hirsh.cube");
    cubes[cube_type::Spin_Density].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_s_rho.cube");

    log2 << "\nCalculating:" << endl;
    if (opt.properties.hdef || opt.properties.def || opt.properties.hirsh)
        log2 << "Rho, ";
    if (opt.properties.hdef || opt.properties.hirsh)
        log2 << "Spherical Rho, ";
    if (opt.properties.def)
        log2 << "Static deformation density, ";
    if (opt.properties.hdef)
        log2 << "Hirshfeld deformation density, ";
    if (opt.properties.hirsh)
        log2 << "Hirshfeld density, ";
    if (opt.properties.lap)
        log2 << "Laplacian, ";
    if (opt.properties.eli)
        log2 << "ELI, ";
    if (opt.properties.elf)
        log2 << "ELF, ";
    if (opt.properties.rdg)
        log2 << "RDG, ";
    if (opt.properties.esp)
        log2 << "ESP, ";
    if (opt.properties.MO_numbers.size() != 0)
        log2 << "MOs, ";
    if (opt.properties.s_rho)
        log2 << "Spin density, ";
    log2 << endl;

    log2 << "Calculating for " << fixed << setprecision(0) << opt.properties.NbSteps[0] * opt.properties.NbSteps[1] * opt.properties.NbSteps[2] << " Gridpoints." << endl;

    Calc_Rho(cubes[cube_type::Rho], wavy, opt.properties.radius, log2, opt.cif != "");

    if (opt.properties.integral_accuracy != -1) {
        log2 << "Refining grid files to integral accuracy of " << opt.properties.integral_accuracy << " ..." << flush;
        vec2 d(16, vec(wavy.get_ncen(), 0.0));
        vec phi(wavy.get_nmo(true), 0.0);
        cubes[cube_type::Rho].adaptive_refine([&wavy](const d3 &pos) { return wavy.compute_dens(pos); }, opt.properties.integral_accuracy, 8, 3);
    }

    if (opt.properties.MO_numbers.size() != 0)
        for (int i = 0; i < opt.properties.MO_numbers.size(); i++)
        {
            log2 << "Calcualting MO: " << opt.properties.MO_numbers[i] << endl;
            cubes[cube_type::MO_val].set_zero();
            cubes[cube_type::MO_val].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_MO_" + to_string(opt.properties.MO_numbers[i]) + ".cube");
            Calc_MO(cubes[cube_type::MO_val], opt.properties.MO_numbers[i], wavy, opt.properties.radius, log2, opt.cif != "");
            cubes[cube_type::MO_val].write_file(true);
        }

    wavy.delete_unoccupied_MOs();
    wavy.delete_Qs();

    for (int i = 1; i < cubes.size(); i++) {
        if (cubes[i].get_loaded())
            cubes[i].resize(cubes[cube_type::Rho].get_sizes());
        cubes[i].set_vectors(cubes[cube_type::Rho].get_vectors());
    }

    if (opt.properties.hdef || opt.properties.def || opt.properties.hirsh)
    {
        cubes[cube_type::spherical_density].resize(cubes[cube_type::Rho].get_sizes());
        for (int i = 0; i < 3; i++)
        {
            cubes[cube_type::spherical_density].set_origin(i, opt.properties.MinMax[i]);
            for (int j = 0; j < 3; j++)
                cubes[cube_type::spherical_density].set_vector(i, j, cell_matrix[i][j]);
        }
        if (opt.properties.hdef || opt.properties.hirsh)
        {
            log2 << "Calcualting spherical Rho...";
            Calc_Spherical_Dens(cubes[cube_type::spherical_density], wavy, opt.properties.radius, log2, opt.cif != "");
            log2 << " ...done!" << endl;
        }

        if (opt.properties.def)
        {
            log2 << "Calculating static deformation density...";
            if (opt.properties.hdef)
                Calc_Static_Def(cubes, wavy, opt.properties.radius, log2, opt.cif != "");
            else
                Calc_Static_Def(cubes, wavy, opt.properties.radius, log2, opt.cif != "");
            log2 << " ...done!" << endl;
        }

        if (opt.properties.hdef)
        {
            for (int a = 0; a < wavy.get_ncen(); a++)
            {
                log2 << "Calcualting Hirshfeld deformation density for atom: " << a << endl;
                cubes[cube_type::HDEF].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_HDEF_" + to_string(a) + ".cube");
                Calc_Hirshfeld(cubes, wavy, opt.properties.radius, a, log2, opt.cif != "");
                cubes[cube_type::HDEF].write_file(true);
                cubes[cube_type::HDEF].set_zero();
            }
        }

        if (opt.properties.hirsh)
        {
            log2 << "Calcualting Hirshfeld density for atom: " << opt.properties.hirsh_number << endl;
            Calc_Hirshfeld_atom(cubes, wavy, opt.properties.radius, opt.properties.hirsh_number, log2, opt.cif != "");
            log2 << "..done!" << endl;
        }
    }

    if (opt.properties.lap || opt.properties.eli || opt.properties.elf || opt.properties.rdg || opt.properties.esp)
        Calc_Prop(cubes, wavy, opt.properties.radius, log2, opt.no_date, opt.cif != "");

    if (opt.properties.s_rho)
        Calc_S_Rho(cubes[cube_type::Spin_Density], wavy, log2, opt.no_date);

    log2 << "Writing cubes to Disk..." << flush;
    if (opt.properties.rdg)
    {
        cubes[cube_type::Rho].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_signed_rho.cube");
        cubes[cube_type::Rho].write_file(true);
        cubes[cube_type::Rho].set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rho.cube");
        cubes[cube_type::Rho].write_file(true, true);
    }
    else if (opt.properties.lap || opt.properties.eli || opt.properties.elf || opt.properties.esp)
        cubes[cube_type::Rho].write_file(true);
    if (opt.properties.rdg)
        cubes[cube_type::RDG].write_file(true);
    if (opt.properties.lap)
        cubes[cube_type::Lap].write_file(true);
    if (opt.properties.elf)
        cubes[cube_type::Elf].write_file(true);
    if (opt.properties.eli)
        cubes[cube_type::Eli].write_file(true);
    if (opt.properties.s_rho)
        cubes[cube_type::Spin_Density].write_file(true);
    if (opt.properties.def)
    {
        cubes[cube_type::DEF].write_file(true);
        cubes[cube_type::Rho].write_file(true);
    }
    if (opt.properties.hirsh)
        cubes[cube_type::Hirsh].write_file(true);

    log2 << " done!" << endl;

    if (opt.properties.esp)
    {
        log2 << "Calculating ESP..." << flush;
        WFN temp = wavy;
        temp.delete_unoccupied_MOs();
        temp.delete_Qs();
        Calc_ESP(cubes[cube_type::ESP], temp, opt.properties.radius, opt.no_date, log2);
        log2 << "Writing cube to Disk..." << flush;
        cubes[cube_type::ESP].write_file(true);
        log2 << "  done!" << endl;
    }
    // return output tostd::cout
    std::cout.rdbuf(_coutbuf);
    log2.close();
    std::cout << "Properties calculation done!" << std::endl;
}

void do_combine_mo(options &opt)
{
    using namespace std;
    WFN wavy1(e_origin::wfn);
    WFN wavy2(e_origin::wfn);
    WFN wavy3(e_origin::wfn);
    wavy1.read_wfn(opt.combine_mo[0], false, std::cout);
    wavy2.read_wfn(opt.combine_mo[1], false, std::cout);
    for (int i = 0; i < wavy1.get_ncen(); i++)
    {
        wavy3.push_back_atom(wavy1.get_atom(i));
    }
    for (int i = 0; i < wavy2.get_ncen(); i++)
    {
        wavy3.push_back_atom(wavy2.get_atom(i));
    }
    std::cout << "In total we have " << wavy3.get_ncen() << " atoms" << endl;

    readxyzMinMax_fromWFN(wavy1, opt.properties, true);
    properties_options prop2 = opt.properties;
    readxyzMinMax_fromWFN(wavy2, prop2, true);

    std::cout << "Read input\nCalculating for MOs ";
    for (int v1 = 0; v1 < opt.cmo1.size(); v1++)
    {
        std::cout << opt.cmo1[v1] << " ";
    }
    std::cout << "of fragment 1 and MOs ";
    for (int v1 = 0; v1 < opt.cmo2.size(); v1++)
    {
        std::cout << opt.cmo2[v1] << " ";
    }
    std::cout << "of fragment 2" << endl;
    std::array<double, 6> MinMax = { 100, 100, 100, -100, -100, -100 };
    std::array<int, 3> steps = { 0, 0, 0 };
    for (int i = 0; i < 3; i++)
    {
        if (opt.properties.MinMax[i] < MinMax[i])
            MinMax[i] = opt.properties.MinMax[i];
        if (opt.properties.MinMax[i + 3] > MinMax[i + 3])
            MinMax[i + 3] = opt.properties.MinMax[i + 3];
    }
    for (int i = 0; i < 3; i++)
    {
        if (prop2.MinMax[i] < MinMax[i])
            MinMax[i] = prop2.MinMax[i];
        if (prop2.MinMax[i + 3] > MinMax[i + 3])
            MinMax[i + 3] = prop2.MinMax[i + 3];
        steps[i] = (int)ceil(constants::bohr2ang(MinMax[i + 3] - MinMax[i]) / 0.1);
    }
    int counter = 0;
    cube total(steps, 0, true);
    cube MO1(steps, 0, true);
    MO1.give_parent_wfn(wavy3);
    MO1.set_na(wavy3.get_ncen());
    cube MO2(steps, 0, true);
    svec fns;
    for (int i = 0; i < 3; i++)
    {
        MO1.set_origin(i, MinMax[i]);
        MO1.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
        total.set_origin(i, MinMax[i]);
        total.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
        MO2.set_origin(i, MinMax[i]);
        MO2.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
    }
    for (int v1 = 0; v1 < opt.cmo1.size(); v1++)
    {
        MO1.set_zero();
        Calc_MO(MO1, opt.cmo1[v1] - 1, wavy1, 40, std::cout);
        for (int j = 0; j < opt.cmo2.size(); j++)
        {
            counter++;
            std::cout << "Running: " << counter << " of " << opt.cmo2.size() * opt.cmo1.size() << endl;
            string filename("");
            MO2.set_zero();
            Calc_MO(MO2, opt.cmo2[j] - 1, wavy2, 40, std::cout);
            std::cout << "writing files..." << flush;
            filename = wavy1.get_path().stem().string() + "_" + std::to_string(opt.cmo1[v1]) + "+" + wavy2.get_path().stem().string() + "_" + std::to_string(opt.cmo2[j]) + ".cube";
            fns.push_back(filename);
            total.set_zero();
            total = MO1;
            total += MO2;
            total.write_file(filename, false);
            filename = wavy1.get_path().stem().string() + "_" + std::to_string(opt.cmo1[v1]) + "-" + wavy2.get_path().stem().string() + "_" + std::to_string(opt.cmo2[j]) + ".cube";
            fns.push_back(filename);
            total.set_zero();
            total = MO1;
            total -= MO2;
            total.write_file(filename, false);
            std::cout << " ... done!" << endl;
        }
    }
    ofstream vmd("read_files.vmd");
    vmd << "mol addrep 0\nmol new {" + fns[0] + "} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 }\n";
    vmd << "animate style Loop\n";
    for (int i = 1; i < fns.size(); i++)
        vmd << "mol addfile {" + fns[i] + "} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 } 0\n";
    vmd << "animate style Loop\ndisplay projection Orthographic\ndisplay depthcue off\n";
    vmd << "axes location Off\ndisplay rendermode GLSL\ncolor Display Background white\ncolor Element P purple\n";
    vmd << "color Element Ni green\ncolor Element C gray\nmol modstyle 0 0 CPK 1.000000 0.300000 12.000000 12.000000\n";
    vmd << "mol modcolor 0 0 Element\nmol color Element\nmol representation CPK 1.000000 0.300000 22.000000 22.000000\n";
    vmd << "mol selection all\nmol material Transparent\nmol addrep 0\nmol modstyle 1 0 Isosurface 0.020000 0 0 0 1 1\n";
    vmd << "mol modcolor 1 0 ColorID 0\nmol selection all\nmol material Transparent\nmol addrep 0\nmol modstyle 2 0 Isosurface -0.020000 0 0 0 1 1\nmol modcolor 2 0 ColorID 1\n";
    vmd << "mol selection all\nmol material Transparent\n";
    vmd.flush();
    vmd.close();
}

static void Calc_Hirshfeld_atom_2(
    cube &CubeHirsh,
    cube &CubeRho,
    cube &CubeSpherical,
    WFN &wavy,
    int _atom,
    std::ostream &file)
{
    (void)file;
    using namespace std;
    Thakkar atom(wavy.get_atom_charge(_atom));

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < CubeHirsh.get_size(0); i++)
    {
        for (int j = 0; j < CubeHirsh.get_size(1); j++)
            for (int k = 0; k < CubeHirsh.get_size(2); k++)
            {

                const d3 PosGrid{ i * CubeHirsh.get_vector(0, 0) + j * CubeHirsh.get_vector(0, 1) + k * CubeHirsh.get_vector(0, 2) + CubeHirsh.get_origin(0),
                                        i * CubeHirsh.get_vector(1, 0) + j * CubeHirsh.get_vector(1, 1) + k * CubeHirsh.get_vector(1, 2) + CubeHirsh.get_origin(1),
                                        i * CubeHirsh.get_vector(2, 0) + j * CubeHirsh.get_vector(2, 1) + k * CubeHirsh.get_vector(2, 2) + CubeHirsh.get_origin(2) };

                // bool skip = true;

                const double dens_choice = atom.get_radial_density(array_length(PosGrid, wavy.get_atom_pos(_atom)));
                const double temp_val = CubeSpherical.get_value(i, j, k);
                if (temp_val != 0)
                    CubeHirsh.set_value(i, j, k, (dens_choice / temp_val * CubeRho.get_value(i, j, k)));
            }
    }
};

enum class dipole_types
{
    atom,
    geometry,
    hirshfeld,
    vdW,
    Unknown
};

dipole_types stringTodipole_types(const std::string &str)
{
    static const std::unordered_map<std::string, dipole_types> stringToEnumMap = {
        {"atom", dipole_types::atom},
        {"geometry", dipole_types::geometry},
        {"hirshfeld", dipole_types::hirshfeld},
        {"vdW", dipole_types::vdW} };

    auto it = stringToEnumMap.find(str);
    if (it != stringToEnumMap.end())
    {
        return it->second;
    }
    else
    {
        return dipole_types::Unknown;
    }
}

vec calc_dipole_for_atom(WFN &wavy, const int &i, cube &Hirshfeld_atom, vec &charges, std::string type = "atom")
{
    double mu_x = 0, mu_y = 0, mu_z = 0;
    double scratch = 0;
    const d3 ax = wavy.get_atom_pos(i);
    double dv = Hirshfeld_atom.get_dv();
    // const int c = wavy.get_atom_charge(i);
    double charge = 0;
    vec origin{ 0, 0, 0 };
    vec bound_atoms;
    for (int j = 0; j < wavy.get_ncen(); j++)
    {
        if (i == j)
            continue;
        const double dist = array_length(ax, wavy.get_atom_pos(j));
        const double svdW = constants::covalent_radii[wavy.get_atom_charge(i)] + constants::covalent_radii[wavy.get_atom_charge(j)];
        if (dist < 1.1 * svdW)
        {
            bound_atoms.push_back(j);
        }
    }
    const double v[9] = { Hirshfeld_atom.get_vector(0, 0), Hirshfeld_atom.get_vector(0, 1), Hirshfeld_atom.get_vector(0, 2),
                         Hirshfeld_atom.get_vector(1, 0), Hirshfeld_atom.get_vector(1, 1), Hirshfeld_atom.get_vector(1, 2),
                         Hirshfeld_atom.get_vector(2, 0), Hirshfeld_atom.get_vector(2, 1), Hirshfeld_atom.get_vector(2, 2) };
    switch (stringTodipole_types(type))
    {
    case dipole_types::atom:
        origin = { ax[0], ax[1], ax[2] };
        break;
    case dipole_types::geometry:
        err_not_impl_f("geometry position not yet implemented", std::cout);
        origin = { 0, 0, 0 };
        break;
    case dipole_types::hirshfeld:
        err_not_impl_f("hirshfeld centers not yet implemented", std::cout);
        break;
    case dipole_types::vdW:
        err_not_impl_f("vdW radius basis not implemented", std::cout);
        break;
    case dipole_types::Unknown:
        err_not_impl_f("Unknown dipole type", std::cout);
        break;
    }
#pragma omp parallel for reduction(+ : mu_x, mu_y, mu_z, charge) private(scratch)
    for (int x = 0; x < Hirshfeld_atom.get_size(0); x++)
    {
        for (int y = 0; y < Hirshfeld_atom.get_size(1); y++)
        {
            for (int z = 0; z < Hirshfeld_atom.get_size(2); z++)
            {
                const d3 PosGrid{
                    x * v[0] + y * v[1] + z * v[2] + Hirshfeld_atom.get_origin(0),
                    x * v[3] + y * v[4] + z * v[5] + Hirshfeld_atom.get_origin(1),
                    x * v[6] + y * v[7] + z * v[8] + Hirshfeld_atom.get_origin(2) };
                scratch = Hirshfeld_atom.get_value(x, y, z) * dv;
                charge += scratch;
                mu_x += (PosGrid[0] - origin[0]) * scratch;
                mu_y += (PosGrid[1] - origin[1]) * scratch;
                mu_z += (PosGrid[2] - origin[2]) * scratch;
            }
        }
    }
    return { mu_x, mu_y, mu_z, charge };
}

void dipole_moments(options &opt, std::ostream &log2)
{
    using namespace std;
    log2 << NoSpherA2_message(opt.no_date);
    if (!opt.no_date)
        log2 << build_date;
    err_checkf(opt.wfn != "", "Error, no wfn file specified!", log2);
    WFN wavy(opt.wfn);
    if (opt.debug)
        log2 << "Starting calculation of dipole moment" << endl;

    if (opt.debug)
        log2 << opt.cif << " " << opt.properties.resolution << " " << opt.properties.radius << endl;
    readxyzMinMax_fromWFN(wavy, opt.properties, true);
    if (opt.debug)
    {
        log2 << "Resolution: " << opt.properties.resolution << endl;
        log2 << "MinMax:" << endl;
        for (int i = 0; i < 6; i++)
            log2 << setw(14) << scientific << opt.properties.MinMax[i];
        log2 << endl;
        log2 << "Steps:" << endl;
        for (int i = 0; i < 3; i++)
            log2 << setw(14) << scientific << opt.properties.NbSteps[i];
        log2 << endl;
    }
    cube Rho(opt.properties.NbSteps, wavy.get_ncen(), true);
    cube SPHER(opt.properties.NbSteps, wavy.get_ncen(), true);

    Rho.give_parent_wfn(wavy);
    SPHER.give_parent_wfn(wavy);
    vec stepsizes{ (opt.properties.MinMax[3] - opt.properties.MinMax[0]) / opt.properties.NbSteps[0],
                  (opt.properties.MinMax[4] - opt.properties.MinMax[1]) / opt.properties.NbSteps[1],
                  (opt.properties.MinMax[5] - opt.properties.MinMax[2]) / opt.properties.NbSteps[2] };

    for (int i = 0; i < 3; i++)
    {
        Rho.set_origin(i, opt.properties.MinMax[i]);
        SPHER.set_origin(i, opt.properties.MinMax[i]);
        Rho.set_vector(i, i, stepsizes[i]);
        SPHER.set_vector(i, i, stepsizes[i]);
    }
    if (opt.debug)
        log2 << "Origins etc are set up" << endl;
    Rho.set_comment1("Calculated density using NoSpherA2");
    SPHER.set_comment1("Calculated Atomic Hirshfeld deformation density values using NoSpherA2");
    Rho.set_comment2("from " + wavy.get_path().string());
    SPHER.set_comment2("from" + wavy.get_path().string());
    Rho.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rho.cube");
    SPHER.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_spher.cube");
    vector<cube> Hirsh(wavy.get_ncen(), Rho);
    vec charges(wavy.get_ncen(), 0);

    log2 << "Calculating for " << fixed << setprecision(0) << opt.properties.NbSteps[0] * opt.properties.NbSteps[1] * opt.properties.NbSteps[2] << " Gridpoints." << endl;

    log2 << "Calcualting Rho...";
    Calc_Rho(Rho, wavy, opt.properties.radius, log2, false);
    log2 << " ...done!\nCalcualting spherical Rho...";
    Calc_Spherical_Dens(SPHER, wavy, opt.properties.radius, log2);
    log2 << " ...done!" << endl;
    vec2 dipole_moments;
    dipole_moments.reserve(wavy.get_ncen());
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        Hirsh[i].calc_dv();
        Hirsh[i].give_parent_wfn(wavy);
        Hirsh[i].set_zero();
        log2 << "Calcualting Hirshfeld density for atom: " << i << endl;
        Calc_Hirshfeld_atom_2(Hirsh[i], Rho, SPHER, wavy, i, log2);
        charges[i] = Hirsh[i].sum();
        log2 << "..done!" << endl;
    }
    for (int i = 0; i < wavy.get_ncen(); i++)
        dipole_moments.emplace_back(calc_dipole_for_atom(wavy, i, Hirsh[i], charges));
    log2 << " atom   |  dipole moment x,        y,         z" << endl
        << "======================================" << endl;
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        log2 << setw(3) << i << " (" << constants::atnr2letter(wavy.get_atom_charge(i)) << ") | " << scientific << setprecision(6) << setw(14) << dipole_moments[i][0] << ", " << setw(14) << dipole_moments[i][1] << ", " << setw(14) << dipole_moments[i][2] << endl;
    }
    std::cout << "\n\nProperties calculation done!" << std::endl;
}

vec2 dipole_moments(WFN &wavy, cube &SPHER, const properties_options &opts, int threads, std::ostream &log2, bool debug)
{
    using namespace std;
    if (debug)
        log2 << "Starting calculation of dipole moment" << endl;
    cube Rho(opts.NbSteps, wavy.get_ncen(), true);

    Rho.give_parent_wfn(wavy);
    vec stepsizes{ (opts.MinMax[3] - opts.MinMax[0]) / opts.NbSteps[0],
                  (opts.MinMax[4] - opts.MinMax[1]) / opts.NbSteps[1],
                  (opts.MinMax[5] - opts.MinMax[2]) / opts.NbSteps[2] };

    for (int i = 0; i < 3; i++)
    {
        Rho.set_origin(i, opts.MinMax[i]);
        Rho.set_vector(i, i, stepsizes[i]);
    }
    if (debug)
        log2 << "Origins etc are set up" << endl;
    Rho.set_comment1("Calculated density using NoSpherA2");
    Rho.set_comment2("from " + wavy.get_path().string());
    Rho.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rho.cube");
    vector<cube> Hirsh(wavy.get_ncen(), Rho);
    vec charges(wavy.get_ncen(), 0);

    log2 << "Calculating for " << fixed << setprecision(0) << opts.n_grid_points() << " Gridpoints." << endl;

    log2 << "Calcualting Rho...";
    Calc_Rho(Rho, wavy, opts.radius, log2, false);
    log2 << " ...done!\nCalcualting spherical Rho...";
    Calc_Spherical_Dens(SPHER, wavy, opts.radius, log2);
    log2 << " ...done!" << endl;
    vec2 dipole_moments;
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        Hirsh[i].calc_dv();
        Hirsh[i].give_parent_wfn(wavy);
        Hirsh[i].set_zero();
        log2 << "Calcualting Hirshfeld density for atom: " << i << endl;
        Calc_Hirshfeld_atom_2(Hirsh[i], Rho, SPHER, wavy, i, log2);
        charges[i] = Hirsh[i].sum();
        log2 << "..done!" << endl;
    }
    for (int i = 0; i < wavy.get_ncen(); i++)
        dipole_moments.push_back(calc_dipole_for_atom(wavy, i, Hirsh[i], charges));
    log2 << "...done!" << endl;
    log2 << " atom   |    charge    | dipole moment x,        y,         z" << endl
        << "===================================================" << endl;
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        log2 << setw(3) << i << " (" << constants::atnr2letter(wavy.get_atom_charge(i)) << ") |" << scientific << setprecision(6) << setw(13) << dipole_moments[i][3] - wavy.get_atom_charge(i) << " | " << scientific << setprecision(6) << setw(14) << dipole_moments[i][0] << ", " << setw(14) << dipole_moments[i][1] << ", " << setw(14) << dipole_moments[i][2] << endl;
    }
    return dipole_moments;
}

void polarizabilities(options &opt, std::ostream &log2)
{
    using namespace std;
    std::vector<WFN> wavy;
    for (int i = 0; i < 7; i++)
    {
        wavy.emplace_back(e_origin::NOT_YET_DEFINED);
        wavy[i].read_known_wavefunction_format(opt.pol_wfns[i], log2, opt.debug);
    }

    // check that all WFN have the same number of atoms

    if (opt.debug)
        log2 << "Starting calculation of Polarizabilities" << endl;

    if (opt.debug)
        log2 << opt.properties.resolution << " " << opt.properties.radius << endl;
    readxyzMinMax_fromWFN(wavy[0], opt.properties, true);
    if (opt.debug)
    {
        log2 << "Resolution: " << opt.properties.resolution << endl;
        log2 << "MinMax:" << endl;
        for (int i = 0; i < 6; i++)
            log2 << setw(14) << scientific << opt.properties.MinMax[i];
        log2 << endl;
        log2 << "Steps:" << endl;
        for (int i = 0; i < 3; i++)
            log2 << setw(14) << scientific << opt.properties.NbSteps[i];
        log2 << endl;
    }
    cube SPHER(opt.properties.NbSteps, wavy[0].get_ncen(), true);

    SPHER.give_parent_wfn(wavy[0]);
    vec stepsizes{ (opt.properties.MinMax[3] - opt.properties.MinMax[0]) / opt.properties.NbSteps[0],
                  (opt.properties.MinMax[4] - opt.properties.MinMax[1]) / opt.properties.NbSteps[1],
                  (opt.properties.MinMax[5] - opt.properties.MinMax[2]) / opt.properties.NbSteps[2] };

    for (int i = 0; i < 3; i++)
    {
        SPHER.set_origin(i, opt.properties.MinMax[i]);
        SPHER.set_vector(i, i, stepsizes[i]);
    }
    if (opt.debug)
        log2 << "Origins etc are set up" << endl;
    SPHER.set_comment1("Calculated Atomic Hirshfeld deformation density values using NoSpherA2");
    SPHER.set_comment2("from" + wavy[0].get_path().string());
    SPHER.set_path((wavy[0].get_path().parent_path() / wavy[0].get_path().stem()).string() + "_spher.cube");

    log2 << "Calculating for " << fixed << setprecision(0) << opt.properties.n_grid_points() << " Gridpoints." << endl;

    log2 << "Calcualting spherical Rho...";
    Calc_Spherical_Dens(SPHER, wavy[0], opt.properties.radius, log2, false);
    log2 << " ...done!" << endl;
    vec3 dipoles(7); // 0, +x, -x, +y, -y, +z, -z
    for (int i = 0; i < 7; i++)
    {
        dipoles[i] = dipole_moments(wavy[i], SPHER, opt.properties, opt.threads, log2, opt.debug);
    }
    vec3 polarizabilities(wavy[0].get_ncen());
    for (int i = 0; i < wavy[0].get_ncen(); i++)
    {
        polarizabilities[i] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
        vec dx = { (dipoles[1][i][0] - dipoles[2][i][0]),
                  (dipoles[3][i][0] - dipoles[4][i][0]),
                  (dipoles[5][i][0] - dipoles[6][i][0]) };
        vec dy = { (dipoles[1][i][1] - dipoles[2][i][1]),
                  (dipoles[3][i][1] - dipoles[4][i][1]),
                  (dipoles[5][i][1] - dipoles[6][i][1]) };
        vec dz = { (dipoles[1][i][2] - dipoles[2][i][2]),
                  (dipoles[3][i][2] - dipoles[4][i][2]),
                  (dipoles[5][i][2] - dipoles[6][i][2]) };
        polarizabilities[i][0][0] = dx[0] / 2 / opt.efield;
        polarizabilities[i][0][1] = dx[1] / 2 / opt.efield;
        polarizabilities[i][0][2] = dx[2] / 2 / opt.efield;
        polarizabilities[i][1][0] = dy[0] / 2 / opt.efield;
        polarizabilities[i][1][1] = dy[1] / 2 / opt.efield;
        polarizabilities[i][1][2] = dy[2] / 2 / opt.efield;
        polarizabilities[i][2][0] = dz[0] / 2 / opt.efield;
        polarizabilities[i][2][1] = dz[1] / 2 / opt.efield;
        polarizabilities[i][2][2] = dz[2] / 2 / opt.efield;
    }
    // print the results per atom
    log2 << "Polarizabilities:\n atom   |    charge    |       xx,            xy,            xz,            yx,            yy,            yz,            zx,            zy,            zz" << endl
        << "========|==============|=======================================================================================================================================" << endl;
    for (int i = 0; i < wavy[0].get_ncen(); i++)
    {
        log2 << setw(3) << i << " (" << constants::atnr2letter(wavy[0].get_atom_charge(i)) << ") |"
            << scientific << setprecision(6) << setw(13) << dipoles[0][i][3] - wavy[0].get_atom_charge(i) << " |"
            << setw(14) << polarizabilities[i][0][0] << ","
            << setw(14) << polarizabilities[i][0][1] << ","
            << setw(14) << polarizabilities[i][0][2] << ","
            << setw(14) << polarizabilities[i][1][0] << ","
            << setw(14) << polarizabilities[i][1][1] << ","
            << setw(14) << polarizabilities[i][1][2] << ","
            << setw(14) << polarizabilities[i][2][0] << ","
            << setw(14) << polarizabilities[i][2][1] << ","
            << setw(14) << polarizabilities[i][2][2] << endl;
    }
    std::cout << "\n\nProperties calculation done!" << std::endl;
}

// end here
