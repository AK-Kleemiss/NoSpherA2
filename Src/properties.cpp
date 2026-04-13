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
    const vector<atom> focus_atom{ wavy.get_atoms()[ignore_atom] };

    evaluate_cube_in_radius_mapped(
        CubeHDEF,
        wrap,
        focus_atom,
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
    const vector<atom> focus_atom{ wavy.get_atoms()[ignore_atom] };

    evaluate_cube_in_radius_mapped(
        Cubes[cube_type::Hirsh],
        wrap,
        focus_atom,
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

    Cubes[cube_type::Rho].evaluate_on_grid(
        [&](const d3 &PosGrid, const i3 &, const i3 &mapped_idx)
        {
            if (!is_within_radius(PosGrid, atoms, radius_bohr))
                return 0.0;

            const PropValues values = compute_prop_values(Cubes, wavy, PosGrid);
            accumulate_prop_values(Cubes, mapped_idx, values);
            return values.rho;
        },
        wrap);

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
