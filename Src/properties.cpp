#include "wfn_class.h"
#include "properties.h"
#include "convenience.h"
#include "spherical_density.h"
#include "cell.h"
#include "cube.h"
#include "constants.h"

void Calc_Spherical_Dens(
    cube &CubeSpher,
    WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    ProgressBar *progress = new ProgressBar(CubeSpher.get_size(0), 50, "=", " ", "Calculating Spherical Density");

    vector<Thakkar> atoms;
    for (int a = 0; a < 92; a++) {
        atoms.push_back(Thakkar(a));
        atoms[a].make_interpolator(1.005*1.005*1.005, 1E-7);
    }

    const int low_i = wrap ? -CubeSpher.get_size(0) : 0;
    const int high_i = wrap ? 2 * CubeSpher.get_size(0) : CubeSpher.get_size(0);
    const int low_j = wrap ? -CubeSpher.get_size(1) : 0;
    const int high_j = wrap ? 2 * CubeSpher.get_size(1) : CubeSpher.get_size(1);
    const int low_k = wrap ? -CubeSpher.get_size(2) : 0;
    const int high_k = wrap ? 2 * CubeSpher.get_size(2) : CubeSpher.get_size(2);

#pragma omp parallel for schedule(dynamic)
    for (int i = low_i; i < high_i; i++)
    {
        vec dists(wavy.get_ncen(), 0);
        for (int j = low_j; j < high_j; j++)
            for (int k = low_k; k < high_k; k++)
            {

                const double PosGrid[3]{i * CubeSpher.get_vector(0, 0) + j * CubeSpher.get_vector(0, 1) + k * CubeSpher.get_vector(0, 2) + CubeSpher.get_origin(0),
                                        i * CubeSpher.get_vector(1, 0) + j * CubeSpher.get_vector(1, 1) + k * CubeSpher.get_vector(1, 2) + CubeSpher.get_origin(1),
                                        i * CubeSpher.get_vector(2, 0) + j * CubeSpher.get_vector(2, 1) + k * CubeSpher.get_vector(2, 2) + CubeSpher.get_origin(2)};

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                {
                    dists[a] = sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(i, 0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(i, 1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(i, 2), 2));
                    if (dists[a] < constants::ang2bohr(radius))
                        skip = false;
                }
                if (skip)
                    continue;

                double dens_all = 0.0;
                for (int a = 0; a < wavy.get_ncen(); a++)
                {
                    dens_all += atoms[wavy.get_atom_charge(a) - 1].get_interpolated_density(dists[a]);
                }

                int temp_i, temp_j, temp_k;
                if (i < 0)
                    temp_i = i + CubeSpher.get_size(0);
                else if (i < CubeSpher.get_size(0))
                    temp_i = i;
                else
                    temp_i = i - CubeSpher.get_size(0);

                if (j < 0)
                    temp_j = j + CubeSpher.get_size(1);
                else if (j < CubeSpher.get_size(1))
                    temp_j = j;
                else
                    temp_j = j - CubeSpher.get_size(1);

                if (k < 0)
                    temp_k = k + CubeSpher.get_size(2);
                else if (k < CubeSpher.get_size(2))
                    temp_k = k;
                else
                    temp_k = k - CubeSpher.get_size(2);

                CubeSpher.set_value(temp_i, temp_j, temp_k, CubeSpher.get_value(temp_i, temp_j, temp_k) + dens_all);
            }
        progress->update();
    }
    delete (progress);

    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
};

void Calc_Static_Def(
    cube &CubeDEF,
    cube &CubeRho,
    WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    ProgressBar *progress = new ProgressBar(CubeDEF.get_size(0), 50, "=", " ", "Calculating Deformation");

    vector<Thakkar> atoms;
    for (int a = 0; a < wavy.get_ncen(); a++)
        atoms.push_back(Thakkar(wavy.get_atom_charge(a)));

    const int low_i = wrap ? -CubeDEF.get_size(0) : 0;
    const int high_i = wrap ? 2 * CubeDEF.get_size(0) : CubeDEF.get_size(0);
    const int low_j = wrap ? -CubeDEF.get_size(1) : 0;
    const int high_j = wrap ? 2 * CubeDEF.get_size(1) : CubeDEF.get_size(1);
    const int low_k = wrap ? -CubeDEF.get_size(2) : 0;
    const int high_k = wrap ? 2 * CubeDEF.get_size(2) : CubeDEF.get_size(2);

#pragma omp parallel for schedule(dynamic)
    for (int i = low_i; i < high_i; i++)
    {
        for (int j = low_j; j < high_j; j++)
            for (int k = low_k; k < high_k; k++)
            {

                const double PosGrid[3]{i * CubeDEF.get_vector(0, 0) + j * CubeDEF.get_vector(0, 1) + k * CubeDEF.get_vector(0, 2) + CubeDEF.get_origin(0),
                                        i * CubeDEF.get_vector(1, 0) + j * CubeDEF.get_vector(1, 1) + k * CubeDEF.get_vector(1, 2) + CubeDEF.get_origin(1),
                                        i * CubeDEF.get_vector(2, 0) + j * CubeDEF.get_vector(2, 1) + k * CubeDEF.get_vector(2, 2) + CubeDEF.get_origin(2)};

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(a,0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(a,1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(a,2), 2)) < constants::ang2bohr(radius))
                        skip = false;
                if (skip)
                    continue;

                double dens_all = 0.0;
                double dist;
                for (int a = 0; a < wavy.get_ncen(); a++)
                {
                    dist = sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(a,0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(a,1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(a,2), 2));
                    dens_all += atoms[a].get_radial_density(dist);
                }

                int temp_i, temp_j, temp_k;
                if (i < 0)
                    temp_i = i + CubeDEF.get_size(0);
                else if (i < CubeDEF.get_size(0))
                    temp_i = i;
                else
                    temp_i = i - CubeDEF.get_size(0);

                if (j < 0)
                    temp_j = j + CubeDEF.get_size(1);
                else if (j < CubeDEF.get_size(1))
                    temp_j = j;
                else
                    temp_j = j - CubeDEF.get_size(1);

                if (k < 0)
                    temp_k = k + CubeDEF.get_size(2);
                else if (k < CubeDEF.get_size(2))
                    temp_k = k;
                else
                    temp_k = k - CubeDEF.get_size(2);

                dens_all -= CubeRho.get_value(temp_i, temp_j, temp_k);
                CubeDEF.set_value(temp_i, temp_j, temp_k, CubeDEF.get_value(temp_i, temp_j, temp_k) - dens_all);
            }
        progress->update();
    }
    delete (progress);

    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
};

void Calc_Static_Def(
    cube &CubeDEF,
    cube &CubeRho,
    cube &CubeSpher,
    WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    ProgressBar *progress = new ProgressBar(CubeDEF.get_size(0), 50, "=", " ", "Calculating Deformation");

    const int low_i = wrap ? -CubeDEF.get_size(0) : 0;
    const int high_i = wrap ? 2 * CubeDEF.get_size(0) : CubeDEF.get_size(0);
    const int low_j = wrap ? -CubeDEF.get_size(1) : 0;
    const int high_j = wrap ? 2 * CubeDEF.get_size(1) : CubeDEF.get_size(1);
    const int low_k = wrap ? -CubeDEF.get_size(2) : 0;
    const int high_k = wrap ? 2 * CubeDEF.get_size(2) : CubeDEF.get_size(2);

#pragma omp parallel for schedule(dynamic)
    for (int i = low_i; i < high_i; i++)
    {
        for (int j = low_j; j < high_j; j++)
            for (int k = low_k; k < high_k; k++)
            {

                const double PosGrid[3]{i * CubeDEF.get_vector(0, 0) + j * CubeDEF.get_vector(0, 1) + k * CubeDEF.get_vector(0, 2) + CubeDEF.get_origin(0),
                                        i * CubeDEF.get_vector(1, 0) + j * CubeDEF.get_vector(1, 1) + k * CubeDEF.get_vector(1, 2) + CubeDEF.get_origin(1),
                                        i * CubeDEF.get_vector(2, 0) + j * CubeDEF.get_vector(2, 1) + k * CubeDEF.get_vector(2, 2) + CubeDEF.get_origin(2)};

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(a,0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(a,1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(a,2), 2)) < constants::ang2bohr(radius))
                        skip = false;
                if (skip)
                    continue;

                int temp_i, temp_j, temp_k;
                if (i < 0)
                    temp_i = i + CubeDEF.get_size(0);
                else if (i < CubeDEF.get_size(0))
                    temp_i = i;
                else
                    temp_i = i - CubeDEF.get_size(0);

                if (j < 0)
                    temp_j = j + CubeDEF.get_size(1);
                else if (j < CubeDEF.get_size(1))
                    temp_j = j;
                else
                    temp_j = j - CubeDEF.get_size(1);

                if (k < 0)
                    temp_k = k + CubeDEF.get_size(2);
                else if (k < CubeDEF.get_size(2))
                    temp_k = k;
                else
                    temp_k = k - CubeDEF.get_size(2);

                double rho = CubeRho.get_value(temp_i, temp_j, temp_k);
                double spher = CubeSpher.get_value(temp_i, temp_j, temp_k);
                double temp = rho - spher;
                CubeDEF.set_value(temp_i, temp_j, temp_k, CubeDEF.get_value(temp_i, temp_j, temp_k) + temp);
            }
        progress->update();
    }
    delete (progress);

    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
};

void Calc_Hirshfeld(
    cube &CubeHDEF,
    cube &CubeRho,
    WFN &wavy,
    double radius,
    int ignore_atom,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    ProgressBar *progress = new ProgressBar(CubeHDEF.get_size(0), 50, "=", " ", "Calculating Values");

    vector<Thakkar> atoms;
    for (int a = 0; a < wavy.get_ncen(); a++)
        atoms.push_back(Thakkar(wavy.get_atom_charge(a)));

    const int low_i = wrap ? -CubeHDEF.get_size(0) : 0;
    const int high_i = wrap ? 2 * CubeHDEF.get_size(0) : CubeHDEF.get_size(0);
    const int low_j = wrap ? -CubeHDEF.get_size(1) : 0;
    const int high_j = wrap ? 2 * CubeHDEF.get_size(1) : CubeHDEF.get_size(1);
    const int low_k = wrap ? -CubeHDEF.get_size(2) : 0;
    const int high_k = wrap ? 2 * CubeHDEF.get_size(2) : CubeHDEF.get_size(2);

#pragma omp parallel for schedule(dynamic)
    for (int i = low_i; i < high_i; i++)
    {
        for (int j = low_j; j < high_j; j++)
            for (int k = low_k; k < high_k; k++)
            {

                std::array<double, 3> PosGrid = CubeRho.get_pos(i, j, k);

                bool skip = true;
                double dist = sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(ignore_atom, 0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(ignore_atom, 1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(ignore_atom, 2), 2));
                if (dist < constants::ang2bohr(radius))
                    skip = false;
                if (skip)
                    continue;

                double dens_choice = 0.0;
                double dens_all = 0.0;
                double temp;
                for (int a = 0; a < wavy.get_ncen(); a++)
                {
                    dist = sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(a,0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(a,1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(a,2), 2));
                    temp = atoms[a].get_radial_density(dist);
                    if (ignore_atom == a)
                        dens_choice = temp;
                    dens_all += temp;
                }

                int temp_i, temp_j, temp_k;
                if (i < 0)
                    temp_i = i + CubeHDEF.get_size(0);
                else if (i < CubeHDEF.get_size(0))
                    temp_i = i;
                else
                    temp_i = i - CubeHDEF.get_size(0);

                if (j < 0)
                    temp_j = j + CubeHDEF.get_size(1);
                else if (j < CubeHDEF.get_size(1))
                    temp_j = j;
                else
                    temp_j = j - CubeHDEF.get_size(1);

                if (k < 0)
                    temp_k = k + CubeHDEF.get_size(2);
                else if (k < CubeHDEF.get_size(2))
                    temp_k = k;
                else
                    temp_k = k - CubeHDEF.get_size(2);

                dens_all = dens_choice / dens_all * CubeRho.get_value(temp_i, temp_j, temp_k);
                CubeHDEF.set_value(temp_i, temp_j, temp_k, CubeHDEF.get_value(temp_i, temp_j, temp_k) + dens_all - dens_choice);
            }
        progress->update();
    }
    delete (progress);

    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
};

void Calc_Hirshfeld(
    cube &CubeHDEF,
    cube &CubeRho,
    cube &CubeSpherical,
    WFN &wavy,
    double radius,
    int ignore_atom,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    ProgressBar *progress = new ProgressBar(CubeHDEF.get_size(0), 50, "=", " ", "Calculating Values");
    Thakkar atom(wavy.get_atom_charge(ignore_atom));

    const int low_i = wrap ? -CubeHDEF.get_size(0) : 0;
    const int high_i = wrap ? 2 * CubeHDEF.get_size(0) : CubeHDEF.get_size(0);
    const int low_j = wrap ? -CubeHDEF.get_size(1) : 0;
    const int high_j = wrap ? 2 * CubeHDEF.get_size(1) : CubeHDEF.get_size(1);
    const int low_k = wrap ? -CubeHDEF.get_size(2) : 0;
    const int high_k = wrap ? 2 * CubeHDEF.get_size(2) : CubeHDEF.get_size(2);

#pragma omp parallel for schedule(dynamic)
    for (int i = low_i; i < high_i; i++)
    {
        for (int j = low_j; j < high_j; j++)
            for (int k = low_k; k < high_k; k++)
            {

                std::array<double, 3> PosGrid = CubeRho.get_pos(i, j, k);

                bool skip = true;
                double dist = sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(ignore_atom,0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(ignore_atom, 1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(ignore_atom, 2), 2));
                if (dist < constants::ang2bohr(radius))
                    skip = false;
                if (skip)
                    continue;

                double dens_choice = atom.get_radial_density(dist);
                int temp_i, temp_j, temp_k;
                if (i < 0)
                    temp_i = i + CubeHDEF.get_size(0);
                else if (i < CubeHDEF.get_size(0))
                    temp_i = i;
                else
                    temp_i = i - CubeHDEF.get_size(0);

                if (j < 0)
                    temp_j = j + CubeHDEF.get_size(1);
                else if (j < CubeHDEF.get_size(1))
                    temp_j = j;
                else
                    temp_j = j - CubeHDEF.get_size(1);

                if (k < 0)
                    temp_k = k + CubeHDEF.get_size(2);
                else if (k < CubeHDEF.get_size(2))
                    temp_k = k;
                else
                    temp_k = k - CubeHDEF.get_size(2);
                CubeHDEF.set_value(temp_i, temp_j, temp_k, CubeHDEF.get_value(temp_i, temp_j, temp_k) + (dens_choice / CubeSpherical.get_value(temp_i, temp_j, temp_k) * CubeRho.get_value(temp_i, temp_j, temp_k)) - dens_choice);
            }
        progress->update();
    }
    delete (progress);

    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
};

void Calc_Hirshfeld_atom(
    cube &CubeHirsh,
    cube &CubeRho,
    cube &CubeSpherical,
    WFN &wavy,
    double radius,
    int ignore_atom,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    ProgressBar *progress = new ProgressBar(CubeHirsh.get_size(0), 50, "=", " ", "Calculating Values");
    Thakkar atom(wavy.get_atom_charge(ignore_atom));

    const int low_i = wrap ? -CubeHirsh.get_size(0) : 0;
    const int high_i = wrap ? 2 * CubeHirsh.get_size(0) : CubeHirsh.get_size(0);
    const int low_j = wrap ? -CubeHirsh.get_size(1) : 0;
    const int high_j = wrap ? 2 * CubeHirsh.get_size(1) : CubeHirsh.get_size(1);
    const int low_k = wrap ? -CubeHirsh.get_size(2) : 0;
    const int high_k = wrap ? 2 * CubeHirsh.get_size(2) : CubeHirsh.get_size(2);

#pragma omp parallel for schedule(dynamic)
    for (int i = low_i; i < high_i; i++)
    {
        for (int j = low_j; j < high_j; j++)
            for (int k = low_k; k < high_k; k++)
            {

                std::array<double, 3> PosGrid = CubeRho.get_pos(i, j, k);

                bool skip = true;
                double dist = sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(ignore_atom, 0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(ignore_atom, 1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(ignore_atom, 2), 2));
                if (dist < constants::ang2bohr(radius))
                    skip = false;
                if (skip)
                    continue;

                double dens_choice = atom.get_radial_density(dist);
                int temp_i, temp_j, temp_k;
                if (i < 0)
                    temp_i = i + CubeHirsh.get_size(0);
                else if (i < CubeHirsh.get_size(0))
                    temp_i = i;
                else
                    temp_i = i - CubeHirsh.get_size(0);

                if (j < 0)
                    temp_j = j + CubeHirsh.get_size(1);
                else if (j < CubeHirsh.get_size(1))
                    temp_j = j;
                else
                    temp_j = j - CubeHirsh.get_size(1);

                if (k < 0)
                    temp_k = k + CubeHirsh.get_size(2);
                else if (k < CubeHirsh.get_size(2))
                    temp_k = k;
                else
                    temp_k = k - CubeHirsh.get_size(2);
                CubeHirsh.set_value(temp_i, temp_j, temp_k, CubeHirsh.get_value(temp_i, temp_j, temp_k) + (dens_choice / CubeSpherical.get_value(temp_i, temp_j, temp_k) * CubeRho.get_value(temp_i, temp_j, temp_k)));
            }
        progress->update();
    }
    delete (progress);

    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
};

void Calc_Rho(
    cube &CubeRho,
    WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    ProgressBar *progress = new ProgressBar(CubeRho.get_size(0), 50, "=", " ", "Calculating Values");

    const int low_i = wrap ? -CubeRho.get_size(0) : 0;
    const int high_i = wrap ? 2 * CubeRho.get_size(0) : CubeRho.get_size(0);
    const int low_j = wrap ? -CubeRho.get_size(1) : 0;
    const int high_j = wrap ? 2 * CubeRho.get_size(1) : CubeRho.get_size(1);
    const int low_k = wrap ? -CubeRho.get_size(2) : 0;
    const int high_k = wrap ? 2 * CubeRho.get_size(2) : CubeRho.get_size(2);

#pragma omp parallel for schedule(dynamic)
    for (int i = low_i; i < high_i; i++)
    {
        for (int j = low_j; j < high_j; j++)
            for (int k = low_k; k < high_k; k++)
            {

                std::array<double, 3> PosGrid = CubeRho.get_pos(i, j, k);
                double Rho = 0;

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(a,0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(a,1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(a,2), 2)) < constants::ang2bohr(radius))
                        skip = false;
                if (skip)
                    continue;

                Rho = wavy.compute_dens(PosGrid[0], PosGrid[1], PosGrid[2]);

                int temp_i, temp_j, temp_k;
                if (i < 0)
                    temp_i = i + CubeRho.get_size(0);
                else if (i < CubeRho.get_size(0))
                    temp_i = i;
                else
                    temp_i = i - CubeRho.get_size(0);

                if (j < 0)
                    temp_j = j + CubeRho.get_size(1);
                else if (j < CubeRho.get_size(1))
                    temp_j = j;
                else
                    temp_j = j - CubeRho.get_size(1);

                if (k < 0)
                    temp_k = k + CubeRho.get_size(2);
                else if (k < CubeRho.get_size(2))
                    temp_k = k;
                else
                    temp_k = k - CubeRho.get_size(2);

                CubeRho.set_value(temp_i, temp_j, temp_k, CubeRho.get_value(temp_i, temp_j, temp_k) + Rho);
            }
        progress->update();
    }
    delete (progress);

    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
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
        vec2 d(5);
        for (int i = 0; i < 5; i++)
            d[i].resize(wavy.get_ncen(), 0.0);
        const int n = wavy.get_nmo(true);
        vec phi(n, 0.0);
        // #pragma omp for schedule(dynamic)
        for (int i = 0; i < CubeRho.get_size(0); i++)
        {
            for (int j = 0; j < CubeRho.get_size(1); j++)
                for (int k = 0; k < CubeRho.get_size(2); k++)
                {

                    std::array<double, 3> PosGrid = CubeRho.get_pos(i, j, k);
                    CubeRho.set_value(i, j, k, wavy.compute_dens(PosGrid[0], PosGrid[1], PosGrid[2], d, phi));
                }
            progress->update();
        }
    }
    delete (progress);

    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
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
            {

                std::array<double, 3> PosGrid = CubeMO.get_pos(i, j, k);
                CubeMO.set_value(i, j, k, wavy.compute_MO_spherical(PosGrid[0], PosGrid[1], PosGrid[2], MO));
            }
        if (!nodate)
            progress->update();
    }
    if (!nodate)
    {
        delete (progress);

        _time_point end = get_time();
        if (get_sec(start, end) < 60)
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
        else if (get_sec(start, end) < 3600)
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
        else
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
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
            {

                std::array<double, 3> PosGrid = Cube_S_Rho.get_pos(i, j, k);
                Cube_S_Rho.set_value(i, j, k, wavy.compute_spin_dens(PosGrid[0], PosGrid[1], PosGrid[2], d, phi));
            }
        if (!nodate)
            progress->update();
    }
    if (!nodate)
    {
        delete (progress);

        _time_point end = get_time();
        if (get_sec(start, end) < 60)
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
        else if (get_sec(start, end) < 3600)
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
        else
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
    }
};

void Calc_Prop(
    cube &CubeRho,
    cube &CubeRDG,
    cube &CubeElf,
    cube &CubeEli,
    cube &CubeLap,
    cube &CubeESP,
    WFN &wavy,
    double radius,
    std::ostream &file,
    bool test,
    bool wrap)
{
    using namespace std;
    _time_point start = get_time();
    ProgressBar *progress = NULL;
    if (!test)
        progress = new ProgressBar(CubeRho.get_size(0), 50, "=", " ", "Calculating Values");

    const int low_i = wrap ? -CubeRho.get_size(0) : 0;
    const int high_i = wrap ? 2 * CubeRho.get_size(0) : CubeRho.get_size(0);
    const int low_j = wrap ? -CubeRho.get_size(1) : 0;
    const int high_j = wrap ? 2 * CubeRho.get_size(1) : CubeRho.get_size(1);
    const int low_k = wrap ? -CubeRho.get_size(2) : 0;
    const int high_k = wrap ? 2 * CubeRho.get_size(2) : CubeRho.get_size(2);

#pragma omp parallel for schedule(dynamic)
    for (int i = low_i; i < high_i; i++)
    {
        for (int j = low_j; j < high_j; j++)
            for (int k = low_k; k < high_k; k++)
            {

                std::array<double, 3> PosGrid = CubeRho.get_pos(i, j, k);
                double Rho = 0,
                    Grad = 0,
                    Elf = 0,
                    Eli = 0,
                    Lap = 0,
                    Hess[9]{0, 0, 0, 0, 0, 0, 0, 0, 0};

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(a,0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(a,1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(a,2), 2)) < constants::ang2bohr(radius))
                        skip = false;
                if (skip)
                    continue;

                if (CubeESP.get_loaded() && !CubeRDG.get_loaded())
                    Rho = wavy.compute_dens(
                        PosGrid[0], PosGrid[1], PosGrid[2]);

                if (CubeRDG.get_loaded() && CubeLap.get_loaded() && (CubeElf.get_loaded() || CubeEli.get_loaded()))
                    wavy.computeValues(
                        PosGrid,
                        Rho,
                        Grad,
                        Hess,
                        Elf,
                        Eli,
                        Lap);
                else if ((CubeElf.get_loaded() && CubeEli.get_loaded()) && !CubeRDG.get_loaded() && !CubeLap.get_loaded())
                    wavy.computeELIELF(
                        PosGrid,
                        Elf,
                        Eli);
                else if (CubeElf.get_loaded() && !CubeEli.get_loaded() && !CubeRDG.get_loaded() && !CubeLap.get_loaded())
                    wavy.computeELF(
                        PosGrid,
                        Elf);
                else if (!CubeElf.get_loaded() && CubeEli.get_loaded() && !CubeRDG.get_loaded() && !CubeLap.get_loaded())
                    wavy.computeELI(
                        PosGrid,
                        Eli);
                else if (CubeElf.get_loaded() && CubeEli.get_loaded() && CubeLap.get_loaded() && !CubeRDG.get_loaded())
                    wavy.computeLapELIELF(
                        PosGrid,
                        Elf,
                        Eli,
                        Lap);
                else if (!CubeElf.get_loaded() && CubeEli.get_loaded() && CubeLap.get_loaded() && !CubeRDG.get_loaded())
                    wavy.computeLapELI(
                        PosGrid,
                        Eli,
                        Lap);
                else
                    wavy.computeValues(
                        PosGrid,
                        Rho,
                        Grad,
                        Hess,
                        Elf,
                        Eli,
                        Lap);

                if (CubeRDG.get_loaded())
                    Rho = get_lambda_1(Hess) < 0 ? -Rho : Rho;

                int temp_i, temp_j, temp_k;
                if (i < 0)
                    temp_i = i + CubeRho.get_size(0);
                else if (i < CubeRho.get_size(0))
                    temp_i = i;
                else
                    temp_i = i - CubeRho.get_size(0);

                if (j < 0)
                    temp_j = j + CubeRho.get_size(1);
                else if (j < CubeRho.get_size(1))
                    temp_j = j;
                else
                    temp_j = j - CubeRho.get_size(1);

                if (k < 0)
                    temp_k = k + CubeRho.get_size(2);
                else if (k < CubeRho.get_size(2))
                    temp_k = k;
                else
                    temp_k = k - CubeRho.get_size(2);

                CubeRho.set_value(temp_i, temp_j, temp_k, CubeRho.get_value(temp_i, temp_j, temp_k) + Rho);
                if (CubeRDG.get_loaded())
                {
                    if (isnan(Grad))
                        Grad = 0;
                    if (isinf(Grad))
                        Grad = 0;
                    CubeRDG.set_value(temp_i, temp_j, temp_k, CubeRDG.get_value(temp_i, temp_j, temp_k) + Grad);
                }
                if (CubeLap.get_loaded())
                {
                    if (isnan(Lap))
                        Lap = 0;
                    if (isinf(Lap))
                        Lap = 0;
                    CubeLap.set_value(temp_i, temp_j, temp_k, CubeLap.get_value(temp_i, temp_j, temp_k) + Lap);
                }
                if (CubeElf.get_loaded())
                {
                    if (isnan(Elf))
                        Elf = 0;
                    if (isinf(Elf))
                        Elf = 0;
                    CubeElf.set_value(temp_i, temp_j, temp_k, CubeElf.get_value(temp_i, temp_j, temp_k) + Elf);
                }
                if (CubeEli.get_loaded())
                {
                    if (isnan(Eli))
                        Eli = 0;
                    if (isinf(Eli))
                        Eli = 0;
                    CubeEli.set_value(temp_i, temp_j, temp_k, CubeEli.get_value(temp_i, temp_j, temp_k) + Eli);
                }
            }
        if (!test)
            progress->update();
    }
    if (!test)
    {
        delete (progress);
        _time_point end = get_time();
        if (get_sec(start, end) < 60)
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
        else if (get_sec(start, end) < 3600)
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
        else
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
    }
};

void Calc_ESP(
    cube &CubeESP,
    WFN &wavy,
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
        d2[i].resize(wavy.get_ncen());
        for (int j = 0; j < wavy.get_ncen(); j++)
        {
            if (i == j)
            {
                d2[i][j] = 0;
                continue;
            }
            d2[i][j] = pow(wavy.get_atom_coordinate(i,0) - wavy.get_atom_coordinate(j,0), 2) + pow(wavy.get_atom_coordinate(i,1) - wavy.get_atom_coordinate(j,1), 2) + pow(wavy.get_atom_coordinate(i,2) - wavy.get_atom_coordinate(j,2), 2);
        }
    }

    ProgressBar *progress = NULL;
    if (!no_date)
        progress = new ProgressBar(CubeESP.get_size(0), 50, "=", " ", "Calculating ESP");

    const int low_i = wrap ? -CubeESP.get_size(0) : 0;
    const int high_i = wrap ? 2 * CubeESP.get_size(0) : CubeESP.get_size(0);
    const int low_j = wrap ? -CubeESP.get_size(1) : 0;
    const int high_j = wrap ? 2 * CubeESP.get_size(1) : CubeESP.get_size(1);
    const int low_k = wrap ? -CubeESP.get_size(2) : 0;
    const int high_k = wrap ? 2 * CubeESP.get_size(2) : CubeESP.get_size(2);
    double temp;
    int temp_i, temp_j, temp_k;

#pragma omp parallel for schedule(dynamic) private(temp, temp_i, temp_j, temp_k)
    for (int i = low_i; i < high_i; i++)
    {
        for (int j = low_j; j < high_j; j++)
            for (int k = low_k; k < high_k; k++)
            {
                std::array<double, 3> PosGrid = CubeESP.get_pos(i, j, k);

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(a,0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(a,1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(a,2), 2)) < constants::ang2bohr(radius))
                        skip = false;
                if (skip)
                    continue;

                if (i < 0)
                    temp_i = i + CubeESP.get_size(0);
                else if (i < CubeESP.get_size(0))
                    temp_i = i;
                else
                    temp_i = i - CubeESP.get_size(0);

                if (j < 0)
                    temp_j = j + CubeESP.get_size(1);
                else if (j < CubeESP.get_size(1))
                    temp_j = j;
                else
                    temp_j = j - CubeESP.get_size(1);

                if (k < 0)
                    temp_k = k + CubeESP.get_size(2);
                else if (k < CubeESP.get_size(2))
                    temp_k = k;
                else
                    temp_k = k - CubeESP.get_size(2);

                temp = wavy.computeESP(PosGrid, d2);
                CubeESP.set_value(temp_i, temp_j, temp_k, CubeESP.get_value(temp_i, temp_j, temp_k) + temp);
                // CubeESP.set_value(i, j, k, computeESP(PosGrid, d2, wavy));
            }
        if (!no_date)
            progress->update();
    }
    if (!no_date)
    {
        delete (progress);

        _time_point end = get_time();
        if (get_sec(start, end) < 60)
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
        else if (get_sec(start, end) < 3600)
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
        else
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
    }
};

void Calc_MO(
    cube &CubeMO,
    int mo,
    WFN &wavy,
    double radius,
    std::ostream &file,
    bool wrap)
{
    using namespace std;
    err_checkf(mo <= wavy.get_nmo(), to_string(mo) + " bigger MO selected than " + to_string(wavy.get_nmo()) + " contained in the wavefunctions!", file);
    _time_point start = get_time();

    ProgressBar *progress = new ProgressBar(CubeMO.get_size(0), 50, "=", " ", "Calculating MO");

    const int low_i = wrap ? -CubeMO.get_size(0) : 0;
    const int high_i = wrap ? 2 * CubeMO.get_size(0) : CubeMO.get_size(0);
    const int low_j = wrap ? -CubeMO.get_size(1) : 0;
    const int high_j = wrap ? 2 * CubeMO.get_size(1) : CubeMO.get_size(1);
    const int low_k = wrap ? -CubeMO.get_size(2) : 0;
    const int high_k = wrap ? 2 * CubeMO.get_size(2) : CubeMO.get_size(2);

#pragma omp parallel for schedule(dynamic)
    for (int i = low_i; i < high_i; i++)
    {
        for (int j = low_j; j < high_j; j++)
            for (int k = low_k; k < high_k; k++)
            {

                std::array<double, 3> PosGrid = CubeMO.get_pos(i, j, k);
                double MO = 0;

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(a,0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(a,1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(a,2), 2)) < constants::ang2bohr(radius))
                        skip = false;
                if (skip)
                    continue;

                MO = wavy.computeMO(
                    PosGrid,
                    mo);

                int temp_i, temp_j, temp_k;
                if (i < 0)
                    temp_i = i + CubeMO.get_size(0);
                else if (i < CubeMO.get_size(0))
                    temp_i = i;
                else
                    temp_i = i - CubeMO.get_size(0);

                if (j < 0)
                    temp_j = j + CubeMO.get_size(1);
                else if (j < CubeMO.get_size(1))
                    temp_j = j;
                else
                    temp_j = j - CubeMO.get_size(1);

                if (k < 0)
                    temp_k = k + CubeMO.get_size(2);
                else if (k < CubeMO.get_size(2))
                    temp_k = k;
                else
                    temp_k = k - CubeMO.get_size(2);

                CubeMO.set_value(temp_i, temp_j, temp_k, CubeMO.get_value(temp_i, temp_j, temp_k) + MO);
            }
        progress->update();
    }
    delete (progress);

    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
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
    if (opt.all_mos)
        for (int mo = 0; mo < wavy.get_nmo(); mo++)
            opt.MOs.push_back(mo);
    if (opt.debug)
        log2 << "Size of MOs: " << opt.MOs.size() << endl;

    vec2 cell_matrix;
    cell_matrix.resize(3);
    for (int i = 0; i < 3; i++)
        cell_matrix[i].resize(3, 0.0);
    if (opt.debug)
    {
        log2 << "cif|resolution|res(bohr)|radius|rad(bohr): " << opt.cif << "|" << opt.resolution << "|" << constants::ang2bohr(opt.resolution) << "|" << opt.radius << "|" << constants::ang2bohr(opt.radius) << endl;
        for (int a = 0; a < wavy.get_ncen(); a++)
            log2 << "Atom " << a << " at " << wavy.get_atom_coordinate(a,0) << " " << wavy.get_atom_coordinate(a,1) << " " << wavy.get_atom_coordinate(a,2) << endl;
    }
    if (opt.cif != "")
        readxyzMinMax_fromCIF(opt.cif, opt.MinMax, opt.NbSteps, cell_matrix, opt.resolution);
    else
    {
        readxyzMinMax_fromWFN(wavy, opt.MinMax, opt.NbSteps, opt.radius, opt.resolution, true);
        for (int i = 0; i < 3; i++)
            cell_matrix[i][i] = constants::ang2bohr(opt.resolution);
    }
    if (opt.debug)
    {
        log2 << "MinMax: ";
        for (int i = 0; i < 6; i++)
            log2 << setw(14) << scientific << opt.MinMax[i];
        log2 << endl;
        log2 << "Steps: ";
        for (int i = 0; i < 3; i++)
            log2 << setw(14) << scientific << opt.NbSteps[i];
        log2 << endl;
        log2 << "Cell Matrix:" << endl;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
                log2 << setw(14) << scientific << cell_matrix[i][j];
            log2 << endl;
        }
    }
    cube Rho(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), true);
    cube RDG(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), opt.rdg);
    cube Elf(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), opt.elf);
    cube Eli(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), opt.eli);
    cube Lap(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), opt.lap);
    cube ESP(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), opt.esp);
    cube MO(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), true);
    cube HDEF(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), opt.hdef);
    cube DEF(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), opt.def);
    cube Hirsh(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), opt.hirsh);
    cube S_Rho(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), opt.hirsh);

    Rho.give_parent_wfn(wavy);
    RDG.give_parent_wfn(wavy);
    Elf.give_parent_wfn(wavy);
    Eli.give_parent_wfn(wavy);
    Lap.give_parent_wfn(wavy);
    ESP.give_parent_wfn(wavy);
    MO.give_parent_wfn(wavy);
    HDEF.give_parent_wfn(wavy);
    DEF.give_parent_wfn(wavy);
    Hirsh.give_parent_wfn(wavy);
    S_Rho.give_parent_wfn(wavy);

    for (int i = 0; i < 3; i++)
    {
        Rho.set_origin(i, opt.MinMax[i]);
        RDG.set_origin(i, opt.MinMax[i]);
        Elf.set_origin(i, opt.MinMax[i]);
        Eli.set_origin(i, opt.MinMax[i]);
        Lap.set_origin(i, opt.MinMax[i]);
        ESP.set_origin(i, opt.MinMax[i]);
        MO.set_origin(i, opt.MinMax[i]);
        HDEF.set_origin(i, opt.MinMax[i]);
        DEF.set_origin(i, opt.MinMax[i]);
        Hirsh.set_origin(i, opt.MinMax[i]);
        S_Rho.set_origin(i, opt.MinMax[i]);
        for (int j = 0; j < 3; j++)
        {
            Rho.set_vector(i, j, cell_matrix[i][j]);
            RDG.set_vector(i, j, cell_matrix[i][j]);
            Elf.set_vector(i, j, cell_matrix[i][j]);
            Eli.set_vector(i, j, cell_matrix[i][j]);
            Lap.set_vector(i, j, cell_matrix[i][j]);
            ESP.set_vector(i, j, cell_matrix[i][j]);
            MO.set_vector(i, j, cell_matrix[i][j]);
            HDEF.set_vector(i, j, cell_matrix[i][j]);
            DEF.set_vector(i, j, cell_matrix[i][j]);
            Hirsh.set_vector(i, j, cell_matrix[i][j]);
            S_Rho.set_vector(i, j, cell_matrix[i][j]);
        }
    }
    if (opt.debug)
        log2 << "Origins etc. are set up" << endl;
    Rho.set_comment1("Calculated density using NoSpherA2");
    RDG.set_comment1("Calculated reduced density gradient using NoSpherA2");
    Elf.set_comment1("Calculated electron localization function using NoSpherA2");
    Eli.set_comment1("Calculated same-spin electron localizability indicator using NoSpherA2");
    Lap.set_comment1("Calculated laplacian of electron density using NoSpherA2");
    ESP.set_comment1("Calculated electrostatic potential using NoSpherA2");
    MO.set_comment1("Calcualted MO values using NoSpherA2");
    HDEF.set_comment1("Calculated Atomic Hirshfeld deformation density values using NoSpherA2");
    DEF.set_comment1("Calculated static deformation density values using NoSpherA2");
    Hirsh.set_comment1("Calculated Hirshfeld atom density values using NoSpherA2");
    S_Rho.set_comment1("Calculated spin density using NoSpherA2");
    Rho.set_comment2("from " + wavy.get_path().string());
    RDG.set_comment2("from " + wavy.get_path().string());
    Elf.set_comment2("from " + wavy.get_path().string());
    Eli.set_comment2("from " + wavy.get_path().string());
    Lap.set_comment2("from " + wavy.get_path().string());
    ESP.set_comment2("from " + wavy.get_path().string());
    MO.set_comment2("from" + wavy.get_path().string());
    HDEF.set_comment2("from" + wavy.get_path().string());
    DEF.set_comment2("from" + wavy.get_path().string());
    Hirsh.set_comment2("from" + wavy.get_path().string());
    S_Rho.set_comment2("from" + wavy.get_path().string());
    Rho.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rho.cube");
    RDG.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rdg.cube");
    Elf.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_elf.cube");
    Eli.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_eli.cube");
    Lap.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_lap.cube");
    ESP.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_esp.cube");
    DEF.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_def.cube");
    Hirsh.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_hirsh.cube");
    S_Rho.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_s_rho.cube");

    if (opt.debug)
    {
        log2 << "Status: " << opt.hdef << opt.def << opt.hirsh << opt.lap << opt.eli << opt.elf << opt.rdg << opt.esp << endl;
        log2 << "Everything is set up; starting calculation..." << endl;
    }
    else
    {
        log2 << "\nCalculating:" << endl;
        if (opt.hdef || opt.def || opt.hirsh)
            log2 << "Rho, ";
        if (opt.hdef || opt.hirsh)
            log2 << "Spherical Rho, ";
        if (opt.def)
            log2 << "Static deformation density, ";
        if (opt.hdef)
            log2 << "Hirshfeld deformation density, ";
        if (opt.hirsh)
            log2 << "Hirshfeld density, ";
        if (opt.lap)
            log2 << "Laplacian, ";
        if (opt.eli)
            log2 << "ELI, ";
        if (opt.elf)
            log2 << "ELF, ";
        if (opt.rdg)
            log2 << "RDG, ";
        if (opt.esp)
            log2 << "ESP, ";
        if (opt.MOs.size() != 0)
            log2 << "MOs, ";
        if (opt.s_rho)
            log2 << "Spin density, ";
        log2 << endl;
    }

    log2 << "Calculating for " << fixed << setprecision(0) << opt.NbSteps[0] * opt.NbSteps[1] * opt.NbSteps[2] << " Gridpoints." << endl;

    if (opt.MOs.size() != 0)
        for (int i = 0; i < opt.MOs.size(); i++)
        {
            log2 << "Calcualting MO: " << opt.MOs[i] << endl;
            MO.set_zero();
            MO.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_MO_" + to_string(opt.MOs[i]) + ".cube");
            Calc_MO(MO, opt.MOs[i], wavy, opt.radius, log2, opt.cif != "");
            MO.write_file(true);
        }

    wavy.delete_unoccupied_MOs();

    if (opt.hdef || opt.def || opt.hirsh)
    {
        log2 << "Calcualting Rho...";
        Calc_Rho(Rho, wavy, opt.radius, log2, opt.cif != "");
        log2 << " ...done!" << endl;
        cube temp(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), opt.hdef || opt.hirsh);
        for (int i = 0; i < 3; i++)
        {
            temp.set_origin(i, opt.MinMax[i]);
            for (int j = 0; j < 3; j++)
                temp.set_vector(i, j, cell_matrix[i][j]);
        }
        if (opt.hdef || opt.hirsh)
        {
            log2 << "Calcualting spherical Rho...";
            Calc_Spherical_Dens(temp, wavy, opt.radius, log2, opt.cif != "");
            log2 << " ...done!" << endl;
        }

        if (opt.def)
        {
            log2 << "Calculating static deformation density...";
            if (opt.hdef)
                Calc_Static_Def(DEF, Rho, temp, wavy, opt.radius, log2, opt.cif != "");
            else
                Calc_Static_Def(DEF, Rho, wavy, opt.radius, log2, opt.cif != "");
            log2 << " ...done!" << endl;
        }

        if (opt.hdef)
        {
            for (int a = 0; a < wavy.get_ncen(); a++)
            {
                log2 << "Calcualting Hirshfeld deformation density for atom: " << a << endl;
                HDEF.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_HDEF_" + to_string(a) + ".cube");
                Calc_Hirshfeld(HDEF, Rho, temp, wavy, opt.radius, a, log2, opt.cif != "");
                HDEF.write_file(true);
                HDEF.set_zero();
            }
        }

        if (opt.hirsh)
        {
            log2 << "Calcualting Hirshfeld density for atom: " << opt.hirsh_number << endl;
            Calc_Hirshfeld_atom(Hirsh, Rho, temp, wavy, opt.radius, opt.hirsh_number, log2, opt.cif != "");
            log2 << "..done!" << endl;
        }
    }

    if (opt.lap || opt.eli || opt.elf || opt.rdg || opt.esp)
        Calc_Prop(Rho, RDG, Elf, Eli, Lap, ESP, wavy, opt.radius, log2, opt.no_date, opt.cif != "");

    if (opt.s_rho)
        Calc_S_Rho(S_Rho, wavy, log2, opt.no_date);

    log2 << "Writing cubes to Disk..." << flush;
    if (opt.rdg)
    {
        Rho.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_signed_rho.cube");
        Rho.write_file(true);
        Rho.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rho.cube");
        Rho.write_file(true, true);
    }
    else if (opt.lap || opt.eli || opt.elf || opt.esp)
        Rho.write_file(true);
    if (opt.rdg)
        RDG.write_file(true);
    if (opt.lap)
        Lap.write_file(true);
    if (opt.elf)
        Elf.write_file(true);
    if (opt.eli)
        Eli.write_file(true);
    if (opt.s_rho)
        S_Rho.write_file(true);
    if (opt.def)
    {
        DEF.write_file(true);
        Rho.write_file(true);
    }
    if (opt.hirsh)
        Hirsh.write_file(true);

    log2 << " done!" << endl;

    if (opt.esp)
    {
        log2 << "Calculating ESP..." << flush;
        WFN temp = wavy;
        temp.delete_unoccupied_MOs();
        Calc_ESP(ESP, temp, opt.radius, opt.no_date, log2);
        log2 << "Writing cube to Disk..." << flush;
        ESP.write_file(true);
        log2 << "  done!" << endl;
    }
    // return output to cout
    std::cout.rdbuf(_coutbuf);
    log2.close();
    std::cout << "Properties calculation done!" << std::endl;
}

void do_combine_mo(options &opt)
{
    using namespace std;
    WFN wavy1(2);
    WFN wavy2(2);
    WFN wavy3(2);
    wavy1.read_wfn(opt.combine_mo[0], false, cout);
    wavy2.read_wfn(opt.combine_mo[1], false, cout);
    for (int i = 0; i < wavy1.get_ncen(); i++)
    {
        wavy3.push_back_atom(wavy1.get_atom(i));
    }
    for (int i = 0; i < wavy2.get_ncen(); i++)
    {
        wavy3.push_back_atom(wavy2.get_atom(i));
    }
    cout << "In total we have " << wavy3.get_ncen() << " atoms" << endl;

    double MinMax1[6];
    int steps1[3];
    readxyzMinMax_fromWFN(wavy1, MinMax1, steps1, opt.radius, opt.resolution, true);
    double MinMax2[6];
    int steps2[3];
    readxyzMinMax_fromWFN(wavy2, MinMax2, steps2, opt.radius, opt.resolution, true);

    cout << "Read input\nCalculating for MOs ";
    for (int v1 = 0; v1 < opt.cmo1.size(); v1++)
    {
        cout << opt.cmo1[v1] << " ";
    }
    cout << "of fragment 1 and MOs ";
    for (int v1 = 0; v1 < opt.cmo2.size(); v1++)
    {
        cout << opt.cmo2[v1] << " ";
    }
    cout << "of fragment 2" << endl;
    double MinMax[6]{100, 100, 100, -100, -100, -100};
    int steps[3]{0, 0, 0};
    for (int i = 0; i < 3; i++)
    {
        if (MinMax1[i] < MinMax[i])
            MinMax[i] = MinMax1[i];
        if (MinMax1[i + 3] > MinMax[i + 3])
            MinMax[i + 3] = MinMax1[i + 3];
    }
    for (int i = 0; i < 3; i++)
    {
        if (MinMax2[i] < MinMax[i])
            MinMax[i] = MinMax2[i];
        if (MinMax2[i + 3] > MinMax[i + 3])
            MinMax[i + 3] = MinMax2[i + 3];
        steps[i] = (int)ceil(constants::bohr2ang(MinMax[i + 3] - MinMax[i]) / 0.1);
    }
    int counter = 0;
    cube total(steps[0], steps[1], steps[2], 0, true);
    cube MO1(steps[0], steps[1], steps[2], 0, true);
    MO1.give_parent_wfn(wavy3);
    MO1.set_na(wavy3.get_ncen());
    cube MO2(steps[0], steps[1], steps[2], 0, true);
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
            cout << "Running: " << counter << " of " << opt.cmo2.size() * opt.cmo1.size() << endl;
            string filename("");
            MO2.set_zero();
            Calc_MO(MO2, opt.cmo2[j] - 1, wavy2, 40, std::cout);
            cout << "writing files..." << flush;
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
            cout << " ... done!" << endl;
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

                const double PosGrid[3]{i * CubeHirsh.get_vector(0, 0) + j * CubeHirsh.get_vector(0, 1) + k * CubeHirsh.get_vector(0, 2) + CubeHirsh.get_origin(0),
                                        i * CubeHirsh.get_vector(1, 0) + j * CubeHirsh.get_vector(1, 1) + k * CubeHirsh.get_vector(1, 2) + CubeHirsh.get_origin(1),
                                        i * CubeHirsh.get_vector(2, 0) + j * CubeHirsh.get_vector(2, 1) + k * CubeHirsh.get_vector(2, 2) + CubeHirsh.get_origin(2)};

                // bool skip = true;
                double dist = sqrt(pow(PosGrid[0] - wavy.get_atom_coordinate(_atom, 0), 2) + pow(PosGrid[1] - wavy.get_atom_coordinate(_atom, 1), 2) + pow(PosGrid[2] - wavy.get_atom_coordinate(_atom, 2), 2));

                double dens_choice = atom.get_radial_density(dist);
                double temp_val = CubeSpherical.get_value(i, j, k);
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
        {"vdW", dipole_types::vdW}};

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
    const double ax = wavy.get_atom_coordinate(i, 0), ay = wavy.get_atom_coordinate(i, 1), az = wavy.get_atom_coordinate(i, 2), dv = Hirshfeld_atom.get_dv();
    // const int c = wavy.get_atom_charge(i);
    double charge = 0;
    vec origin{0, 0, 0};
    vec bound_atoms;
    for (int j = 0; j < wavy.get_ncen(); j++)
    {
        if (i == j)
            continue;
        double dist = sqrt(pow(ax - wavy.get_atom_coordinate(j, 0), 2) + pow(ay - wavy.get_atom_coordinate(j, 1), 2) + pow(az - wavy.get_atom_coordinate(j, 2), 2));
        double svdW = constants::covalent_radii[wavy.get_atom_charge(i)] + constants::covalent_radii[wavy.get_atom_charge(j)];
        if (dist < 1.1 * svdW)
        {
            bound_atoms.push_back(j);
        }
    }
    const double v[9] = {Hirshfeld_atom.get_vector(0, 0), Hirshfeld_atom.get_vector(0, 1), Hirshfeld_atom.get_vector(0, 2),
                         Hirshfeld_atom.get_vector(1, 0), Hirshfeld_atom.get_vector(1, 1), Hirshfeld_atom.get_vector(1, 2),
                         Hirshfeld_atom.get_vector(2, 0), Hirshfeld_atom.get_vector(2, 1), Hirshfeld_atom.get_vector(2, 2)};
    switch (stringTodipole_types(type))
    {
    case dipole_types::atom:
        origin = {ax, ay, az};
        break;
    case dipole_types::geometry:
        err_not_impl_f("geometry position not yet implemented", std::cout);
        origin = {0, 0, 0};
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
                const double PosGrid[3]{
                    x * v[0] + y * v[1] + z * v[2] + Hirshfeld_atom.get_origin(0),
                    x * v[3] + y * v[4] + z * v[5] + Hirshfeld_atom.get_origin(1),
                    x * v[6] + y * v[7] + z * v[8] + Hirshfeld_atom.get_origin(2)};
                scratch = Hirshfeld_atom.get_value(x, y, z) * dv;
                charge += scratch;
                mu_x += (PosGrid[0] - origin[0]) * scratch;
                mu_y += (PosGrid[1] - origin[1]) * scratch;
                mu_z += (PosGrid[2] - origin[2]) * scratch;
            }
        }
    }
    return {mu_x, mu_y, mu_z, charge};
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
        log2 << opt.cif << " " << opt.resolution << " " << opt.radius << endl;
    readxyzMinMax_fromWFN(wavy, opt.MinMax, opt.NbSteps, opt.radius, opt.resolution, true);
    if (opt.debug)
    {
        log2 << "Resolution: " << opt.resolution << endl;
        log2 << "MinMax:" << endl;
        for (int i = 0; i < 6; i++)
            log2 << setw(14) << scientific << opt.MinMax[i];
        log2 << endl;
        log2 << "Steps:" << endl;
        for (int i = 0; i < 3; i++)
            log2 << setw(14) << scientific << opt.NbSteps[i];
        log2 << endl;
    }
    cube Rho(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), true);
    cube SPHER(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy.get_ncen(), true);

    Rho.give_parent_wfn(wavy);
    SPHER.give_parent_wfn(wavy);
    vec stepsizes{(opt.MinMax[3] - opt.MinMax[0]) / opt.NbSteps[0],
                  (opt.MinMax[4] - opt.MinMax[1]) / opt.NbSteps[1],
                  (opt.MinMax[5] - opt.MinMax[2]) / opt.NbSteps[2]};

    for (int i = 0; i < 3; i++)
    {
        Rho.set_origin(i, opt.MinMax[i]);
        SPHER.set_origin(i, opt.MinMax[i]);
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

    log2 << "Calculating for " << fixed << setprecision(0) << opt.NbSteps[0] * opt.NbSteps[1] * opt.NbSteps[2] << " Gridpoints." << endl;

    log2 << "Calcualting Rho...";
    Calc_Rho(Rho, wavy,  opt.radius, log2, false);
    log2 << " ...done!\nCalcualting spherical Rho...";
    Calc_Spherical_Dens(SPHER, wavy, opt.radius, log2);
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
    log2 << " atom   |  dipole moment x,        y,         z" << endl
         << "======================================" << endl;
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        log2 << setw(3) << i << " (" << constants::atnr2letter(wavy.get_atom_charge(i)) << ") | " << scientific << setprecision(6) << setw(14) << dipole_moments[i][0] << ", " << setw(14) << dipole_moments[i][1] << ", " << setw(14) << dipole_moments[i][2] << endl;
    }
    std::cout << "\n\nProperties calculation done!" << std::endl;
}

vec2 dipole_moments(WFN &wavy, cube &SPHER, double *MinMax, int *NbSteps, int threads, double radius, std::ostream &log2, bool debug)
{
    using namespace std;
    if (debug)
        log2 << "Starting calculation of dipole moment" << endl;
    cube Rho(NbSteps[0], NbSteps[1], NbSteps[2], wavy.get_ncen(), true);

    Rho.give_parent_wfn(wavy);
    vec stepsizes{(MinMax[3] - MinMax[0]) / NbSteps[0],
                  (MinMax[4] - MinMax[1]) / NbSteps[1],
                  (MinMax[5] - MinMax[2]) / NbSteps[2]};

    for (int i = 0; i < 3; i++)
    {
        Rho.set_origin(i, MinMax[i]);
        Rho.set_vector(i, i, stepsizes[i]);
    }
    if (debug)
        log2 << "Origins etc are set up" << endl;
    Rho.set_comment1("Calculated density using NoSpherA2");
    Rho.set_comment2("from " + wavy.get_path().string());
    Rho.set_path((wavy.get_path().parent_path() / wavy.get_path().stem()).string() + "_rho.cube");
    vector<cube> Hirsh(wavy.get_ncen(), Rho);
    vec charges(wavy.get_ncen(), 0);

    log2 << "Calculating for " << fixed << setprecision(0) << NbSteps[0] * NbSteps[1] * NbSteps[2] << " Gridpoints." << endl;

    log2 << "Calcualting Rho...";
    Calc_Rho(Rho, wavy, radius, log2, false);
    log2 << " ...done!\nCalcualting spherical Rho...";
    Calc_Spherical_Dens(SPHER, wavy, radius, log2);
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
        wavy.push_back(WFN(0));
        wavy[i].read_known_wavefunction_format(opt.pol_wfns[i], log2, opt.debug);
    }

    // check that all WFN have the same number of atoms

    if (opt.debug)
        log2 << "Starting calculation of Polarizabilities" << endl;

    if (opt.debug)
        log2 << opt.resolution << " " << opt.radius << endl;
    readxyzMinMax_fromWFN(wavy[0], opt.MinMax, opt.NbSteps, opt.radius, opt.resolution, true);
    if (opt.debug)
    {
        log2 << "Resolution: " << opt.resolution << endl;
        log2 << "MinMax:" << endl;
        for (int i = 0; i < 6; i++)
            log2 << setw(14) << scientific << opt.MinMax[i];
        log2 << endl;
        log2 << "Steps:" << endl;
        for (int i = 0; i < 3; i++)
            log2 << setw(14) << scientific << opt.NbSteps[i];
        log2 << endl;
    }
    cube SPHER(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), true);

    SPHER.give_parent_wfn(wavy[0]);
    vec stepsizes{(opt.MinMax[3] - opt.MinMax[0]) / opt.NbSteps[0],
                  (opt.MinMax[4] - opt.MinMax[1]) / opt.NbSteps[1],
                  (opt.MinMax[5] - opt.MinMax[2]) / opt.NbSteps[2]};

    for (int i = 0; i < 3; i++)
    {
        SPHER.set_origin(i, opt.MinMax[i]);
        SPHER.set_vector(i, i, stepsizes[i]);
    }
    if (opt.debug)
        log2 << "Origins etc are set up" << endl;
    SPHER.set_comment1("Calculated Atomic Hirshfeld deformation density values using NoSpherA2");
    SPHER.set_comment2("from" + wavy[0].get_path().string());
    SPHER.set_path((wavy[0].get_path().parent_path() / wavy[0].get_path().stem()).string() + "_spher.cube");

    log2 << "Calculating for " << fixed << setprecision(0) << opt.NbSteps[0] * opt.NbSteps[1] * opt.NbSteps[2] << " Gridpoints." << endl;

    log2 << "Calcualting spherical Rho...";
    Calc_Spherical_Dens(SPHER, wavy[0], opt.radius, log2, false);
    log2 << " ...done!" << endl;
    vec3 dipoles(7); // 0, +x, -x, +y, -y, +z, -z
    for (int i = 0; i < 7; i++)
    {
        dipoles[i] = dipole_moments(wavy[i], SPHER, opt.MinMax, opt.NbSteps, opt.threads, opt.radius, log2, opt.debug);
    }
    vec3 polarizabilities(wavy[0].get_ncen());
    for (int i = 0; i < wavy[0].get_ncen(); i++)
    {
        polarizabilities[i] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        vec dx = {(dipoles[1][i][0] - dipoles[2][i][0]),
                  (dipoles[3][i][0] - dipoles[4][i][0]),
                  (dipoles[5][i][0] - dipoles[6][i][0])};
        vec dy = {(dipoles[1][i][1] - dipoles[2][i][1]),
                  (dipoles[3][i][1] - dipoles[4][i][1]),
                  (dipoles[5][i][1] - dipoles[6][i][1])};
        vec dz = {(dipoles[1][i][2] - dipoles[2][i][2]),
                  (dipoles[3][i][2] - dipoles[4][i][2]),
                  (dipoles[5][i][2] - dipoles[6][i][2])};
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