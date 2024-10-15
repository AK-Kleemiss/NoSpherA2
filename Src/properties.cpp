#include "wfn_class.h"
#include "convenience.h"
#include "spherical_density.h"
#include "cell.h"
#include "cube.h"
#include "constants.h"

void Calc_Spherical_Dens(
    cube &CubeSpher,
    WFN &wavy,
    int cpus,
    double radius,
    std::ostream &file)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif

    time_point start = get_time();

    progress_bar *progress = new progress_bar{file, 50u, "Calculating Spherical Density"};
    const int step = (int)max(floor(3 * CubeSpher.get_size(0) / 20.0), 1.0);

    vector<Thakkar> atoms;
    for (int a = 0; a < wavy.get_ncen(); a++)
        atoms.push_back(Thakkar(wavy.get_atom_charge(a)));

#pragma omp parallel for schedule(dynamic)
    for (int i = -CubeSpher.get_size(0); i < 2 * CubeSpher.get_size(0); i++)
    {
        for (int j = -CubeSpher.get_size(1); j < 2 * CubeSpher.get_size(1); j++)
            for (int k = -CubeSpher.get_size(2); k < 2 * CubeSpher.get_size(2); k++)
            {

                const double PosGrid[3]{i * CubeSpher.get_vector(0, 0) + j * CubeSpher.get_vector(0, 1) + k * CubeSpher.get_vector(0, 2) + CubeSpher.get_origin(0),
                                        i * CubeSpher.get_vector(1, 0) + j * CubeSpher.get_vector(1, 1) + k * CubeSpher.get_vector(1, 2) + CubeSpher.get_origin(1),
                                        i * CubeSpher.get_vector(2, 0) + j * CubeSpher.get_vector(2, 1) + k * CubeSpher.get_vector(2, 2) + CubeSpher.get_origin(2)};

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                {
                    if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < constants::ang2bohr(radius))
                        skip = false;
                }
                if (skip)
                    continue;

                double dens_all = 0.0;
                double dist;
                for (int a = 0; a < wavy.get_ncen(); a++)
                {
                    dist = sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2));
                    dens_all += atoms[a].get_radial_density(dist);
                    ;
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
        if (i != 0 && i % step == 0
#ifdef _OPENMP
            && omp_get_thread_num() == 0
#endif
        )
            progress->write((i + CubeSpher.get_size(0)) / static_cast<double>(3 * CubeSpher.get_size(0)));
    }
    delete (progress);

    time_point end = get_time();
    if (get_sec(start, end) < 60)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
};

void Calc_Spherical_Dens_no_trans(
  cube& CubeSpher,
  WFN& wavy,
  int cpus,
  double radius,
  std::ostream& file)
{
    using namespace std;
#ifdef _OPENMP
  if (cpus != -1)
  {
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif

  time_point start = get_time();

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating Spherical Density" };
  const int step = (int)max(floor(CubeSpher.get_size(0) / 20.0), 1.0);

  vector<Thakkar> atoms;
  for (int a = 0; a < wavy.get_ncen(); a++)
    atoms.push_back(Thakkar(wavy.get_atom_charge(a)));

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < CubeSpher.get_size(0); i++)
  {
    for (int j = 0; j < CubeSpher.get_size(1); j++)
      for (int k = 0; k < CubeSpher.get_size(2); k++)
      {

        const double PosGrid[3]{ i * CubeSpher.get_vector(0, 0) + j * CubeSpher.get_vector(0, 1) + k * CubeSpher.get_vector(0, 2) + CubeSpher.get_origin(0),
                                i * CubeSpher.get_vector(1, 0) + j * CubeSpher.get_vector(1, 1) + k * CubeSpher.get_vector(1, 2) + CubeSpher.get_origin(1),
                                i * CubeSpher.get_vector(2, 0) + j * CubeSpher.get_vector(2, 1) + k * CubeSpher.get_vector(2, 2) + CubeSpher.get_origin(2) };

        bool skip = true;
        for (int a = 0; a < wavy.get_ncen(); a++)
        {
          if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < constants::ang2bohr(radius))
            skip = false;
        }
        if (skip)
          continue;

        double dens_all = 0.0;
        double dist;
        for (int a = 0; a < wavy.get_ncen(); a++)
        {
          dist = sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2));
          dens_all += atoms[a].get_radial_density(dist);
          ;
        }

        CubeSpher.set_value(i, j, k, dens_all);
      }
    if (i != 0 && i % step == 0
#ifdef _OPENMP
      && omp_get_thread_num() == 0
#endif
      )
      progress->write((i + CubeSpher.get_size(0)) / static_cast<double>(3 * CubeSpher.get_size(0)));
  }
  delete (progress);

  time_point end = get_time();
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
    int cpus,
    double radius,
    std::ostream &file)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif

    time_point start = get_time();

    progress_bar *progress = new progress_bar{file, 50u, "Calculating Deformation Density"};
    const int step = (int)max(floor(3 * CubeDEF.get_size(0) / 20.0), 1.0);

    vector<Thakkar> atoms;
    for (int a = 0; a < wavy.get_ncen(); a++)
        atoms.push_back(Thakkar(wavy.get_atom_charge(a)));

#pragma omp parallel for schedule(dynamic)
    for (int i = -CubeDEF.get_size(0); i < 2 * CubeDEF.get_size(0); i++)
    {
        for (int j = -CubeDEF.get_size(1); j < 2 * CubeDEF.get_size(1); j++)
            for (int k = -CubeDEF.get_size(2); k < 2 * CubeDEF.get_size(2); k++)
            {

                const double PosGrid[3]{i * CubeDEF.get_vector(0, 0) + j * CubeDEF.get_vector(0, 1) + k * CubeDEF.get_vector(0, 2) + CubeDEF.get_origin(0),
                                        i * CubeDEF.get_vector(1, 0) + j * CubeDEF.get_vector(1, 1) + k * CubeDEF.get_vector(1, 2) + CubeDEF.get_origin(1),
                                        i * CubeDEF.get_vector(2, 0) + j * CubeDEF.get_vector(2, 1) + k * CubeDEF.get_vector(2, 2) + CubeDEF.get_origin(2)};

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < constants::ang2bohr(radius))
                        skip = false;
                if (skip)
                    continue;

                double dens_all = 0.0;
                double dist;
                for (int a = 0; a < wavy.get_ncen(); a++)
                {
                    dist = sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2));
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
        if (i != 0 && i % step == 0
#ifdef _OPENMP
            && omp_get_thread_num() == 0
#endif
        )
            progress->write((i + CubeDEF.get_size(0)) / static_cast<double>(3 * CubeDEF.get_size(0)));
    }
    delete (progress);

    time_point end = get_time();
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
    int cpus,
    double radius,
    std::ostream &file)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif

    time_point start = get_time();

    progress_bar *progress = new progress_bar{file, 50u, "Calculating Deformation Density"};
    const int step = (int)max(floor(3 * CubeDEF.get_size(0) / 20.0), 1.0);

#pragma omp parallel for schedule(dynamic)
    for (int i = -CubeDEF.get_size(0); i < 2 * CubeDEF.get_size(0); i++)
    {
        for (int j = -CubeDEF.get_size(1); j < 2 * CubeDEF.get_size(1); j++)
            for (int k = -CubeDEF.get_size(2); k < 2 * CubeDEF.get_size(2); k++)
            {

                const double PosGrid[3]{i * CubeDEF.get_vector(0, 0) + j * CubeDEF.get_vector(0, 1) + k * CubeDEF.get_vector(0, 2) + CubeDEF.get_origin(0),
                                        i * CubeDEF.get_vector(1, 0) + j * CubeDEF.get_vector(1, 1) + k * CubeDEF.get_vector(1, 2) + CubeDEF.get_origin(1),
                                        i * CubeDEF.get_vector(2, 0) + j * CubeDEF.get_vector(2, 1) + k * CubeDEF.get_vector(2, 2) + CubeDEF.get_origin(2)};

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < constants::ang2bohr(radius))
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
        if (i != 0 && i % step == 0
#ifdef _OPENMP
            && omp_get_thread_num() == 0
#endif
        )
            progress->write((i + CubeDEF.get_size(0)) / static_cast<double>(3 * CubeDEF.get_size(0)));
    }
    delete (progress);

    time_point end = get_time();
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
    int cpus,
    double radius,
    int ignore_atom,
    std::ostream &file)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif

    time_point start = get_time();

    progress_bar *progress = new progress_bar{file, 50u, "Calculating Values"};
    const int step = (int)max(floor(3 * CubeHDEF.get_size(0) / 20.0), 1.0);

    vector<Thakkar> atoms;
    for (int a = 0; a < wavy.get_ncen(); a++)
        atoms.push_back(Thakkar(wavy.get_atom_charge(a)));

#pragma omp parallel for schedule(dynamic)
    for (int i = -CubeHDEF.get_size(0); i < 2 * CubeHDEF.get_size(0); i++)
    {
        for (int j = -CubeHDEF.get_size(1); j < 2 * CubeHDEF.get_size(1); j++)
            for (int k = -CubeHDEF.get_size(2); k < 2 * CubeHDEF.get_size(2); k++)
            {

                const double PosGrid[3]{i * CubeHDEF.get_vector(0, 0) + j * CubeHDEF.get_vector(0, 1) + k * CubeHDEF.get_vector(0, 2) + CubeHDEF.get_origin(0),
                                        i * CubeHDEF.get_vector(1, 0) + j * CubeHDEF.get_vector(1, 1) + k * CubeHDEF.get_vector(1, 2) + CubeHDEF.get_origin(1),
                                        i * CubeHDEF.get_vector(2, 0) + j * CubeHDEF.get_vector(2, 1) + k * CubeHDEF.get_vector(2, 2) + CubeHDEF.get_origin(2)};

                bool skip = true;
                if (sqrt(pow(PosGrid[0] - wavy.atoms[ignore_atom].x, 2) + pow(PosGrid[1] - wavy.atoms[ignore_atom].y, 2) + pow(PosGrid[2] - wavy.atoms[ignore_atom].z, 2)) < constants::ang2bohr(radius))
                    skip = false;
                if (skip)
                    continue;

                double dens_choice = 0.0;
                double dens_all = 0.0;
                double dist, temp;
                for (int a = 0; a < wavy.get_ncen(); a++)
                {
                    dist = sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2));
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
        if (i != 0 && i % step == 0
#ifdef _OPENMP
            && omp_get_thread_num() == 0
#endif
        )
            progress->write((i + CubeHDEF.get_size(0)) / static_cast<double>(3 * CubeHDEF.get_size(0)));
    }
    delete (progress);

    time_point end = get_time();
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
    int cpus,
    double radius,
    int ignore_atom,
    std::ostream &file)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif

    time_point start = get_time();

    progress_bar *progress = new progress_bar{file, 50u, "Calculating Values"};
    const int step = (int)max(floor(3 * CubeHDEF.get_size(0) / 20.0), 1.0);

    Thakkar atom(wavy.get_atom_charge(ignore_atom));

#pragma omp parallel for schedule(dynamic)
    for (int i = -CubeHDEF.get_size(0); i < 2 * CubeHDEF.get_size(0); i++)
    {
        for (int j = -CubeHDEF.get_size(1); j < 2 * CubeHDEF.get_size(1); j++)
            for (int k = -CubeHDEF.get_size(2); k < 2 * CubeHDEF.get_size(2); k++)
            {

                const double PosGrid[3]{i * CubeHDEF.get_vector(0, 0) + j * CubeHDEF.get_vector(0, 1) + k * CubeHDEF.get_vector(0, 2) + CubeHDEF.get_origin(0),
                                        i * CubeHDEF.get_vector(1, 0) + j * CubeHDEF.get_vector(1, 1) + k * CubeHDEF.get_vector(1, 2) + CubeHDEF.get_origin(1),
                                        i * CubeHDEF.get_vector(2, 0) + j * CubeHDEF.get_vector(2, 1) + k * CubeHDEF.get_vector(2, 2) + CubeHDEF.get_origin(2)};

                bool skip = true;
                double dist = sqrt(pow(PosGrid[0] - wavy.atoms[ignore_atom].x, 2) + pow(PosGrid[1] - wavy.atoms[ignore_atom].y, 2) + pow(PosGrid[2] - wavy.atoms[ignore_atom].z, 2));
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
        if (i != 0 && i % step == 0
#ifdef _OPENMP
            && omp_get_thread_num() == 0
#endif
        )
            progress->write((i + CubeHDEF.get_size(0)) / static_cast<double>(3 * CubeHDEF.get_size(0)));
    }
    delete (progress);

    time_point end = get_time();
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
    int cpus,
    double radius,
    int ignore_atom,
    std::ostream &file)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif

    time_point start = get_time();

    progress_bar *progress = new progress_bar{file, 50u, "Calculating Values"};
    const int step = (int)max(floor(3 * CubeHirsh.get_size(0) / 20.0), 1.0);

    Thakkar atom(wavy.get_atom_charge(ignore_atom));

#pragma omp parallel for schedule(dynamic)
    for (int i = -CubeHirsh.get_size(0); i < 2 * CubeHirsh.get_size(0); i++)
    {
        for (int j = -CubeHirsh.get_size(1); j < 2 * CubeHirsh.get_size(1); j++)
            for (int k = -CubeHirsh.get_size(2); k < 2 * CubeHirsh.get_size(2); k++)
            {

                const double PosGrid[3]{i * CubeHirsh.get_vector(0, 0) + j * CubeHirsh.get_vector(0, 1) + k * CubeHirsh.get_vector(0, 2) + CubeHirsh.get_origin(0),
                                        i * CubeHirsh.get_vector(1, 0) + j * CubeHirsh.get_vector(1, 1) + k * CubeHirsh.get_vector(1, 2) + CubeHirsh.get_origin(1),
                                        i * CubeHirsh.get_vector(2, 0) + j * CubeHirsh.get_vector(2, 1) + k * CubeHirsh.get_vector(2, 2) + CubeHirsh.get_origin(2)};

                bool skip = true;
                double dist = sqrt(pow(PosGrid[0] - wavy.atoms[ignore_atom].x, 2) + pow(PosGrid[1] - wavy.atoms[ignore_atom].y, 2) + pow(PosGrid[2] - wavy.atoms[ignore_atom].z, 2));
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
        if (i != 0 && i % step == 0
#ifdef _OPENMP
            && omp_get_thread_num() == 0
#endif
        )
            progress->write((i + CubeHirsh.get_size(0)) / static_cast<double>(3 * CubeHirsh.get_size(0)));
    }
    delete (progress);

    time_point end = get_time();
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
    int cpus,
    double radius,
    std::ostream &file)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif

    time_point start = get_time();

    progress_bar *progress = new progress_bar{file, 50u, "Calculating Values"};
    const int step = (int)max(floor(CubeRho.get_size(0) * 3 / 20.0), 1.0);

#pragma omp parallel for schedule(dynamic)
    for (int i = -CubeRho.get_size(0); i < 2 * CubeRho.get_size(0); i++)
    {
        for (int j = -CubeRho.get_size(1); j < 2 * CubeRho.get_size(1); j++)
            for (int k = -CubeRho.get_size(2); k < 2 * CubeRho.get_size(2); k++)
            {

                double PosGrid[3]{
                    i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0),
                    i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1),
                    i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2)},
                    Rho = 0;

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < constants::ang2bohr(radius))
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
        if (i != 0 && i % step == 0)
            progress->write((i + CubeRho.get_size(0)) / static_cast<double>(CubeRho.get_size(0) * 3));
    }
    delete (progress);

    time_point end = get_time();
    if (get_sec(start, end) < 60)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
};

void Calc_Rho_no_trans(
    cube &CubeRho,
    WFN &wavy,
    int cpus,
    double radius,
    std::ostream &file,
    bool print)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif

    time_point start = get_time();

    progress_bar* progress = nullptr;
    if (print)
        progress = new progress_bar{file, 50u, "Calculating Values"};
    const int step = (int)max(floor(CubeRho.get_size(0) / 20.0), 1.0);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < CubeRho.get_size(0); i++)
    {
        for (int j = 0; j < CubeRho.get_size(1); j++)
            for (int k = 0; k < CubeRho.get_size(2); k++)
            {

                double PosGrid[3]{
                    i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0),
                    i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1),
                    i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2)},
                    Rho = 0;

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < constants::ang2bohr(radius))
                        skip = false;
                if (skip)
                    continue;

                Rho = wavy.compute_dens(PosGrid[0], PosGrid[1], PosGrid[2]);

                CubeRho.set_value(i, j, k, Rho);
            }
        if (print && i != 0 && i % step == 0 && progress != nullptr)
            progress->write((i + CubeRho.get_size(0)) / static_cast<double>(CubeRho.get_size(0) * 3));
    }
    if(progress != nullptr)
        delete (progress);

    time_point end = get_time();
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
    int cpus,
    std::ostream &file)
{
    using namespace std;
    time_point start = get_time();

    progress_bar *progress = new progress_bar{file, 50u, "Calculating Rho"};
    const int step = (int)max(floor(CubeRho.get_size(0) * 3 / 20.0), 1.0);

#pragma omp parallel shared(CubeRho) num_threads(cpus)
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

                    double PosGrid[3]{
                        i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0),
                        i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1),
                        i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2)};

                    CubeRho.set_value(i, j, k, wavy.compute_dens(PosGrid[0], PosGrid[1], PosGrid[2], d, phi));
                }
            if (i != 0 && i % step == 0)
                progress->write((i) / static_cast<double>(CubeRho.get_size(0)));
        }
    }
    delete (progress);

    time_point end = get_time();
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
    int cpus,
    int MO,
    std::ostream &file,
    bool nodate)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif

    time_point start = get_time();

    progress_bar* progress = NULL;
    if (!nodate)
        progress = new progress_bar{file, 50u, "Calculating Values"};
    const int step = (int)max(floor(CubeMO.get_size(0) / 20.0), 1.0);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < CubeMO.get_size(0); i++)
    {
        for (int j = 0; j < CubeMO.get_size(1); j++)
            for (int k = 0; k < CubeMO.get_size(2); k++)
            {

                double PosGrid[3]{
                    i * CubeMO.get_vector(0, 0) + j * CubeMO.get_vector(0, 1) + k * CubeMO.get_vector(0, 2) + CubeMO.get_origin(0),
                    i * CubeMO.get_vector(1, 0) + j * CubeMO.get_vector(1, 1) + k * CubeMO.get_vector(1, 2) + CubeMO.get_origin(1),
                    i * CubeMO.get_vector(2, 0) + j * CubeMO.get_vector(2, 1) + k * CubeMO.get_vector(2, 2) + CubeMO.get_origin(2)};

                CubeMO.set_value(i, j, k, wavy.compute_MO_spherical(PosGrid[0], PosGrid[1], PosGrid[2], MO));
            }
        if (i != 0 && i % step == 0&& !nodate)
            progress->write((i) / static_cast<double>(CubeMO.get_size(0)));
    }
    if (!nodate) {
        delete (progress);

        time_point end = get_time();
        if (get_sec(start, end) < 60)
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
        else if (get_sec(start, end) < 3600)
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
        else
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
    }
};

void Calc_S_Rho(
    cube& Cube_S_Rho,
    WFN& wavy,
    int cpus,
    std::ostream& file,
    bool& nodate)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif

    time_point start = get_time();
    progress_bar* progress = NULL;
    if(!nodate)
        progress = new progress_bar{ file, 50u, "Calculating Values" };
    const int step = (int)max(floor(Cube_S_Rho.get_size(0) / 20.0), 1.0);

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

                double PosGrid[3]{
                    i * Cube_S_Rho.get_vector(0, 0) + j * Cube_S_Rho.get_vector(0, 1) + k * Cube_S_Rho.get_vector(0, 2) + Cube_S_Rho.get_origin(0),
                    i * Cube_S_Rho.get_vector(1, 0) + j * Cube_S_Rho.get_vector(1, 1) + k * Cube_S_Rho.get_vector(1, 2) + Cube_S_Rho.get_origin(1),
                    i * Cube_S_Rho.get_vector(2, 0) + j * Cube_S_Rho.get_vector(2, 1) + k * Cube_S_Rho.get_vector(2, 2) + Cube_S_Rho.get_origin(2) };

                Cube_S_Rho.set_value(i, j, k, wavy.compute_spin_dens(PosGrid[0], PosGrid[1], PosGrid[2], d, phi));
            }
        if (i != 0 && i % step == 0 && !nodate)
            progress->write((i) / static_cast<double>(Cube_S_Rho.get_size(0)));
    }
    if (!nodate) {
        delete (progress);

        time_point end = get_time();
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
    int cpus,
    double radius,
    std::ostream &file,
    bool test)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif

    time_point start = get_time();

    progress_bar *progress = NULL;
    if (!test)
        progress = new progress_bar{file, 50u, "Calculating Values"};
    const int step = (int)max(floor(CubeRho.get_size(0) * 3 / 20.0), 1.0);

#pragma omp parallel for schedule(dynamic)
    for (int i = -CubeRho.get_size(0); i < 2 * CubeRho.get_size(0); i++)
    {
        for (int j = -CubeRho.get_size(1); j < 2 * CubeRho.get_size(1); j++)
            for (int k = -CubeRho.get_size(2); k < 2 * CubeRho.get_size(2); k++)
            {

                double PosGrid[3]{
                    i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0),
                    i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1),
                    i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2)},
                    Rho = 0,
                    Grad = 0,
                    Elf = 0,
                    Eli = 0,
                    Lap = 0,
                    Hess[9]{0, 0, 0, 0, 0, 0, 0, 0, 0};

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < constants::ang2bohr(radius))
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
        {
            if (i != 0 && i % step == 0)
                progress->write((i + CubeRho.get_size(0)) / static_cast<double>(CubeRho.get_size(0) * 3));
        }
    }
    if (!test)
    {
        delete (progress);
        time_point end = get_time();
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
    int cpus,
    double radius,
    bool no_date,
    std::ostream &file)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif
    time_point start = get_time();

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
            d2[i][j] = pow(wavy.atoms[i].x - wavy.atoms[j].x, 2) + pow(wavy.atoms[i].y - wavy.atoms[j].y, 2) + pow(wavy.atoms[i].z - wavy.atoms[j].z, 2);
        }
    }

    progress_bar *progress = NULL;
    if (!no_date)
        progress = new progress_bar{file, 50u, "Calculating ESP"};
    const int step = (int)max(floor(CubeESP.get_size(0) * 3 / 20.0), 1.0);

#pragma omp parallel for schedule(dynamic)
    for (int i = -CubeESP.get_size(0); i < 2 * CubeESP.get_size(0); i++)
    {
        double temp;
        int temp_i, temp_j, temp_k;

        for (int j = -CubeESP.get_size(1); j < 2 * CubeESP.get_size(1); j++)
            for (int k = -CubeESP.get_size(2); k < 2 * CubeESP.get_size(2); k++)
            {
                const double PosGrid[3]{
                    i * CubeESP.get_vector(0, 0) + j * CubeESP.get_vector(0, 1) + k * CubeESP.get_vector(0, 2) + CubeESP.get_origin(0),
                    i * CubeESP.get_vector(1, 0) + j * CubeESP.get_vector(1, 1) + k * CubeESP.get_vector(1, 2) + CubeESP.get_origin(1),
                    i * CubeESP.get_vector(2, 0) + j * CubeESP.get_vector(2, 1) + k * CubeESP.get_vector(2, 2) + CubeESP.get_origin(2)};

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < constants::ang2bohr(radius))
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
        {
            if (i != 0 && i % step == 0)
                progress->write((i + CubeESP.get_size(0)) / static_cast<double>(CubeESP.get_size(0) * 3));
        }
    }
    if (!no_date)
    {
        delete (progress);

        time_point end = get_time();
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
    int cpus,
    double radius,
    std::ostream &file)
{
    using namespace std;
    err_checkf(mo <= wavy.get_nmo(), to_string(mo) + " bigger MO selected than " + to_string(wavy.get_nmo()) + " contained in the wavefunctions!", file);
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif
    time_point start = get_time();

    progress_bar *progress = new progress_bar{file, 50u, "Calculating MO"};
    const int step = (int)max(floor(CubeMO.get_size(0) * 3 / 20.0), 1.0);

#pragma omp parallel for schedule(dynamic)
    for (int i = -CubeMO.get_size(0); i < 2 * CubeMO.get_size(0); i++)
    {
        for (int j = -CubeMO.get_size(1); j < 2 * CubeMO.get_size(1); j++)
            for (int k = -CubeMO.get_size(2); k < 2 * CubeMO.get_size(2); k++)
            {

                const double PosGrid[3]{
                    i * CubeMO.get_vector(0, 0) + j * CubeMO.get_vector(0, 1) + k * CubeMO.get_vector(0, 2) + CubeMO.get_origin(0),
                    i * CubeMO.get_vector(1, 0) + j * CubeMO.get_vector(1, 1) + k * CubeMO.get_vector(1, 2) + CubeMO.get_origin(1),
                    i * CubeMO.get_vector(2, 0) + j * CubeMO.get_vector(2, 1) + k * CubeMO.get_vector(2, 2) + CubeMO.get_origin(2)};
                double MO = 0;

                bool skip = true;
                for (int a = 0; a < wavy.get_ncen(); a++)
                    if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < constants::ang2bohr(radius))
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
        if (i != 0 && i % step == 0)
            progress->write((i + CubeMO.get_size(0)) / static_cast<double>(CubeMO.get_size(0) * 3));
    }
    delete (progress);

    time_point end = get_time();
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
        log2 << build_date();
    }
    log2.flush();

    err_checkf(opt.wfn != "", "Error, no wfn file specified!", log2);
    WFN wavy(0);
    wavy.read_known_wavefunction_format(opt.wfn, log2, opt.debug);
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
        cell_matrix[i].resize(3);
    if (opt.debug)
        log2 << opt.cif << " " << opt.resolution << " " << opt.radius << endl;
    if (opt.cif != "")
        readxyzMinMax_fromCIF(opt.cif, opt.MinMax, opt.NbSteps, cell_matrix, opt.resolution, log2, opt.debug);
    else {
        readxyzMinMax_fromWFN(wavy, opt.MinMax, opt.NbSteps, opt.radius, opt.resolution);
        for (int i = 0; i < 3; i++)
            cell_matrix[i][i] = opt.resolution;
    }
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
        log2 << "Origins etc are set up" << endl;
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
    Rho.set_comment2("from " + wavy.get_path());
    RDG.set_comment2("from " + wavy.get_path());
    Elf.set_comment2("from " + wavy.get_path());
    Eli.set_comment2("from " + wavy.get_path());
    Lap.set_comment2("from " + wavy.get_path());
    ESP.set_comment2("from " + wavy.get_path());
    MO.set_comment2("from" + wavy.get_path());
    HDEF.set_comment2("from" + wavy.get_path());
    DEF.set_comment2("from" + wavy.get_path());
    Hirsh.set_comment2("from" + wavy.get_path());
    S_Rho.set_comment2("from" + wavy.get_path());
    Rho.path = get_basename_without_ending(wavy.get_path()) + "_rho.cube";
    RDG.path = get_basename_without_ending(wavy.get_path()) + "_rdg.cube";
    Elf.path = get_basename_without_ending(wavy.get_path()) + "_elf.cube";
    Eli.path = get_basename_without_ending(wavy.get_path()) + "_eli.cube";
    Lap.path = get_basename_without_ending(wavy.get_path()) + "_lap.cube";
    ESP.path = get_basename_without_ending(wavy.get_path()) + "_esp.cube";
    DEF.path = get_basename_without_ending(wavy.get_path()) + "_def.cube";
    Hirsh.path = get_basename_without_ending(wavy.get_path()) + "_hirsh.cube";
    S_Rho.path = get_basename_without_ending(wavy.get_path()) + "_s_rho.cube";

    if (opt.debug)
    {
        log2 << "Status: " << opt.hdef << opt.def << opt.hirsh << opt.lap << opt.eli << opt.elf << opt.rdg << opt.esp << endl;
        log2 << "Everything is set up; starting calculation..." << endl;
    }

    log2 << "Calculating for " << fixed << setprecision(0) << opt.NbSteps[0] * opt.NbSteps[1] * opt.NbSteps[2] << " Gridpoints." << endl;

    if (opt.MOs.size() != 0)
        for (int i = 0; i < opt.MOs.size(); i++)
        {
            log2 << "Calcualting MO: " << opt.MOs[i] << endl;
            MO.set_zero();
            MO.path = get_basename_without_ending(wavy.get_path()) + "_MO_" + to_string(opt.MOs[i]) + ".cube";
            Calc_MO(MO, opt.MOs[i], wavy, opt.threads, opt.radius, log2);
            MO.write_file(true);
        }

    if (opt.hdef || opt.def || opt.hirsh)
    {
        log2 << "Calcualting Rho...";
        Calc_Rho(Rho, wavy, opt.threads, opt.radius, log2);
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
            Calc_Spherical_Dens(temp, wavy, opt.threads, opt.radius, log2);
            log2 << " ...done!" << endl;
        }

        if (opt.def)
        {
            log2 << "Calculating static deformation density...";
            if (opt.hdef)
                Calc_Static_Def(DEF, Rho, temp, wavy, opt.threads, opt.radius, log2);
            else
                Calc_Static_Def(DEF, Rho, wavy, opt.threads, opt.radius, log2);
            log2 << " ...done!" << endl;
        }

        if (opt.hdef)
        {
            for (int a = 0; a < wavy.get_ncen(); a++)
            {
                log2 << "Calcualting Hirshfeld deformation density for atom: " << a << endl;
                HDEF.path = get_basename_without_ending(wavy.get_path()) + "_HDEF_" + to_string(a) + ".cube";
                Calc_Hirshfeld(HDEF, Rho, temp, wavy, opt.threads, opt.radius, a, log2);
                HDEF.write_file(true);
                HDEF.set_zero();
            }
        }

        if (opt.hirsh)
        {
            log2 << "Calcualting Hirshfeld density for atom: " << opt.hirsh_number << endl;
            Calc_Hirshfeld_atom(Hirsh, Rho, temp, wavy, opt.threads, opt.radius, opt.hirsh_number, log2);
            log2 << "..done!" << endl;
        }
    }

    if (opt.lap || opt.eli || opt.elf || opt.rdg || opt.esp)
        Calc_Prop(Rho, RDG, Elf, Eli, Lap, ESP, wavy, opt.threads, opt.radius, log2, opt.no_date);

    if (opt.s_rho)
        Calc_S_Rho(S_Rho, wavy, opt.threads, log2, opt.no_date);

    log2 << "Writing cubes to Disk..." << flush;
    if (opt.rdg)
    {
        Rho.path = get_basename_without_ending(wavy.get_path()) + "_signed_rho.cube";
        Rho.write_file(true);
        Rho.path = get_basename_without_ending(wavy.get_path()) + "_rho.cube";
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
        Calc_ESP(ESP, temp, opt.threads, opt.radius, opt.no_date, log2);
        log2 << "Writing cube to Disk..." << flush;
        ESP.write_file(true);
        log2 << "  done!" << endl;
    }
    // return output to cout
    std::cout.rdbuf(_coutbuf);
    log2.close();
    std::cout<< "Properties calculation done!" << std::endl;
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
        Calc_MO(MO1, opt.cmo1[v1] - 1, wavy1, -1, 400, std::cout);
        for (int j = 0; j < opt.cmo2.size(); j++)
        {
            counter++;
            cout << "Running: " << counter << " of " << opt.cmo2.size() * opt.cmo1.size() << endl;
            string filename("");
            MO2.set_zero();
            Calc_MO(MO2, opt.cmo2[j] - 1, wavy2, -1, 400, std::cout);
            cout << "writing files..." << flush;
            filename = get_basename_without_ending(wavy1.get_path()) + "_" + std::to_string(opt.cmo1[v1]) + "+" + get_basename_without_ending(wavy2.get_path()) + "_" + std::to_string(opt.cmo2[j]) + ".cube";
            fns.push_back(filename);
            total.set_zero();
            total = MO1;
            total += MO2;
            total.write_file(filename, false);
            filename = get_basename_without_ending(wavy1.get_path()) + "_" + std::to_string(opt.cmo1[v1]) + "-" + get_basename_without_ending(wavy2.get_path()) + "_" + std::to_string(opt.cmo2[j]) + ".cube";
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
  cube& CubeHirsh,
  cube& CubeRho,
  cube& CubeSpherical,
  WFN& wavy,
  int cpus,
  int _atom,
  std::ostream& file)
{
    (void)file;
    using namespace std;
#ifdef _OPENMP
  if (cpus != -1)
  {
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif
  Thakkar atom(wavy.get_atom_charge(_atom));

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < CubeHirsh.get_size(0); i++)
  {
    for (int j = 0; j < CubeHirsh.get_size(1); j++)
      for (int k = 0; k < CubeHirsh.get_size(2); k++)
      {

        const double PosGrid[3]{ i * CubeHirsh.get_vector(0, 0) + j * CubeHirsh.get_vector(0, 1) + k * CubeHirsh.get_vector(0, 2) + CubeHirsh.get_origin(0),
                                i * CubeHirsh.get_vector(1, 0) + j * CubeHirsh.get_vector(1, 1) + k * CubeHirsh.get_vector(1, 2) + CubeHirsh.get_origin(1),
                                i * CubeHirsh.get_vector(2, 0) + j * CubeHirsh.get_vector(2, 1) + k * CubeHirsh.get_vector(2, 2) + CubeHirsh.get_origin(2) };

        //bool skip = true;
        double dist = sqrt(pow(PosGrid[0] - wavy.atoms[_atom].x, 2) + pow(PosGrid[1] - wavy.atoms[_atom].y, 2) + pow(PosGrid[2] - wavy.atoms[_atom].z, 2));

        double dens_choice = atom.get_radial_density(dist);
        double temp_val = CubeSpherical.get_value(i, j, k);
        if (temp_val != 0)
          CubeHirsh.set_value(i, j, k, (dens_choice / temp_val * CubeRho.get_value(i, j, k)));
      }
  }
};

vec calc_dipole_for_atom(WFN& wavy, const int& i, cube& Hirshfeld_atom) {
  double mu_x = 0, mu_y = 0, mu_z = 0;
  double scratch = 0;
  const double ax = wavy.get_atom_coordinate(i, 0), ay = wavy.get_atom_coordinate(i, 1), az = wavy.get_atom_coordinate(i, 2), dv = Hirshfeld_atom.get_dv();
  //const int c = wavy.get_atom_charge(i);
  double charge = 0;
  const double v[9] = { Hirshfeld_atom.get_vector(0, 0), Hirshfeld_atom.get_vector(0, 1), Hirshfeld_atom.get_vector(0, 2),
												 Hirshfeld_atom.get_vector(1, 0), Hirshfeld_atom.get_vector(1, 1), Hirshfeld_atom.get_vector(1, 2),
												 Hirshfeld_atom.get_vector(2, 0), Hirshfeld_atom.get_vector(2, 1), Hirshfeld_atom.get_vector(2, 2) };
#pragma omp parallel for reduction(+:mu_x, mu_y, mu_z, charge) private(scratch)
  for (int x = 0; x < Hirshfeld_atom.get_size(0); x++) {
    for (int y = 0; y < Hirshfeld_atom.get_size(1); y++) {
      for (int z = 0; z < Hirshfeld_atom.get_size(2); z++) {
        const double PosGrid[3]{
                    x * v[0] + y * v[1] + z * v[2] + Hirshfeld_atom.get_origin(0),
                    x * v[3] + y * v[4] + z * v[5] + Hirshfeld_atom.get_origin(1),
                    x * v[6] + y * v[7] + z * v[8] + Hirshfeld_atom.get_origin(2) };
        //scratch = -Hirshfeld_atom.get_value(x, y, z) * dv + c;
        //scratch = c;
        scratch = Hirshfeld_atom.get_value(x, y, z) * dv;
        charge += scratch;
        mu_x += (PosGrid[0] - ax) * scratch;
        mu_y += (PosGrid[1] - ay) * scratch;
        mu_z += (PosGrid[2] - az) * scratch;
      }
    }
  }
  return { mu_x, mu_y, mu_z ,charge};
}

void dipole_moments(options& opt, std::ostream& log2)
{
  using namespace std;
  log2 << NoSpherA2_message(opt.no_date);
  if (!opt.no_date)
    log2 << build_date();
  err_checkf(opt.wfn != "", "Error, no wfn file specified!", log2);
  WFN wavy(0);
  wavy.read_known_wavefunction_format(opt.wfn, log2, opt.debug);
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
  vec stepsizes{ (opt.MinMax[3] - opt.MinMax[0]) / opt.NbSteps[0],
  (opt.MinMax[4] - opt.MinMax[1]) / opt.NbSteps[1],
  (opt.MinMax[5] - opt.MinMax[2]) / opt.NbSteps[2] };

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
  Rho.set_comment2("from " + wavy.get_path());
  SPHER.set_comment2("from" + wavy.get_path());
  Rho.path = get_basename_without_ending(wavy.get_path()) + "_rho.cube";
  SPHER.path = get_basename_without_ending(wavy.get_path()) + "_spher.cube";
  cube Hirsh = Rho;
  Hirsh.calc_dv();
  Hirsh.give_parent_wfn(wavy);

  log2 << "Calculating for " << fixed << setprecision(0) << opt.NbSteps[0] * opt.NbSteps[1] * opt.NbSteps[2] << " Gridpoints." << endl;

  log2 << "Calcualting Rho...";
  Calc_Rho_no_trans(Rho, wavy, opt.threads, opt.radius, log2, false);
  log2 << " ...done!\nCalcualting spherical Rho...";
  Calc_Spherical_Dens_no_trans(SPHER, wavy, opt.threads, opt.radius, log2);
  log2 << " ...done!" << endl;
  vec2 dipole_moments;
  for (int i = 0; i < wavy.get_ncen(); i++) {
      Hirsh.set_zero();
      log2 << "Calcualting Hirshfeld density for atom: " << i << endl;
      Calc_Hirshfeld_atom_2(Hirsh, Rho, SPHER, wavy, opt.threads, i, log2);
      //Hirsh.path = get_basename_without_ending(wavy.get_path()) + "_h" + to_string(i) + ".cube";
      //Hirsh.write_file(true);
      log2 << "..done!" << endl;
      dipole_moments.push_back(calc_dipole_for_atom(wavy, i, Hirsh));
  }
  log2 << " atom   |  dipole moment x,        y,         z" << endl << "======================================" << endl;
  for (int i = 0; i < wavy.get_ncen(); i++) {
    log2 << setw(3) << i << " (" << constants::atnr2letter(wavy.get_atom_charge(i)) << ") | " << scientific << setprecision(6) << setw(14) << dipole_moments[i][0] << ", " << setw(14) << dipole_moments[i][1] << ", " << setw(14) << dipole_moments[i][2] << endl;
  }
  std::cout << "\n\nProperties calculation done!" << std::endl;
}

vec2 dipole_moments(WFN& wavy, cube& SPHER, double* MinMax, int* NbSteps, int threads, double radius, std::ostream& log2, bool debug)
{
    using namespace std;
  if (debug)
    log2 << "Starting calculation of dipole moment" << endl;
  cube Rho(NbSteps[0], NbSteps[1], NbSteps[2], wavy.get_ncen(), true);

  Rho.give_parent_wfn(wavy);
  vec stepsizes{ (MinMax[3] - MinMax[0]) / NbSteps[0],
  (MinMax[4] - MinMax[1]) / NbSteps[1],
  (MinMax[5] - MinMax[2]) / NbSteps[2] };

  for (int i = 0; i < 3; i++)
  {
    Rho.set_origin(i, MinMax[i]);
    Rho.set_vector(i, i, stepsizes[i]);
  }
  if (debug)
    log2 << "Origins etc are set up" << endl;
  Rho.set_comment1("Calculated density using NoSpherA2");
  Rho.set_comment2("from " + wavy.get_path());
  Rho.path = get_basename_without_ending(wavy.get_path()) + "_rho.cube";
  cube Hirsh = Rho;
  Hirsh.calc_dv();
  Hirsh.give_parent_wfn(wavy);

  log2 << "Calculating electron density for " << fixed << setprecision(0) << NbSteps[0] * NbSteps[1] * NbSteps[2] << " Gridpoints." << endl;
  Calc_Rho_no_trans(Rho, wavy, threads, radius, log2, false);
  vec2 dipole_moments;
  log2 << "Calcualting Hirshfeld density for atom: " << flush;
  for (int i = 0; i < wavy.get_ncen(); i++) {
      Hirsh.set_zero();
      log2 << i << " " << flush;
      Calc_Hirshfeld_atom_2(Hirsh, Rho, SPHER, wavy, threads, i, log2);
      dipole_moments.push_back(calc_dipole_for_atom(wavy, i, Hirsh));
  }
  log2 << "...done!" << endl;
  log2 << " atom   |    charge    | dipole moment x,        y,         z" << endl << "===================================================" << endl;
  for (int i = 0; i < wavy.get_ncen(); i++) {
    log2 << setw(3) << i <<" (" << constants::atnr2letter(wavy.get_atom_charge(i))<< ") |" << scientific << setprecision(6) << setw(13) << dipole_moments[i][3]-wavy.atoms[i].charge << " | " << scientific << setprecision(6) << setw(14) << dipole_moments[i][0] << ", " << setw(14) << dipole_moments[i][1] << ", " << setw(14) << dipole_moments[i][2] << endl;
  }
  return dipole_moments;
}

void polarizabilities(options& opt, std::ostream& log2)
{
    using namespace std;
  std::vector<WFN> wavy;
  for (int i = 0; i < 7; i++) {
    wavy.push_back(WFN(0));
    wavy[i].read_known_wavefunction_format(opt.pol_wfns[i], log2, opt.debug);
  }

  //check that all WFN have the same number of atoms

  
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
  vec stepsizes{ (opt.MinMax[3] - opt.MinMax[0]) / opt.NbSteps[0],
  (opt.MinMax[4] - opt.MinMax[1]) / opt.NbSteps[1],
  (opt.MinMax[5] - opt.MinMax[2]) / opt.NbSteps[2] };

  for (int i = 0; i < 3; i++)
  {
    SPHER.set_origin(i, opt.MinMax[i]);
    SPHER.set_vector(i, i, stepsizes[i]);
  }
  if (opt.debug)
    log2 << "Origins etc are set up" << endl;
  SPHER.set_comment1("Calculated Atomic Hirshfeld deformation density values using NoSpherA2");
  SPHER.set_comment2("from" + wavy[0].get_path());
  SPHER.path = get_basename_without_ending(wavy[0].get_path()) + "_spher.cube";

  log2 << "Calculating for " << fixed << setprecision(0) << opt.NbSteps[0] * opt.NbSteps[1] * opt.NbSteps[2] << " Gridpoints." << endl;

  log2 << "Calcualting spherical Rho...";
  Calc_Spherical_Dens_no_trans(SPHER, wavy[0], opt.threads, opt.radius, log2);
  log2 << " ...done!" << endl;
  vec3 dipoles(7); // 0, +x, -x, +y, -y, +z, -z
  for (int i = 0; i < 7; i++) {
    dipoles[i] = dipole_moments(wavy[i], SPHER, opt.MinMax, opt.NbSteps, opt.threads, opt.radius, log2, opt.debug);
  }
  vec3 polarizabilities(wavy[0].get_ncen());
  for (int i = 0; i < wavy[0].get_ncen(); i++) {
    polarizabilities[i] = { { 0, 0, 0 }, {0,0,0}, {0,0,0} };
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
  for (int i = 0; i < wavy[0].get_ncen(); i++) {
		log2 << setw(3) << i << " (" << constants::atnr2letter(wavy[0].get_atom_charge(i)) << ") |" 
      << scientific << setprecision(6) << setw(13) << dipoles[0][i][3]-wavy[0].atoms[i].charge << " |"
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