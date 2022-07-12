#include "wfn_class.h"
#include "convenience.h"
#include "spherical_density.h"
#include "cell.h"
#include "cube.h"

using namespace std;

void Calc_Spherical_Dens(
  cube& CubeSpher,
  WFN& wavy,
  int cpus,
  double radius,
  ofstream& file
)
{
#ifdef _OPENMP
  if (cpus != -1) {
    omp_set_num_threads(cpus);
    omp_set_dynamic(0);
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif

  time_t start;
  time(&start);

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating Spherical Density" };
  const int step = max(floor(3 * CubeSpher.get_size(0) / 20), 1.0);

  vector<Thakkar> atoms;
  for (int a = 0; a < wavy.get_ncen(); a++)
    atoms.push_back(Thakkar(wavy.get_atom_charge(a)));

#pragma omp parallel for schedule(dynamic)
  for (int i = -CubeSpher.get_size(0); i < 2 * CubeSpher.get_size(0); i++) {
    for (int j = -CubeSpher.get_size(1); j < 2 * CubeSpher.get_size(1); j++)
      for (int k = -CubeSpher.get_size(2); k < 2 * CubeSpher.get_size(2); k++) {

        double PosGrid[3];

        PosGrid[0] = i * CubeSpher.get_vector(0, 0) + j * CubeSpher.get_vector(0, 1) + k * CubeSpher.get_vector(0, 2) + CubeSpher.get_origin(0);
        PosGrid[1] = i * CubeSpher.get_vector(1, 0) + j * CubeSpher.get_vector(1, 1) + k * CubeSpher.get_vector(1, 2) + CubeSpher.get_origin(1);
        PosGrid[2] = i * CubeSpher.get_vector(2, 0) + j * CubeSpher.get_vector(2, 1) + k * CubeSpher.get_vector(2, 2) + CubeSpher.get_origin(2);

        bool skip = true;
        for (int a = 0; a < wavy.get_ncen(); a++) {
          if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius / 0.52)
            skip = false;
        }
        if (skip)
          continue;

        double dens_all = 0.0;
        double dist;
        for (int a = 0; a < wavy.get_ncen(); a++) {
          dist = sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2));
          dens_all += atoms[a].get_radial_density(dist);;
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
    if (i != 0 && i % step == 0 && omp_get_thread_num() == 0)
      progress->write((i + CubeSpher.get_size(0)) / double(3 * CubeSpher.get_size(0)));
  }
  delete(progress);

  time_t end;
  time(&end);
  if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
  else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
  else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Static_Def(
  cube& CubeDEF,
  cube& CubeRho,
  WFN& wavy,
  int cpus,
  double radius,
  ofstream& file
)
{
#ifdef _OPENMP
  if (cpus != -1) {
    omp_set_num_threads(cpus);
    omp_set_dynamic(0);
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif

  time_t start;
  time(&start);

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating Deformation Density" };
  const int step = max(floor(3 * CubeDEF.get_size(0) / 20), 1.0);

  vector<Thakkar> atoms;
  for (int a = 0; a < wavy.get_ncen(); a++)
    atoms.push_back(Thakkar(wavy.get_atom_charge(a)));

#pragma omp parallel for schedule(dynamic)
  for (int i = -CubeDEF.get_size(0); i < 2 * CubeDEF.get_size(0); i++) {
    for (int j = -CubeDEF.get_size(1); j < 2 * CubeDEF.get_size(1); j++)
      for (int k = -CubeDEF.get_size(2); k < 2 * CubeDEF.get_size(2); k++) {

        double PosGrid[3];

        PosGrid[0] = i * CubeDEF.get_vector(0, 0) + j * CubeDEF.get_vector(0, 1) + k * CubeDEF.get_vector(0, 2) + CubeDEF.get_origin(0);
        PosGrid[1] = i * CubeDEF.get_vector(1, 0) + j * CubeDEF.get_vector(1, 1) + k * CubeDEF.get_vector(1, 2) + CubeDEF.get_origin(1);
        PosGrid[2] = i * CubeDEF.get_vector(2, 0) + j * CubeDEF.get_vector(2, 1) + k * CubeDEF.get_vector(2, 2) + CubeDEF.get_origin(2);

        bool skip = true;
        for (int a = 0; a < wavy.get_ncen(); a++)
          if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius / 0.52)
            skip = false;
        if (skip)
          continue;

        double dens_all = 0.0;
        double dist;
        for (int a = 0; a < wavy.get_ncen(); a++) {
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
    if (i != 0 && i % step == 0 && omp_get_thread_num() == 0)
      progress->write((i + CubeDEF.get_size(0)) / double(3 * CubeDEF.get_size(0)));
  }
  delete(progress);

  time_t end;
  time(&end);
  if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
  else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
  else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Static_Def(
  cube& CubeDEF,
  cube& CubeRho,
  cube& CubeSpher,
  WFN& wavy,
  int cpus,
  double radius,
  ofstream& file
)
{
#ifdef _OPENMP
  if (cpus != -1) {
    omp_set_num_threads(cpus);
    omp_set_dynamic(0);
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif

  time_t start;
  time(&start);

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating Deformation Density" };
  const int step = max(floor(3 * CubeDEF.get_size(0) / 20), 1.0);

#pragma omp parallel for schedule(dynamic)
  for (int i = -CubeDEF.get_size(0); i < 2 * CubeDEF.get_size(0); i++) {
    for (int j = -CubeDEF.get_size(1); j < 2 * CubeDEF.get_size(1); j++)
      for (int k = -CubeDEF.get_size(2); k < 2 * CubeDEF.get_size(2); k++) {

        double PosGrid[3];

        PosGrid[0] = i * CubeDEF.get_vector(0, 0) + j * CubeDEF.get_vector(0, 1) + k * CubeDEF.get_vector(0, 2) + CubeDEF.get_origin(0);
        PosGrid[1] = i * CubeDEF.get_vector(1, 0) + j * CubeDEF.get_vector(1, 1) + k * CubeDEF.get_vector(1, 2) + CubeDEF.get_origin(1);
        PosGrid[2] = i * CubeDEF.get_vector(2, 0) + j * CubeDEF.get_vector(2, 1) + k * CubeDEF.get_vector(2, 2) + CubeDEF.get_origin(2);

        bool skip = true;
        for (int a = 0; a < wavy.get_ncen(); a++)
          if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius / 0.52)
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
    if (i != 0 && i % step == 0 && omp_get_thread_num() == 0)
      progress->write((i + CubeDEF.get_size(0)) / double(3 * CubeDEF.get_size(0)));
  }
  delete(progress);

  time_t end;
  time(&end);
  if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
  else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
  else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Hirshfeld(
  cube& CubeHDEF,
  cube& CubeRho,
  WFN& wavy,
  int cpus,
  double radius,
  int ignore,
  ofstream& file
)
{
#ifdef _OPENMP
  if (cpus != -1) {
    omp_set_num_threads(cpus);
    omp_set_dynamic(0);
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif

  time_t start;
  time(&start);

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating Values" };
  const int step = max(floor(3 * CubeHDEF.get_size(0) / 20), 1.0);

  vector<Thakkar> atoms;
  for (int a = 0; a < wavy.get_ncen(); a++)
    atoms.push_back(Thakkar(wavy.get_atom_charge(a)));

#pragma omp parallel for schedule(dynamic)
  for (int i = -CubeHDEF.get_size(0); i < 2 * CubeHDEF.get_size(0); i++) {
    for (int j = -CubeHDEF.get_size(1); j < 2 * CubeHDEF.get_size(1); j++)
      for (int k = -CubeHDEF.get_size(2); k < 2 * CubeHDEF.get_size(2); k++) {

        double PosGrid[3];

        PosGrid[0] = i * CubeHDEF.get_vector(0, 0) + j * CubeHDEF.get_vector(0, 1) + k * CubeHDEF.get_vector(0, 2) + CubeHDEF.get_origin(0);
        PosGrid[1] = i * CubeHDEF.get_vector(1, 0) + j * CubeHDEF.get_vector(1, 1) + k * CubeHDEF.get_vector(1, 2) + CubeHDEF.get_origin(1);
        PosGrid[2] = i * CubeHDEF.get_vector(2, 0) + j * CubeHDEF.get_vector(2, 1) + k * CubeHDEF.get_vector(2, 2) + CubeHDEF.get_origin(2);

        bool skip = true;
        if (sqrt(pow(PosGrid[0] - wavy.atoms[ignore].x, 2) + pow(PosGrid[1] - wavy.atoms[ignore].y, 2) + pow(PosGrid[2] - wavy.atoms[ignore].z, 2)) < radius / 0.52)
          skip = false;
        if (skip)
          continue;

        double dens_choice = 0.0;
        double dens_all = 0.0;
        double dist, temp;
        for (int a = 0; a < wavy.get_ncen(); a++) {
          dist = sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2));
          temp = atoms[a].get_radial_density(dist);
          if (ignore == a)
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
    if (i != 0 && i % step == 0 && omp_get_thread_num() == 0)
      progress->write((i + CubeHDEF.get_size(0)) / double(3 * CubeHDEF.get_size(0)));
  }
  delete(progress);

  time_t end;
  time(&end);
  if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
  else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
  else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Hirshfeld(
  cube& CubeHDEF,
  cube& CubeRho,
  cube& CubeSpherical,
  WFN& wavy,
  int cpus,
  double radius,
  int ignore,
  ofstream& file
)
{
#ifdef _OPENMP
  if (cpus != -1) {
    omp_set_num_threads(cpus);
    omp_set_dynamic(0);
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif

  time_t start;
  time(&start);

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating Values" };
  const int step = max(floor(3 * CubeHDEF.get_size(0) / 20), 1.0);

  Thakkar atom(wavy.get_atom_charge(ignore));

#pragma omp parallel for schedule(dynamic)
  for (int i = -CubeHDEF.get_size(0); i < 2 * CubeHDEF.get_size(0); i++) {
    for (int j = -CubeHDEF.get_size(1); j < 2 * CubeHDEF.get_size(1); j++)
      for (int k = -CubeHDEF.get_size(2); k < 2 * CubeHDEF.get_size(2); k++) {

        double PosGrid[3];

        PosGrid[0] = i * CubeHDEF.get_vector(0, 0) + j * CubeHDEF.get_vector(0, 1) + k * CubeHDEF.get_vector(0, 2) + CubeHDEF.get_origin(0);
        PosGrid[1] = i * CubeHDEF.get_vector(1, 0) + j * CubeHDEF.get_vector(1, 1) + k * CubeHDEF.get_vector(1, 2) + CubeHDEF.get_origin(1);
        PosGrid[2] = i * CubeHDEF.get_vector(2, 0) + j * CubeHDEF.get_vector(2, 1) + k * CubeHDEF.get_vector(2, 2) + CubeHDEF.get_origin(2);

        bool skip = true;
        double dist = sqrt(pow(PosGrid[0] - wavy.atoms[ignore].x, 2) + pow(PosGrid[1] - wavy.atoms[ignore].y, 2) + pow(PosGrid[2] - wavy.atoms[ignore].z, 2));
        if (dist < radius / 0.52)
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
    if (i != 0 && i % step == 0 && omp_get_thread_num() == 0)
      progress->write((i + CubeHDEF.get_size(0)) / double(3 * CubeHDEF.get_size(0)));
  }
  delete(progress);

  time_t end;
  time(&end);
  if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
  else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
  else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Hirshfeld_atom(
  cube& CubeHirsh,
  cube& CubeRho,
  cube& CubeSpherical,
  WFN& wavy,
  int cpus,
  double radius,
  int ignore,
  ofstream& file
)
{
#ifdef _OPENMP
  if (cpus != -1) {
    omp_set_num_threads(cpus);
    omp_set_dynamic(0);
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif

  time_t start;
  time(&start);

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating Values" };
  const int step = max(floor(3 * CubeHirsh.get_size(0) / 20), 1.0);

  Thakkar atom(wavy.get_atom_charge(ignore));

#pragma omp parallel for schedule(dynamic)
  for (int i = -CubeHirsh.get_size(0); i < 2 * CubeHirsh.get_size(0); i++) {
    for (int j = -CubeHirsh.get_size(1); j < 2 * CubeHirsh.get_size(1); j++)
      for (int k = -CubeHirsh.get_size(2); k < 2 * CubeHirsh.get_size(2); k++) {

        double PosGrid[3];

        PosGrid[0] = i * CubeHirsh.get_vector(0, 0) + j * CubeHirsh.get_vector(0, 1) + k * CubeHirsh.get_vector(0, 2) + CubeHirsh.get_origin(0);
        PosGrid[1] = i * CubeHirsh.get_vector(1, 0) + j * CubeHirsh.get_vector(1, 1) + k * CubeHirsh.get_vector(1, 2) + CubeHirsh.get_origin(1);
        PosGrid[2] = i * CubeHirsh.get_vector(2, 0) + j * CubeHirsh.get_vector(2, 1) + k * CubeHirsh.get_vector(2, 2) + CubeHirsh.get_origin(2);

        bool skip = true;
        double dist = sqrt(pow(PosGrid[0] - wavy.atoms[ignore].x, 2) + pow(PosGrid[1] - wavy.atoms[ignore].y, 2) + pow(PosGrid[2] - wavy.atoms[ignore].z, 2));
        if (dist < radius / 0.52)
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
    if (i != 0 && i % step == 0 && omp_get_thread_num() == 0)
      progress->write((i + CubeHirsh.get_size(0)) / double(3 * CubeHirsh.get_size(0)));
  }
  delete(progress);

  time_t end;
  time(&end);
  if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
  else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
  else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Rho(
  cube& CubeRho,
  WFN& wavy,
  int cpus,
  double radius,
  ofstream& file
)
{
#ifdef _OPENMP
  if (cpus != -1) {
    omp_set_num_threads(cpus);
    omp_set_dynamic(0);
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif

  time_t start;
  time(&start);

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating Values" };
  const int step = max(floor(CubeRho.get_size(0) * 3 / 20), 1.0);

#pragma omp parallel for schedule(dynamic)
  for (int i = -CubeRho.get_size(0); i < 2 * CubeRho.get_size(0); i++) {
    for (int j = -CubeRho.get_size(1); j < 2 * CubeRho.get_size(1); j++)
      for (int k = -CubeRho.get_size(2); k < 2 * CubeRho.get_size(2); k++) {

        double PosGrid[3],
          Rho;

        PosGrid[0] = i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0);
        PosGrid[1] = i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1);
        PosGrid[2] = i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2);

        bool skip = true;
        for (int a = 0; a < wavy.get_ncen(); a++)
          if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius / 0.52)
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
      progress->write((i + CubeRho.get_size(0)) / double(CubeRho.get_size(0) * 3));
  }
  delete(progress);

  time_t end;
  time(&end);
  if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
  else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
  else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Rho_spherical_harmonics(
  cube& CubeRho,
  WFN& wavy,
  int cpus,
  ofstream& file
)
{
#ifdef _OPENMP
  if (cpus != -1) {
    omp_set_num_threads(cpus);
    //omp_set_dynamic(0);
  }
#endif

  time_t start;
  time(&start);

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating Rho" };
  const int step = max(floor(CubeRho.get_size(0) * 3 / 20), 1.0);

//#pragma omp parallel shared(CubeRho)
  {
    cout << omp_get_thread_num() << endl;
    vector<vector<double>> d(5);
    for (int i = 0; i < 5; i++)
      d[i].resize(wavy.get_ncen(), 0.0);
    const int n = wavy.get_nmo(true);
    vector<double> phi(n, 0.0);
//#pragma omp for schedule(dynamic)
    for (int i = 0; i < CubeRho.get_size(0); i++) {
      for (int j = 0; j < CubeRho.get_size(1); j++)
        for (int k = 0; k < CubeRho.get_size(2); k++) {

          double PosGrid[3];

          PosGrid[0] = i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0);
          PosGrid[1] = i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1);
          PosGrid[2] = i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2);

          CubeRho.set_value(i, j, k, wavy.compute_dens(PosGrid[0], PosGrid[1], PosGrid[2], d, phi));
        }
      if (i != 0 && i % step == 0)
        progress->write((i) / double(CubeRho.get_size(0)));
    }
  }
  delete(progress);

  time_t end;
  time(&end);
  if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
  else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
  else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_MO_spherical_harmonics(
  cube& CubeMO,
  WFN& wavy,
  int cpus,
  int MO,
  ofstream& file
)
{
#ifdef _OPENMP
  if (cpus != -1) {
    omp_set_num_threads(cpus);
    omp_set_dynamic(0);
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif

  time_t start;
  time(&start);

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating Values" };
  const int step = max(floor(CubeMO.get_size(0) * 3 / 20), 1.0);

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < CubeMO.get_size(0); i++) {
    for (int j = 0; j < CubeMO.get_size(1); j++)
      for (int k = 0; k < CubeMO.get_size(2); k++) {

        double PosGrid[3];

        PosGrid[0] = i * CubeMO.get_vector(0, 0) + j * CubeMO.get_vector(0, 1) + k * CubeMO.get_vector(0, 2) + CubeMO.get_origin(0);
        PosGrid[1] = i * CubeMO.get_vector(1, 0) + j * CubeMO.get_vector(1, 1) + k * CubeMO.get_vector(1, 2) + CubeMO.get_origin(1);
        PosGrid[2] = i * CubeMO.get_vector(2, 0) + j * CubeMO.get_vector(2, 1) + k * CubeMO.get_vector(2, 2) + CubeMO.get_origin(2);

        CubeMO.set_value(i, j, k, wavy.compute_MO_spherical(PosGrid[0], PosGrid[1], PosGrid[2], MO));
      }
    if (i != 0 && i % step == 0)
      progress->write((i) / double(CubeMO.get_size(0)));
  }
  delete(progress);

  time_t end;
  time(&end);
  if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
  else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
  else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Prop(
  cube& CubeRho,
  cube& CubeRDG,
  cube& CubeElf,
  cube& CubeEli,
  cube& CubeLap,
  cube& CubeESP,
  WFN& wavy,
  int cpus,
  double radius,
  ofstream& file,
  bool test
)
{
#ifdef _OPENMP
  if (cpus != -1) {
    omp_set_num_threads(cpus);
    omp_set_dynamic(0);
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif

  time_t start;
  time(&start);

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating Values" };
  const int step = max(floor(CubeRho.get_size(0) * 3 / 20), 1.0);

#pragma omp parallel for schedule(dynamic)
  for (int i = -CubeRho.get_size(0); i < 2 * CubeRho.get_size(0); i++) {
    for (int j = -CubeRho.get_size(1); j < 2 * CubeRho.get_size(1); j++)
      for (int k = -CubeRho.get_size(2); k < 2 * CubeRho.get_size(2); k++) {

        double PosGrid[3],
          Rho,
          Grad,
          Elf,
          Eli,
          Lap,
          Hess[9];

        PosGrid[0] = i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0);
        PosGrid[1] = i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1);
        PosGrid[2] = i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2);

        bool skip = true;
        for (int a = 0; a < wavy.get_ncen(); a++)
          if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius / 0.52)
            skip = false;
        if (skip)
          continue;

        if (CubeESP.get_loaded() && !CubeRDG.get_loaded())
          Rho = wavy.compute_dens(
            PosGrid[0], PosGrid[1], PosGrid[2]
          );

        if (CubeRDG.get_loaded() && CubeLap.get_loaded() && (CubeElf.get_loaded() || CubeEli.get_loaded()))
          wavy.computeValues(
            PosGrid,
            Rho,
            Grad,
            Hess,
            Elf,
            Eli,
            Lap
          );
        else if ((CubeElf.get_loaded() && CubeEli.get_loaded()) && !CubeRDG.get_loaded() && !CubeLap.get_loaded())
          wavy.computeELIELF(
            PosGrid,
            Elf,
            Eli
          );
        else if (CubeElf.get_loaded() && !CubeEli.get_loaded() && !CubeRDG.get_loaded() && !CubeLap.get_loaded())
          wavy.computeELF(
            PosGrid,
            Elf
          );
        else if (!CubeElf.get_loaded() && CubeEli.get_loaded() && !CubeRDG.get_loaded() && !CubeLap.get_loaded())
          wavy.computeELI(
            PosGrid,
            Eli
          );
        else if (CubeElf.get_loaded() && CubeEli.get_loaded() && CubeLap.get_loaded() && !CubeRDG.get_loaded())
          wavy.computeLapELIELF(
            PosGrid,
            Elf,
            Eli,
            Lap
          );
        else if (!CubeElf.get_loaded() && CubeEli.get_loaded() && CubeLap.get_loaded() && !CubeRDG.get_loaded())
          wavy.computeLapELI(
            PosGrid,
            Eli,
            Lap
          );
        else
          wavy.computeValues(
            PosGrid,
            Rho,
            Grad,
            Hess,
            Elf,
            Eli,
            Lap
          );

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
        if (CubeRDG.get_loaded()) {
          if (isnan(Grad)) Grad = 0;
          if (isinf(Grad)) Grad = 0;
          CubeRDG.set_value(temp_i, temp_j, temp_k, CubeRDG.get_value(temp_i, temp_j, temp_k) + Grad);
        }
        if (CubeLap.get_loaded()) {
          if (isnan(Lap)) Lap = 0;
          if (isinf(Lap)) Lap = 0;
          CubeLap.set_value(temp_i, temp_j, temp_k, CubeLap.get_value(temp_i, temp_j, temp_k) + Lap);
        }
        if (CubeElf.get_loaded()) {
          if (isnan(Elf)) Elf = 0;
          if (isinf(Elf)) Elf = 0;
          CubeElf.set_value(temp_i, temp_j, temp_k, CubeElf.get_value(temp_i, temp_j, temp_k) + Elf);
        }
        if (CubeEli.get_loaded()) {
          if (isnan(Eli)) Eli = 0;
          if (isinf(Eli)) Eli = 0;
          CubeEli.set_value(temp_i, temp_j, temp_k, CubeEli.get_value(temp_i, temp_j, temp_k) + Eli);
        }
      }
    if (i != 0 && i % step == 0)
      progress->write((i + CubeRho.get_size(0)) / double(CubeRho.get_size(0) * 3));
  }
  delete(progress);
  if (!test) {
    time_t end;
    time(&end);
    if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
    else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
    else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
  }
};

void Calc_ESP(
  cube& CubeESP,
  WFN& wavy,
  int cpus,
  double radius,
  ofstream& file
)
{
#ifdef _OPENMP
  if (cpus != -1) {
    omp_set_num_threads(cpus);
    omp_set_dynamic(0);
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif
  time_t start;
  time(&start);
  vector < vector <double> > d2;
  d2.resize(wavy.get_ncen());
  for (int i = 0; i < wavy.get_ncen(); i++) {
    d2[i].resize(wavy.get_ncen());
    for (int j = 0; j < wavy.get_ncen(); j++) {
      if (i == j) {
        d2[i][j] = 0;
        continue;
      }
      d2[i][j] = pow(wavy.atoms[i].x - wavy.atoms[j].x, 2) + pow(wavy.atoms[i].y - wavy.atoms[j].y, 2) + pow(wavy.atoms[i].z - wavy.atoms[j].z, 2);
    }
  }

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating ESP" };
  const int step = max(floor(CubeESP.get_size(0) * 3 / 20), 1.0);

#pragma omp parallel for schedule(dynamic)
  for (int i = -CubeESP.get_size(0); i < 2 * CubeESP.get_size(0); i++) {
    for (int j = -CubeESP.get_size(1); j < 2 * CubeESP.get_size(1); j++)
      for (int k = -CubeESP.get_size(2); k < 2 * CubeESP.get_size(2); k++) {

        double PosGrid[3];

        PosGrid[0] = i * CubeESP.get_vector(0, 0) + j * CubeESP.get_vector(0, 1) + k * CubeESP.get_vector(0, 2) + CubeESP.get_origin(0);
        PosGrid[1] = i * CubeESP.get_vector(1, 0) + j * CubeESP.get_vector(1, 1) + k * CubeESP.get_vector(1, 2) + CubeESP.get_origin(1);
        PosGrid[2] = i * CubeESP.get_vector(2, 0) + j * CubeESP.get_vector(2, 1) + k * CubeESP.get_vector(2, 2) + CubeESP.get_origin(2);

        bool skip = true;
        for (int a = 0; a < wavy.get_ncen(); a++)
          if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius / 0.52)
            skip = false;
        if (skip)
          continue;

        int temp_i, temp_j, temp_k;
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

        CubeESP.set_value(temp_i, temp_j, temp_k, CubeESP.get_value(temp_i, temp_j, temp_k) + wavy.computeESP(PosGrid, d2));
        //CubeESP.set_value(i, j, k, computeESP(PosGrid, d2, wavy));
      }
    if (i != 0 && i % step == 0)
      progress->write((i + CubeESP.get_size(0)) / double(CubeESP.get_size(0) * 3));
  }
  delete(progress);

  time_t end;
  time(&end);
  if (difftime(end, start) < 60) file << "Time to calculate ESP: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
  else if (difftime(end, start) < 3600) file << "Time to calculate ESP: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
  else file << "Time to calculate ESP: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_MO(
  cube& CubeMO,
  int mo,
  WFN& wavy,
  int cpus,
  double radius,
  ofstream& file
)
{
#ifdef _OPENMP
  if (cpus != -1) {
    omp_set_num_threads(cpus);
    omp_set_dynamic(0);
    if (cpus > 1)
      omp_set_nested(1);
  }
#endif
  time_t start;
  time(&start);

  progress_bar* progress = new progress_bar{ file, 50u, "Calculating MO" };
  const int step = max(floor(CubeMO.get_size(0) * 3 / 20), 1.0);

#pragma omp parallel for schedule(dynamic)
  for (int i = -CubeMO.get_size(0); i < 2 * CubeMO.get_size(0); i++) {
    for (int j = -CubeMO.get_size(1); j < 2 * CubeMO.get_size(1); j++)
      for (int k = -CubeMO.get_size(2); k < 2 * CubeMO.get_size(2); k++) {

        double PosGrid[3];
        double MO;

        PosGrid[0] = i * CubeMO.get_vector(0, 0) + j * CubeMO.get_vector(0, 1) + k * CubeMO.get_vector(0, 2) + CubeMO.get_origin(0);
        PosGrid[1] = i * CubeMO.get_vector(1, 0) + j * CubeMO.get_vector(1, 1) + k * CubeMO.get_vector(1, 2) + CubeMO.get_origin(1);
        PosGrid[2] = i * CubeMO.get_vector(2, 0) + j * CubeMO.get_vector(2, 1) + k * CubeMO.get_vector(2, 2) + CubeMO.get_origin(2);

        bool skip = true;
        for (int a = 0; a < wavy.get_ncen(); a++)
          if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius / 0.52)
            skip = false;
        if (skip)
          continue;

        MO = wavy.computeMO(
          PosGrid,
          mo
        );

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
      progress->write((i + CubeMO.get_size(0)) / double(CubeMO.get_size(0) * 3));
  }
  delete(progress);

  time_t end;
  time(&end);
  if (difftime(end, start) < 60) file << "Time to calculate MO: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
  else if (difftime(end, start) < 3600) file << "Time to calculate MO: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
  else file << "Time to calculate MO: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};
//end here