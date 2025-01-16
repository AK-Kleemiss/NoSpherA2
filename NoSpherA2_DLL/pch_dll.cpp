// pch_dll.cpp: source file corresponding to the pre-compiled header

#include "pch_dll.h"

// When you are using pre-compiled headers, this source file is necessary for compilation to succeed.
// This is the only source file that includes pch_dll.h

std::vector<Triangle> compute_Hirshfeld_suface_i(std::filesystem::path& fn1, std::filesystem::path& fn2, double& resolution, double& radius) {

    if (radius < 2.5)
    {
        std::cout << "Resetting Radius to at least 2.5!" << std::endl;
        radius = 2.5;
    }
    WFN wfn1(fn1, false);
    WFN wfn2(fn2, false);
    double MinMax[6];
    int NbSteps[3];
    readxyzMinMax_fromWFN(wfn1, MinMax, NbSteps, radius, resolution, false);
    cube Hirshfeld_grid(NbSteps[0], NbSteps[1], NbSteps[2], wfn1.get_ncen(), true);
    cube Hirshfeld_grid2(NbSteps[0], NbSteps[1], NbSteps[2], wfn2.get_ncen(), true);
    Hirshfeld_grid.give_parent_wfn(wfn1);
    Hirshfeld_grid2.give_parent_wfn(wfn2);
    double len[3]{ 0, 0, 0 };
    for (int i = 0; i < 3; i++)
    {
        len[i] = (MinMax[3 + i] - MinMax[i]) / NbSteps[i];
    }
    for (int i = 0; i < 3; i++)
    {
        Hirshfeld_grid.set_origin(i, MinMax[i]);
        Hirshfeld_grid2.set_origin(i, MinMax[i]);
        Hirshfeld_grid.set_vector(i, i, len[i]);
        Hirshfeld_grid2.set_vector(i, i, len[i]);
    }
    Hirshfeld_grid.set_comment1("Calculated density using NoSpherA2");
    Hirshfeld_grid.set_comment2("from " + wfn1.get_path().string());
    Hirshfeld_grid2.set_comment1("Calculated density using NoSpherA2");
    Hirshfeld_grid2.set_comment2("from " + wfn2.get_path().string());
    Calc_Spherical_Dens(Hirshfeld_grid, wfn1, radius, std::cout, false);
    Calc_Spherical_Dens(Hirshfeld_grid2, wfn2, radius, std::cout, false);
    cube Total_Dens = Hirshfeld_grid + Hirshfeld_grid2;
    Total_Dens.give_parent_wfn(wfn1);
    cube Hirshfeld_weight = Hirshfeld_grid / Total_Dens;
    Hirshfeld_weight.give_parent_wfn(wfn1);
    std::array<std::array<int, 3>, 3> Colourcode;

    Colourcode[0] = { 255, 0, 0 };
    Colourcode[1] = { 255, 255, 255 };
    Colourcode[2] = { 0, 0, 255 };

    std::vector<Triangle> triangles_i = marchingCubes(Hirshfeld_weight, 0.5, 1);
    std::cout << "Found " << triangles_i.size() << " triangles!" << std::endl;
    double area = 0.0;
    double volume = 0.0;
    double low_lim_di = 1E7;
    double high_lim_di = 0.0;
#pragma omp parallel for reduction(+ : area, volume)
    for (int i = 0; i < triangles_i.size(); i++)
    {
        area += triangles_i[i].calc_area();
        volume += triangles_i[i].calc_inner_volume();
        std::array<double, 3> pos = triangles_i[i].calc_center();
        double d_i = calc_d_i(pos, wfn1);
#pragma omp critical
        {
            if (d_i < low_lim_di)
                low_lim_di = d_i;
            if (d_i > high_lim_di)
                high_lim_di = d_i;
        }
    }
    std::cout << "d_i is scaled from " << low_lim_di << " to " << high_lim_di * 0.9 << std::endl;
#pragma omp parallel for
    for (int i = 0; i < triangles_i.size(); i++)
    {
        get_colour(triangles_i[i], calc_d_i, wfn1, Colourcode, low_lim_di, high_lim_di * 0.9);
    }
    std::cout << "Total area: " << area << std::endl;
    std::cout << "Total volume: " << volume << std::endl;
    std::cout << "Finished!" << std::endl;
    return triangles_i;
}