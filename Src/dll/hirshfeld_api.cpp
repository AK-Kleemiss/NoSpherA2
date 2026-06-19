#include "pch.h"

#include "hirshfeld_api.h"

#include "../core/cube.h"
#include "../core/properties.h"

NOS_API std::vector<Triangle> NOS_CALLCONV compute_Hirshfeld_suface_i(
    const std::filesystem::path& fn1,
    const std::filesystem::path& fn2,
    double resolution,
    double radius)
{
    if (radius < 2.5) {
        std::cout << "Resetting Radius to at least 2.5!" << std::endl;
        radius = 2.5;
    }

    std::cout << "Calculating Hirshfeld Surface for " << fn1 << " and " << fn2 << std::endl;

    // NOTE: This mirrors the legacy implementation exactly, including the hard-coded
    // test file paths. If you want this to use fn1/fn2 directly, say so and I will
    // update it.
    WFN wfn1("D:\\git\\NoSpherA2\\tests\\isosurface\\asu.xyz", true);
    WFN wfn2("D:\\git\\NoSpherA2\\tests\\isosurface\\pack.xyz", true);

    properties_options opts;
    opts.radius = radius;
    opts.resolution = resolution;

    readxyzMinMax_fromWFN(wfn1, opts, false);

    cube Hirshfeld_grid(opts.NbSteps, wfn1.get_ncen(), true);
    cube Hirshfeld_grid2(opts.NbSteps, wfn2.get_ncen(), true);

    Hirshfeld_grid.give_parent_wfn(wfn1);
    Hirshfeld_grid2.give_parent_wfn(wfn2);

    double len[3]{0, 0, 0};
    for (int i = 0; i < 3; i++) {
        len[i] = (opts.MinMax[3 + i] - opts.MinMax[i]) / opts.NbSteps[i];
    }

    for (int i = 0; i < 3; i++) {
        Hirshfeld_grid.set_origin(i, opts.MinMax[i]);
        Hirshfeld_grid2.set_origin(i, opts.MinMax[i]);
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
    Colourcode[0] = {255, 0, 0};
    Colourcode[1] = {255, 255, 255};
    Colourcode[2] = {0, 0, 255};

    std::vector<Triangle> triangles_i = marchingCubes(Hirshfeld_weight, 0.5);

    std::cout << "Found " << triangles_i.size() << " triangles!" << std::endl;

    double area = 0.0;
    double volume = 0.0;
    double low_lim_di = 1E7;
    double high_lim_di = 0.0;

#pragma omp parallel for reduction(+ : area, volume)
    for (int i = 0; i < static_cast<int>(triangles_i.size()); i++) {
        area += triangles_i[i].calc_area();
        volume += triangles_i[i].calc_inner_volume();
        d3 pos = triangles_i[i].calc_center();
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
    for (int i = 0; i < static_cast<int>(triangles_i.size()); i++) {
        get_colour(triangles_i[i], calc_d_i, wfn1, Colourcode, low_lim_di, high_lim_di * 0.9);
    }

    std::cout << "Total area: " << area << std::endl;
    std::cout << "Total volume: " << volume << std::endl;
    std::cout << "Finished!" << std::endl;

    return triangles_i;
}
