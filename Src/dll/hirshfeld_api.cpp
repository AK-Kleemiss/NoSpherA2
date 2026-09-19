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
    std::cout << "Calculating Hirshfeld Surface for " << fn1 << " and " << fn2 << std::endl;
    WFN wfn1(fn1, false);
    WFN wfn2(fn2, false);
    properties_options opts;
    opts.radius = radius;
    opts.resolution = resolution;
    std::vector<Triangle> triangles_i = Hirshfeld_surface(wfn1, wfn2, opts, std::cout);

    std::array<std::array<int, 3>, 3> Colourcode;
    Colourcode[0] = {255, 0, 0};
    Colourcode[1] = {255, 255, 255};
    Colourcode[2] = {0, 0, 255};

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
