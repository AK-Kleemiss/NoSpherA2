#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "constants.h"
#include "b2c.h"

using namespace std;

bool b2c(const cube *cub, const vector<atom> &atoms, bool debug, bool bcp)
{
    iMatrix3 CP(cub->get_size(0), cub->get_size(1), cub->get_size(2));
    ivec2 Liste(3);
    double GradMax, xlength, ylength, zlength;
    xlength = std::sqrt(
        cub->get_vector(0, 0) * cub->get_vector(0, 0) +
        cub->get_vector(1, 0) * cub->get_vector(1, 0) +
        cub->get_vector(2, 0) * cub->get_vector(2, 0));

    ylength = std::sqrt(
        cub->get_vector(0, 1) * cub->get_vector(0, 1) +
        cub->get_vector(1, 1) * cub->get_vector(1, 1) +
        cub->get_vector(2, 1) * cub->get_vector(2, 1));

    zlength = std::sqrt(
        cub->get_vector(0, 2) * cub->get_vector(0, 2) +
        cub->get_vector(1, 2) * cub->get_vector(1, 2) +
        cub->get_vector(2, 2) * cub->get_vector(2, 2));
    if (debug)std::cout << "calculated lengths!" << endl;
    string s, s1, s2, s3;
    double xmin = cub->get_origin(0);//, xstep;
    double ymin = cub->get_origin(1);//, ystep;
    double zmin = cub->get_origin(2);//, zstep;

    if (debug)std::cout << "Resized CP and initialized it; GOING TO MAKE LISTE NOW" << endl;
    int iCP = 0, ListeMax;
    cubepoint Max = { 0,0,0,0.0 };
    vec Maxima;
    vec3 distances;
    svec Labels;
    ivec nrs;
    vector<vector<bvec>> border;
    ivec2 neighbours;
    vector<vector<cubepoint> > BCPs;
    distances.resize(3);
    for (int i = 0; i < 3; i++) {
        distances[i].resize(3);
        for (int j = 0; j < 3; j++) {
            distances[i][j].resize(3);
            for (int k = 0; k < 3; k++) distances[i][j][k] = sqrt(
                pow((1 - i) * xlength, 2)
                + pow((1 - j) * ylength, 2)
                + pow((1 - k) * zlength, 2));
        }
    }
    if (debug)std::cout << "calculated distances!" << endl;
    for (int x = 0; x < cub->get_size(0); x++)
        for (int y = 0; y < cub->get_size(1); y++)
            for (int z = 0; z < cub->get_size(2); z++) {
                if (CP(x, y, z) == 0) {
                    ListeMax = 0;
                    Max.x = x; Max.y = y; Max.z = z;
                    int xs, ys, zs;
                    do {
                        ListeMax++;
                        xs = Max.x, ys = Max.y, zs = Max.z;
                        if (xs < cub->get_size(0) && ys < cub->get_size(1) && zs < cub->get_size(2))
                            Max.value = (cub->get_value(xs, ys, zs));
                        else
                            Max.value = 0.0;
                        if (Liste[0].size() <= ListeMax)
                            Liste[0].push_back(Max.x);  //only append List when needed check if neccesarry
                        else
                            Liste[0][ListeMax - 1] = Max.x;
                        if (Liste[1].size() <= ListeMax)
                            Liste[1].push_back(Max.y);
                        else
                            Liste[1][ListeMax - 1] = Max.y;
                        if (Liste[2].size() <= ListeMax)
                            Liste[2].push_back(Max.z);
                        else
                            Liste[2][ListeMax - 1] = Max.z;
                        GradMax = 0;
                        for (int ix = xs - 1; ix < xs + 2; ix++)
                            for (int iy = ys - 1; iy < ys + 2; iy++)
                                for (int iz = zs - 1; iz < zs + 2; iz++) {
                                    if (ix == xs && iy == ys && iz == zs) continue;
                                    if (ix < 0 || iy < 0 || iz < 0) continue;
                                    if (ix >= cub->get_size(0) || iy >= cub->get_size(1) || iz >= cub->get_size(2)) continue;
                                    /*double dist=sqrt(
                                             ((ix-xs)*xlength)*((ix-xs)*xlength)
                                            +((iy-ys)*ylength)*((iy-ys)*ylength)
                                            +((iz-zs)*zlength)*((iz-zs)*zlength));*/
                                            //    Tests if this voxel has maximum value compared to neighboring voxels
                                    if (((cub->get_value(ix, iy, iz) - (cub->get_value(xs, ys, zs))) / distances[1 + ix - xs][1 + iy - ys][1 + iz - zs]) > GradMax) {
                                        Max.x = ix; Max.y = iy; Max.z = iz;
                                        Max.value = cub->get_value(ix, iy, iz);
                                        GradMax = ((cub->get_value(ix, iy, iz)) - (cub->get_value(xs, ys, zs))) / distances[1 + ix - xs][1 + iy - ys][1 + iz - zs];
                                    }
                                }
                    } while (!(GradMax == 0 || CP(Max.x, Max.y, Max.z) > 0));
                    if (CP(Max.x, Max.y, Max.z) > 0)
                        for (int i = 0; i < ListeMax; i++)
                            CP(Liste[0][i], Liste[1][i], Liste[2][i]) = CP(Max.x, Max.y, Max.z);
                    else {
                        const d3 pos = { Max.x * xlength + cub->get_origin(0),Max.y * ylength + cub->get_origin(1),Max.z * zlength + cub->get_origin(2) };
                        if (debug) {
                            std::cout << "DBUG: Position of CP: ";
                            for (int i = 0; i < 3; i++)std::cout << pos[i] << " ";
                            std::cout << endl;
                        }
                        double min_dist, temp;
                        unsigned int atom = 0;
                        Maxima.push_back(Max.value);
                        iCP++;
                        min_dist = 10000;
                        for (int i = 0; i < ListeMax; i++)
                            CP(Liste[0][i], Liste[1][i], Liste[2][i]) = iCP;
                        for (int i = 0; i < atoms.size(); i++) {
                            temp = array_length(pos, atoms[i].get_pos());
                            //if(debug)std::cout << "DBUG: Distance to atom " << i  << " with label " << atoms[i].get_label() << " is " << temp;
                            if (min_dist > temp) {
                                min_dist = temp;
                                atom = i;
                                //if(debug)std::cout << " *";
                            }
                            //if(debug)std::cout << endl;
                        }
                        Labels.push_back(atoms[atom].get_label());
                        nrs.push_back(atom);
                    }
                    for (int i = 0; i < 3; i++)
                        Liste[i].resize(0); //Revert Liste and free memory
                }
            }
    if (bcp) {
        neighbours.resize(iCP);
        BCPs.resize(iCP);
        border.resize(cub->get_size(0));
        for (int i = 0; i < cub->get_size(0); i++) {
            border[i].resize(cub->get_size(1));
            for (int j = 0; j < cub->get_size(1); j++) {
                border[i][j].resize(cub->get_size(2));
                for (int k = 0; k < cub->get_size(2); k++)
                    border[i][j][k] = false;
            }
        }
        for (int x = 0; x < cub->get_size(0); x++)
            for (int y = 0; y < cub->get_size(1); y++)
                for (int z = 0; z < cub->get_size(2); z++)
                    for (int ix = x - 1; ix < x + 2; ix++)
                        for (int iy = y - 1; iy < y + 2; iy++)
                            for (int iz = z - 1; iz < z + 2; iz++) {
                                if (ix == x && iy == y && iz == z) continue;
                                if (ix <= 0 || iy <= 0 || iz <= 0) continue;
                                if (ix >= cub->get_size(0) || iy >= cub->get_size(1) || iz >= cub->get_size(2)) continue;
                                //    Tests if this voxel has neighboring voxels, that are from a different basin and assigns list of neighbors for each basin
                                if (CP(ix, iy, iz) != CP(x, y, z)) {
                                    border[x][y][z] = true;
                                    int found = false;
                                    for (int i = 0; i < neighbours[CP(x, y, z)].size(); i++)
                                        if (neighbours[CP(x, y, z)][i] == CP(ix, iy, iz))
                                            found = true;
                                    if (!found)
                                        neighbours[CP(x, y, z)].push_back(CP(ix, iy, iz));
                                }
                            }
        //sanity check, all basins must be neighboring each other pairwise!
        for (int b = 0; b < iCP; b++)
            for (int n = 0; n < neighbours[b].size(); n++)
                if (find(neighbours[n].begin(), neighbours[n].end(), b) == neighbours[n].end()) {
                    std::cout << "ERROR: Basins should be neighbours in a pairwise way! Basin " << b << " has neighbour " << n << ", but it does not appear to be the case the other way around!" << endl;
                    return false;
                }
        //possibly better if we assigned max and min values for coordinates of basins to reduce amount of calculations, lets see...
        for (int b = 0; b < iCP; b++) {
            Max = { 0,0,0,0.0 };
            for (int n = 0; n < neighbours[b].size(); n++) {
                for (int x = 0; x < cub->get_size(0); x++)
                    for (int y = 0; y < cub->get_size(1); y++)
                        for (int z = 0; z < cub->get_size(2); z++) {
                            //check that we are on A) A border B) in the basin we want C) the neighboring basin we look at is at this border D) the value is bigger than the previously found maximum
                            if (!border[x][y][z])
                                continue;
                            if (CP(x, y, z) != b)
                                continue;
                            bool found = false;
                            for (int ix = x - 1; ix < x + 2; ix++)
                                for (int iy = y - 1; iy < y + 2; iy++)
                                    for (int iz = z - 1; iz < z + 2; iz++) {
                                        if (ix == x && iy == y && iz == z) continue;
                                        if (ix <= 0 || iy <= 0 || iz <= 0) continue;
                                        if (ix >= cub->get_size(0) || iy >= cub->get_size(1) || iz >= cub->get_size(2)) continue;
                                        if (CP(ix, iy, iz) != n) continue;
                                        else found = true;
                                    }
                            if (found && cub->get_value(x, y, z) > Max.value)
                                Max = { x,y,z,cub->get_value(x,y,z) };
                        }
                BCPs[b].push_back(Max);
            }
        }
        for (int b = 0; b < iCP; b++)
            for (int n = 0; n < neighbours[b].size(); n++) {
                int back_reference = 0;
                for (int i = 0; i < neighbours[neighbours[b][n]].size(); i++)
                    if (neighbours[neighbours[b][n]][i] == b)
                        back_reference = i;
                if (BCPs[b][n].x - BCPs[neighbours[b][n]][back_reference].x > 2 || BCPs[b][n].y - BCPs[neighbours[b][n]][back_reference].y > 2 || BCPs[b][n].z - BCPs[neighbours[b][n]][back_reference].z > 2) {
                    std::cout << "The BCPs on both sides of the ZFS are too far apart..." << endl;
                }
                else {
                    std::cout << "BCP: " << Labels[b] << "-" << Labels[neighbours[b][n]] << " ED: " << (BCPs[b][n].value + BCPs[neighbours[b][n]][back_reference].value) / 2 << " Position: "
                        << (BCPs[b][n].x + BCPs[neighbours[b][n]][back_reference].x) / 2 * cub->get_vector(0, 0)
                        + (BCPs[b][n].y + BCPs[neighbours[b][n]][back_reference].y) / 2 * cub->get_vector(1, 0)
                        + (BCPs[b][n].z + BCPs[neighbours[b][n]][back_reference].z) / 2 * cub->get_vector(2, 0)
                        + xmin << " "
                        << (BCPs[b][n].x + BCPs[neighbours[b][n]][back_reference].x) / 2 * cub->get_vector(0, 1)
                        + (BCPs[b][n].y + BCPs[neighbours[b][n]][back_reference].y) / 2 * cub->get_vector(1, 1)
                        + (BCPs[b][n].z + BCPs[neighbours[b][n]][back_reference].z) / 2 * cub->get_vector(2, 1)
                        + ymin
                        << (BCPs[b][n].x + BCPs[neighbours[b][n]][back_reference].x) / 2 * cub->get_vector(0, 2)
                        + (BCPs[b][n].y + BCPs[neighbours[b][n]][back_reference].y) / 2 * cub->get_vector(1, 2)
                        + (BCPs[b][n].z + BCPs[neighbours[b][n]][back_reference].z) / 2 * cub->get_vector(2, 2)
                        + zmin
                        << endl;
                }
            }
        if (debug)std::cout << "done with BCPs" << endl;
    }
    if (debug)
        std::cout << "done with liste, writing basins now!" << endl;
    ivec basins;
    vec EDS(iCP);
    vec VOL(iCP);
    double dv = cub->get_dv();
    if (debug)std::cout << "dv: " << dv << " iCP: " << iCP << endl;
    for (int x = 0; x < cub->get_size(0); x++)
        for (int y = 0; y < cub->get_size(1); y++)
            for (int z = 0; z < cub->get_size(2); z++)
                for (int a = 0; a < iCP; a++) if (a + 1 == CP(x, y, z)) {
                    EDS[a] += cub->get_value(x, y, z) * dv;            //    Electrons in this voxel
                    VOL[a] += dv;                                    //    Size of the voxel
                }
    if (debug) for (int a = 0; a < iCP; a++)std::cout << a << " eds: " << EDS[a] << " vol: " << VOL[a] << endl;
    unsigned int in = 9999999;
    std::cout << "I found " << iCP << " Basins." << endl;
    for (int a = 0; a < iCP; a++)std::cout << "Integrated value in Basin " << toString<int>(a + 1) << ": " << scientific << setw(14) << setprecision(7) << EDS[a]
        << " Volume: " << VOL[a] << " Maximum: " << Maxima[a] << " possibly atom: " << Labels[a] << "_" << nrs[a] << endl;
    if (debug) {
        std::cout << "DEBUG: Labels_size: " << Labels.size() << endl;
        for (int i = 0; i < Labels.size(); i++)std::cout << Labels[i] << " ";
        std::cout << endl;
    }

    std::cout << "Which of these do you want to include into one cube file?" << endl
        << "Put 0 for end of input. If none are selected (first number is zero), all basins will be written into separate files" << endl;
    while (in != 0) {
        cin >> in;
        if (in != 0) basins.push_back(in);
    }
    string temp;
    string replace("E");
    string path_temp = cub->get_path().replace_extension(".b2c_log").generic_string();
    path_temp.erase(path_temp.find(".cub"), 5);
    ofstream logfile(path_temp.c_str(), ios::out);
    logfile << "Number of Basins: " << iCP << endl;
    double Integral = 0.0;
    temp = cub->get_path().generic_string();
    temp.erase(temp.find(".cub"), 5);
    if (basins.size() > 0) {
        temp += "_" + toString<int>((int)basins.size()) + "_basins";
        temp += ".cube";
        ofstream outfile(temp.c_str(), ios::out);
        outfile << "First comment line?" << endl << "second comment line" << endl << "   " << atoms.size();
        outfile << fixed;
        outfile << setprecision(6);
        outfile << setw(12);
        outfile << xmin << " " << ymin << " " << zmin << endl;
        for (int i = 0; i < 3; i++) {
            outfile << setprecision(0);
            outfile << setw(6);
            outfile << cub->get_size(i) << " ";
            outfile << setw(12);
            outfile << setprecision(6);
            for (int j = 0; j < 3; j++) outfile << cub->get_vector(i, j) << " ";
            outfile << endl;
            outfile.flush();
        }
        for (int i = 0; i < atoms.size(); i++) {
            outfile << "  " << atoms[i].get_charge();
            outfile << setw(12) << setprecision(6);
            outfile << atoms[i].get_charge() << " ";
            outfile << atoms[i].get_coordinate(0) << " " << atoms[i].get_coordinate(1) << " " << atoms[i].get_coordinate(2) << " ";
            outfile << endl;
        }
        outfile.flush();
        stringstream stream;
        for (int x = 0; x < cub->get_size(0); x++) {
            for (int y = 0; y < cub->get_size(1); y++) {
                unsigned int r = 0;
                for (int z = 0; z < cub->get_size(2); z++) {
                    bool include = false;
                    for (int a = 0; a < basins.size(); a++)
                        if (basins[a] == CP(x, y, z))
                            include = true;
                    stream << uppercase << scientific << setw(14) << setprecision(7) << cub->get_value(x, y, z) * include;
                    outfile << stream.str();
                    stream.str("");
                    r++;
                    if (r % 6 == 0) outfile << endl;
                }
                if (r % 6 != 0) outfile << endl;
            }
            outfile.flush();
        }
        for (int a = 0; a < basins.size(); a++) {
            logfile << "Integrated value in Basin " << toString<int>(basins[a]) << ": " << scientific << setw(14) << setprecision(7) << EDS[basins[a] - 1]
                << "Volume: " << VOL[basins[a] - 1] << " Atom: " << Labels[basins[a] - 1] << "_" << nrs[basins[a] - 1] << endl;
            Integral += EDS[basins[a] - 1];
        }
        outfile.close();
    }
    else {
        for (int f = 0; f < iCP; f++) {
            temp = cub->get_path().generic_string() + '_' + Labels[f] + '_' + toString<unsigned int>(nrs[f]) + '_' + toString<int>(f) + ".cube";
            temp.erase(temp.find(".cub"), 5);
            ofstream outfile(temp.c_str(), ios::out);
            outfile << s1 << "\n" << s2 << "\n" << setw(5) << atoms.size();
            outfile << fixed << setprecision(6) << setw(12) << xmin << " " << ymin << " " << zmin << "\n";
            for (int i = 0; i < 3; i++)
                outfile << setprecision(0) << setw(5) << cub->get_size(i) << setw(12) << setprecision(6) << cub->get_vector(i, 0) << setw(12) << setprecision(6) << cub->get_vector(i, 1) << setw(12) << setprecision(6) << cub->get_vector(i, 2) << "\n";
            for (int i = 0; i < atoms.size(); i++) {
                outfile << setw(5) << atoms[i].get_charge() << setw(5) << atoms[i].get_charge() << ".000000";
                outfile << fixed << setw(12) << setprecision(6) << atoms[i].get_coordinate(0);
                outfile << fixed << setw(12) << setprecision(6) << atoms[i].get_coordinate(1);
                outfile << fixed << setw(12) << setprecision(6) << atoms[i].get_coordinate(2);
                outfile << "\n";
            }
            for (int x = 0; x < cub->get_size(0); x++) {
                for (int y = 0; y < cub->get_size(1); y++) {
                    unsigned int r = 0;
                    for (int z = 0; z < cub->get_size(2); z++) {
                        outfile << uppercase << scientific << setw(13) << setprecision(5) << cub->get_value(x, y, z) * (CP(x, y, z) == (f + 1));
                        r++;
                        if (r % 6 == 0) outfile << "\n";
                    }
                    if (r % 6 != 0) outfile << "\n";
                }
            }
            outfile.flush();
            logfile << "Integrated value in Basin " << toString<int>(f + 1) << ": " << scientific << setw(14) << setprecision(7) << EDS[f] << "Volume: " << VOL[f]
                << " Atom: " << Labels[f] << "_" << nrs[f] << endl;
            outfile.close();
            Integral += EDS[f];
        }
    }
    logfile << "Integral over all basins: " << Integral << endl;
    if (bcp)
        for (int b = 0; b < iCP; b++)
            for (int n = 0; n < neighbours[b].size(); n++) {
                int back_reference = 0;
                for (int i = 0; i < neighbours[neighbours[b][n]].size(); i++)
                    if (neighbours[neighbours[b][n]][i] == b)
                        back_reference = i;
                auto pos = cub->get_pos(
                    (BCPs[b][n].x + BCPs[neighbours[b][n]][back_reference].x) / 2,
                    (BCPs[b][n].y + BCPs[neighbours[b][n]][back_reference].y) / 2,
                    (BCPs[b][n].z + BCPs[neighbours[b][n]][back_reference].z) / 2
                );
                logfile << "BCP: " << Labels[b] << "-" << Labels[neighbours[b][n]] << " ED: " << (BCPs[b][n].value + BCPs[neighbours[b][n]][back_reference].value) / 2 << " Position: "
                    << pos[0] << " "
                    << pos[1] << " "
                    << pos[2] << endl;
            }
    logfile.flush();
    logfile.close();
    return true;
};

// Structure to store pre-computed gradient direction for each grid point
struct GradientDirection {
    int next_x, next_y, next_z;  // Next point to follow (-1 if local maximum or invalid)
    double gradient;              // Gradient magnitude
};

std::pair<cubei, std::vector<d4>> topological_cube_analysis(const cube *cub, const vector<atom> &atoms, bool debug, bool bcp, double value_floor, double grad_epsilon, double assignment_radius)
{
    cubei basin_cube({ cub->get_size(0), cub->get_size(1), cub->get_size(2) }, 0, true);
    double xlength, ylength, zlength;
    xlength = std::sqrt(
        cub->get_vector(0, 0) * cub->get_vector(0, 0) +
        cub->get_vector(1, 0) * cub->get_vector(1, 0) +
        cub->get_vector(2, 0) * cub->get_vector(2, 0));

    ylength = std::sqrt(
        cub->get_vector(0, 1) * cub->get_vector(0, 1) +
        cub->get_vector(1, 1) * cub->get_vector(1, 1) +
        cub->get_vector(2, 1) * cub->get_vector(2, 1));

    zlength = std::sqrt(
        cub->get_vector(0, 2) * cub->get_vector(0, 2) +
        cub->get_vector(1, 2) * cub->get_vector(1, 2) +
        cub->get_vector(2, 2) * cub->get_vector(2, 2));

    if (debug)
        std::cout << "calculated lengths!" << endl;

    // Cache grid sizes
    const int size_x = cub->get_size(0);
    const int size_y = cub->get_size(1);
    const int size_z = cub->get_size(2);

    // Pre-compute distance lookup table
    vec3 distances;
    distances.resize(3);
    for (int i = 0; i < 3; i++) {
        distances[i].resize(3);
        for (int j = 0; j < 3; j++) {
            distances[i][j].resize(3);
            for (int k = 0; k < 3; k++) distances[i][j][k] = sqrt(
                pow((1 - i) * xlength, 2)
                + pow((1 - j) * ylength, 2)
                + pow((1 - k) * zlength, 2));
        }
    }

    const double gradient_threshold = std::max(0.0, grad_epsilon);
    const double assignment_radius_bohr = assignment_radius > 0.0 ? constants::ang2bohr(assignment_radius) : -1.0;
    const double assignment_radius2 = assignment_radius_bohr > 0.0 ? assignment_radius_bohr * assignment_radius_bohr : -1.0;
    const bool use_radius_mask = assignment_radius2 > 0.0;

    // Pre-compute radius mask for all grid points (if needed)
    vector<vector<vector<int>>> radius_mask;
    if (use_radius_mask) {
        radius_mask.resize(size_x);
        for (int x = 0; x < size_x; x++) {
            radius_mask[x].resize(size_y);
            for (int y = 0; y < size_y; y++) {
                radius_mask[x][y].resize(size_z, 0);
            }
        }

#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
        for (int x = 0; x < size_x; x++) {
            for (int y = 0; y < size_y; y++) {
                for (int z = 0; z < size_z; z++) {
                    const d3 pos = cub->get_pos(x, y, z);
                    int in_radius = 0;
                    for (const atom &a : atoms) {
                        const d3 apos = a.get_pos();
                        const double dx = pos[0] - apos[0];
                        const double dy = pos[1] - apos[1];
                        const double dz = pos[2] - apos[2];
                        if (dx * dx + dy * dy + dz * dz <= assignment_radius2) {
                            in_radius = 1;
                            break;
                        }
                    }
                    radius_mask[x][y][z] = in_radius;
                }
            }
        }
    }
    std::cout << "PHASE 1: Computing gradient directions..." << endl;

    // Phase 1: Parallel computation of gradient directions
    vector<vector<vector<GradientDirection>>> gradient_field(size_x);
    for (int x = 0; x < size_x; x++) {
        gradient_field[x].resize(size_y);
        for (int y = 0; y < size_y; y++) {
            gradient_field[x][y].resize(size_z);
        }
    }

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int x = 0; x < size_x; x++) {
        for (int y = 0; y < size_y; y++) {
            int *orml = radius_mask[x][y].data(); // Cache radius mask row for performance
            for (int z = 0; z < size_z; z++) {
                GradientDirection &gd = gradient_field[x][y][z];
                gd.next_x = -1;
                gd.next_y = -1;
                gd.next_z = -1;
                gd.gradient = 0.0;

                // Check if point is valid
                const double center_value = cub->get_value(x, y, z);
                if (center_value <= value_floor)
                    continue;
                if (use_radius_mask && orml[z] == 0)
                    continue;

                // Find steepest gradient direction
                double max_grad = 0.0;
                int best_ix = -1, best_iy = -1, best_iz = -1;

                for (int ix = x - 1; ix <= x + 1; ix++) {
                    if (ix < 0 || ix >= size_x) continue;
                    for (int iy = y - 1; iy <= y + 1; iy++) {
                        if (iy < 0 || iy >= size_y) continue;
                        int *rml = radius_mask[ix][iy].data(); // Cache radius mask row for performance
                        for (int iz = z - 1; iz <= z + 1; iz++) {
                            if (ix == x && iy == y && iz == z) continue;
                            if (iz < 0 || iz >= size_z) continue;

                            const double neighbor_value = cub->get_value(ix, iy, iz);
                            if (neighbor_value <= value_floor)
                                continue;
                            if (use_radius_mask && rml[iz] == 0)
                                continue;

                            const double grad = (neighbor_value - center_value) / distances[1 + ix - x][1 + iy - y][1 + iz - z];
                            if (grad > max_grad) {
                                max_grad = grad;
                                best_ix = ix;
                                best_iy = iy;
                                best_iz = iz;
                            }
                        }
                    }
                }

                if (max_grad > gradient_threshold) {
                    gd.next_x = best_ix;
                    gd.next_y = best_iy;
                    gd.next_z = best_iz;
                    gd.gradient = max_grad;
                }
            }
        }
    }
    
    std::cout << "PHASE 2: Following paths to assign basins..." << endl;

    // Phase 2: Sequential basin assignment by following pre-computed paths
    int iCP = 0;
    std::vector<d4> Maxima;
    ivec2 Liste(3);

    for (int x = 0; x < size_x; x++) {
        for (int y = 0; y < size_y; y++) {
            for (int z = 0; z < size_z; z++) {
                if (basin_cube.get_value(x, y, z) != 0)
                    continue;

                // Check if starting point is valid
                const double start_value = cub->get_value(x, y, z);
                if (start_value <= value_floor)
                    continue;
                if (use_radius_mask && !radius_mask[x][y][z])
                    continue;

                // Follow gradient path
                Liste[0].clear();
                Liste[1].clear();
                Liste[2].clear();

                int cx = x, cy = y, cz = z;
                int path_length = 0;
                const int max_path_length = size_x * size_y * size_z; // Safety limit

                while (path_length < max_path_length) {
                    Liste[0].push_back(cx);
                    Liste[1].push_back(cy);
                    Liste[2].push_back(cz);
                    path_length++;

                    // Check if we reached an already-assigned basin
                    int current_basin = basin_cube.get_value(cx, cy, cz);
                    if (current_basin > 0) {
                        // Assign entire path to this basin
                        for (int i = 0; i < path_length; i++) {
                            basin_cube.set_value(Liste[0][i], Liste[1][i], Liste[2][i], current_basin);
                        }
                        break;
                    }

                    // Check if we reached a local maximum
                    const GradientDirection &gd = gradient_field[cx][cy][cz];
                    if (gd.next_x == -1) {
                        // Local maximum - create new basin
                        iCP++;
                        const d3 pos = cub->get_pos(cx, cy, cz);
                        const double max_value = cub->get_value(cx, cy, cz);
                        Maxima.push_back(d4{pos[0], pos[1], pos[2], max_value});

                        if (debug) {
                            std::cout << "DBUG: Position of CP: " << fixed << pos[0] << " " << pos[1] << " " << pos[2] << endl;
                        }

                        // Assign entire path to new basin
                        for (int i = 0; i < path_length; i++) {
                            basin_cube.set_value(Liste[0][i], Liste[1][i], Liste[2][i], iCP);
                        }
                        break;
                    }

                    // Move to next point
                    cx = gd.next_x;
                    cy = gd.next_y;
                    cz = gd.next_z;
                }
            }
        }
    }

    if (debug)
        std::cout << "done with basin assignment!" << endl;
    std::cout << "I found " << iCP << " Basins." << endl;

    // Post-processing: Merge adjacent basins with very similar maxima values
    // This handles cases where high-resolution grids split a single basin across multiple local maxima
    {
        const double merge_value_tolerance = 2e-1;  // Merge if max values differ by less than this
        const double merge_distance_threshold = 27.5; // Merge if basin centers are closer than this (in grid units)
        
        std::vector<int> basin_mapping(iCP + 1);
        for (int i = 0; i <= iCP; i++) basin_mapping[i] = i; // Identity mapping initially
        
        std::cout << "MERGE PHASE: Consolidating adjacent basins with similar maxima..." << endl;
        
        // Find adjacent basins and merge candidates
        for (int b1 = 1; b1 <= iCP; b1++) {
            if (basin_mapping[b1] != b1) continue; // Already merged into another basin
            
            for (int b2 = b1 + 1; b2 <= iCP; b2++) {
                if (basin_mapping[b2] != b2) continue; // Already merged into another basin
                
                const double val1 = Maxima[b1 - 1][3];
                const double val2 = Maxima[b2 - 1][3];
                
                // Check if maxima values are similar enough
                if (std::abs(val1 - val2) > merge_value_tolerance)
                    continue;
                
                // Check if basin centers are close enough
                const double dx = Maxima[b1 - 1][0] - Maxima[b2 - 1][0];
                const double dy = Maxima[b1 - 1][1] - Maxima[b2 - 1][1];
                const double dz = Maxima[b1 - 1][2] - Maxima[b2 - 1][2];
                const double dist = std::sqrt(dx * dx + dy * dy + dz * dz) / xlength;
                
                if (dist <= merge_distance_threshold) {
                    // Merge b2 into b1 (keep larger basin ID to maintain basin indices)
                    int keep_basin = std::max(b1, b2);
                    int merge_basin = std::min(b1, b2);
                    basin_mapping[merge_basin] = keep_basin;
                    
                    if (debug)
                        std::cout << "Merging basin " << merge_basin << " into basin " << keep_basin 
                                  << " (distance: " << dist << " grid units, value diff: " 
                                  << std::abs(val1 - val2) << ")" << endl;
                }
            }
        }
        
        // Apply basin mapping to grid
        int merged_count = 0;
        for (int x = 0; x < size_x; x++) {
            for (int y = 0; y < size_y; y++) {
                for (int z = 0; z < size_z; z++) {
                    int current_basin = basin_cube.get_value(x, y, z);
                    if (current_basin > 0 && current_basin <= iCP) {
                        int target_basin = basin_mapping[current_basin];
                        if (target_basin != current_basin) {
                            basin_cube.set_value(x, y, z, target_basin);
                            merged_count++;
                        }
                    }
                }
            }
        }
        
        if (merged_count > 0) {
            std::cout << "Merged basins: Reassigned " << merged_count << " grid points" << endl;
        }
        
        // Compact basin IDs if some basins were completely merged away
        std::vector<int> basin_active(iCP + 1, 0);
        for (int x = 0; x < size_x; x++) {
            for (int y = 0; y < size_y; y++) {
                for (int z = 0; z < size_z; z++) {
                    int basin = basin_cube.get_value(x, y, z);
                    if (basin > 0 && basin <= iCP)
                        basin_active[basin] = 1;
                }
            }
        }
        
        // Create remapping for compact IDs
        int new_id = 0;
        std::vector<int> compact_map(iCP + 1);
        std::vector<d4> new_Maxima;
        for (int i = 1; i <= iCP; i++) {
            if (basin_active[i]) {
                new_id++;
                compact_map[i] = new_id;
                new_Maxima.push_back(Maxima[i - 1]);
            }
        }
        
        if (new_id < iCP) {
            // Reassign with compact IDs
            for (int x = 0; x < size_x; x++) {
                for (int y = 0; y < size_y; y++) {
                    for (int z = 0; z < size_z; z++) {
                        int basin = basin_cube.get_value(x, y, z);
                        if (basin > 0)
                            basin_cube.set_value(x, y, z, compact_map[basin]);
                    }
                }
            }
            iCP = new_id;
            Maxima = new_Maxima;
            std::cout << "Basin compaction: Reduced from " << basin_mapping.size() - 1 << " to " << iCP << " basins" << endl;
        }
    }

    if (debug)
        std::cout << "done with basin consolidation!" << endl;
    std::cout << "Final basin count: " << iCP << " Basins." << endl;

    for (int d = 0; d < 3; d++) {
        basin_cube.set_origin(d, cub->get_origin(d));
        for (int j = 0; j < 3; j++)
            basin_cube.set_vector(d, j, cub->get_vector(d, j));
    }
    basin_cube.set_path(cub->get_path());
    basin_cube.set_comment1("Topological basin map");
    basin_cube.set_comment2("Value at each grid point = basin index");

    if (debug)
        std::cout << "I found " << iCP << " basins. Returning basin-index cubei." << endl;

    return std::make_pair(basin_cube, Maxima);
};

vec integrate_values_in_basins(const cube *cub, const cubei *basin_cube, svec& basin_label, bool debug)
{
    const int basin_count = basin_cube->max_value();
    vec EDS(basin_count);
    vec VOL(basin_count);
    double dv = cub->get_dv();
    err_checkf(basin_cube->get_size(0) == cub->get_size(0) && basin_cube->get_size(1) == cub->get_size(1) && basin_cube->get_size(2) == cub->get_size(2), "Basin cube and original cube must have the same dimensions!", std::cout);
    err_checkf(basin_cube->get_dv() - cub->get_dv() < 1E-10, "Basin must have reasonable size!", std::cout);
    if (debug)
        std::cout << "dv: " << dv << " iCP: " << basin_count << endl;

#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
    for (int x = 0; x < cub->get_size(0); x++)
        for (int y = 0; y < cub->get_size(1); y++)
            for (int z = 0; z < cub->get_size(2); z++) {
                const int basin_index = basin_cube->get_value(x, y, z) - 1;
                if (basin_index < 0 || basin_index >= basin_count)
                    continue;
                const double contribution = cub->get_value(x, y, z) * dv;
#ifdef _OPENMP
#pragma omp atomic
#endif
                EDS[basin_index] += contribution; // Value in this voxel
#ifdef _OPENMP
#pragma omp atomic
                VOL[basin_index] += dv; // Volume of this voxel
#endif
            }
    for (int a = 0; a < basin_count; a++) {
        if (EDS[a] > 0.001)
            std::cout << "basin: " << setw(4) << a << " label: " << setw(16) << basin_label[a] << " integrated value: " << setw(8) << std::setprecision(4) << std::fixed << EDS[a] << " volume:" << setw(10) << VOL[a] << endl;
    }
    for (int a = 0; a < basin_count; a++) {
        if (EDS[a] < 0.001)
            std::cout << "WARNING: Integrated value in basin " << a + 1 << " is very low (" << EDS[a] << ") and might be inaccurate due to numerical errors!" << endl;
    }
    return EDS;
};

svec assign_labels_to_basins(const vector<d4> &Maxima, const vector<atom> &atoms, bool debug, int type_switch)
{
    svec result(Maxima.size());
    //type switch defines the field that was analyzed and will define how we assign labels to the basins
    switch (type_switch) {
        case 0: //QTAIM case
        for (size_t i = 0; i < Maxima.size(); i++) {
            const d3 pos = {Maxima[i][0], Maxima[i][1], Maxima[i][2]};
            double min_dist = std::numeric_limits<double>::max();
            int atom_index = -1;
            for (size_t j = 0; j < atoms.size(); j++) {
                const d3 apos = atoms[j].get_pos();
                const double dx = pos[0] - apos[0];
                const double dy = pos[1] - apos[1];
                const double dz = pos[2] - apos[2];
                const double dist2 = dx * dx + dy * dy + dz * dz;
                if (dist2 < min_dist) {
                    min_dist = dist2;
                    atom_index = j;
                }
            }
            err_checkf(atom_index >= 0, "No atom found for basin " + toString<size_t>(i) + " at position (" + toString<double>(pos[0]) + ", " + toString<double>(pos[1]) + ", " + toString<double>(pos[2]) + ")!", std::cout);
            result[i] = atoms[atom_index].get_label() + to_string(atom_index);
        }
        break;
        case 1: //ELI case, assign core basins based on atom within, valence basins as the connected basin and bonds as the two closest atoms.
            for (size_t i = 0; i < Maxima.size(); i++) {
                const d3 pos = {Maxima[i][0], Maxima[i][1], Maxima[i][2]};
                double min_dist1 = std::numeric_limits<double>::max();
                double min_dist2 = std::numeric_limits<double>::max();
                int atom_index1 = -1;
                int atom_index2 = -1;
                for (size_t j = 0; j < atoms.size(); j++) {
                    const d3 apos = atoms[j].get_pos();
                    const double dx = pos[0] - apos[0];
                    const double dy = pos[1] - apos[1];
                    const double dz = pos[2] - apos[2];
                    const double dist2 = dx * dx + dy * dy + dz * dz;
                    if (dist2 < min_dist1) {
                        min_dist2 = min_dist1;
                        atom_index2 = atom_index1;
                        min_dist1 = dist2;
                        atom_index1 = j;
                    }
                    else if (dist2 < min_dist2) {
                        min_dist2 = dist2;
                        atom_index2 = j;
                    }
                }
                const double core_dist = atoms[atom_index1].get_charge() < 18 ? 0.125 : atoms[atom_index1].get_charge() < 36 ? 0.35 : 0.7;
                err_checkf(atom_index1 >= 0, "No atom found for basin " + toString<size_t>(i) + " at position (" + toString<double>(pos[0]) + ", " + toString<double>(pos[1]) + ", " + toString<double>(pos[2]) + ")!", std::cout);
                err_checkf(atom_index2 >= 0, "Only one atom found for basin " + toString<size_t>(i) + " at position (" + toString<double>(pos[0]) + ", " + toString<double>(pos[1]) + ", " + toString<double>(pos[2]) + ")!", std::cout);
                double ratio = std::max(1e-5, min_dist1) / std::max(1e-5, min_dist2);
                if (min_dist1 < core_dist && atoms[atom_index1].get_charge() > 2) // If the maximum is very close to an atom, we assume it's a core basin and label it with that atom
                    result[i] = atoms[atom_index1].get_label() + to_string(atom_index1) + " core";
                else if ((ratio < 0.333 || ratio > 3) && atoms[atom_index1].get_charge() > 2) // If the maximum is significantly closer to one atom than to the other, we assume it's a valence basin and label it with the closest atom
                    result[i] = atoms[atom_index1].get_label() + to_string(atom_index1) + " LP";
                else // Otherwise, we assume it's a bond basin and label it with both atoms
                    result[i] = atoms[atom_index1].get_label() + to_string(atom_index1) + "-" + atoms[atom_index2].get_label() + to_string(atom_index2) + " bond";
            }
            break;
        default:
            err_not_impl_f("Label assignment type " + toString<int>(type_switch) + " is not implemented!", std::cout);
    }
    return result;
}
