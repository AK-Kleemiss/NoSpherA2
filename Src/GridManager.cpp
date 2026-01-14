#include "pch.h"
#include "GridManager.h"
#include "spherical_density.h"
#include "integrator.h"

double make_sphericals(
    vec2& dens,
    vec2& dist,
    const ivec& atom_type_list,
    std::ostream& file,
    bool debug,
    double incr_start,
    double min_dist,
    int accuracy)
{
    using namespace std;
    const double incr = pow(incr_start, max(1, accuracy - 1));
    const double lincr = log(incr);
    vector<Thakkar> sphericals; // sphericals is a vector that will store objects of type Thakkar.
    for (int i = 0; i < atom_type_list.size(); i++)
        sphericals.emplace_back(atom_type_list[i]); // push_back is like a copy without making a copy
    // Make radial grids
    if (debug)
    {
        file << "\nSize of atom_type_list:" << setw(5) << atom_type_list.size() << endl;
        for (int i = 0; i < atom_type_list.size(); i++)
        {
            file << "\nCalculating for atomic number " << atom_type_list[i] << endl;
            double current = 1;
            double _dist = min_dist;
            if (accuracy > 3)
                while (current > 1E-10)
                {
                    dist[i].push_back(_dist);
                    current = sphericals[i].get_radial_density(_dist);
                    if (current == -20)
                        return -1000;
                    dens[i].push_back(current);
                    _dist *= incr;
                }
            else
                while (current > 1E-12)
                {
                    dist[i].push_back(_dist);
                    current = sphericals[i].get_radial_density(_dist);
                    if (current == -20)
                        return false;
                    dens[i].push_back(current);
                    _dist *= incr;
                }
            file << "Number of radial density points for atomic number " << atom_type_list[i] << ": " << dens[i].size() << endl;
        }
    }
    else
    {
#pragma omp parallel for
        for (int i = 0; i < atom_type_list.size(); i++)
        {
            double current = 1;
            double _dist = min_dist;
            if (accuracy > 3)
                while (current > 1E-10)
                {
                    dist[i].push_back(_dist);
                    current = sphericals[i].get_radial_density(_dist);
                    dens[i].push_back(current);
                    _dist *= incr;
                }
            else
                while (current > 1E-12)
                {
                    dist[i].push_back(_dist);
                    current = sphericals[i].get_radial_density(_dist);
                    dens[i].push_back(current);
                    _dist *= incr;
                }
        }
    }
    sphericals.clear();
    if (debug)
    {
        file << "Cleared the sphericals!" << endl;
    }
    return lincr;
}

// Add this constexpr function before the GridManager class methods
constexpr LebedevGridParams getLebedevGridParams(int accuracy, int atom_type, int max_l_temp) {
    const bool is_hydrogen = (atom_type == 1);
    const bool low_angular_momentum = (max_l_temp < 3);

    switch (accuracy) {
    case 0:
        if (is_hydrogen) {
            return {
                low_angular_momentum ? constants::lebedev_table[0] : constants::lebedev_table[1],
                low_angular_momentum ? constants::lebedev_table[0] : constants::lebedev_table[1],
                1e-3
            };
        }
        else {
            return {
                low_angular_momentum ? constants::lebedev_table[0] : constants::lebedev_table[1],
                low_angular_momentum ? constants::lebedev_table[0] : constants::lebedev_table[1],
                1e-4
            };
        }

    case 1:
        if (is_hydrogen) {
            return {
                low_angular_momentum ? constants::lebedev_table[2] : constants::lebedev_table[3],
                low_angular_momentum ? constants::lebedev_table[6] : constants::lebedev_table[7],
                1e-5
            };
        }
        else {
            return {
                low_angular_momentum ? constants::lebedev_table[3] : constants::lebedev_table[4],
                low_angular_momentum ? constants::lebedev_table[7] : constants::lebedev_table[8],
                1e-4
            };
        }

    case 2:
        if (is_hydrogen) {
            return {
                low_angular_momentum ? constants::lebedev_table[6] : constants::lebedev_table[7],
                low_angular_momentum ? constants::lebedev_table[10] : constants::lebedev_table[11],
                1e-6
            };
        }
        else {
            return {
                low_angular_momentum ? constants::lebedev_table[7] : constants::lebedev_table[8],
                low_angular_momentum ? constants::lebedev_table[11] : constants::lebedev_table[12],
                1e-5
            };
        }

    case 3:
        if (is_hydrogen) {
            return {
                low_angular_momentum ? constants::lebedev_table[11] : constants::lebedev_table[13],
                low_angular_momentum ? constants::lebedev_table[13] : constants::lebedev_table[15],
                1e-11
            };
        }
        else {
            return {
                low_angular_momentum ? constants::lebedev_table[12] : constants::lebedev_table[14],
                low_angular_momentum ? constants::lebedev_table[14] : constants::lebedev_table[16],
                1e-10
            };
        }

    case 4:
        if (is_hydrogen) {
            return {
                low_angular_momentum ? constants::lebedev_table[13] : constants::lebedev_table[16],
                low_angular_momentum ? constants::lebedev_table[18] : constants::lebedev_table[20],
                1e-15
            };
        }
        else {
            return {
                low_angular_momentum ? constants::lebedev_table[14] : constants::lebedev_table[17],
                low_angular_momentum ? constants::lebedev_table[19] : constants::lebedev_table[21],
                1e-20
            };
        }

    default: // accuracy >= 5
        if (is_hydrogen) {
            return {
                low_angular_momentum ? constants::lebedev_table[28] : constants::lebedev_table[30],
                low_angular_momentum ? constants::lebedev_table[30] : constants::lebedev_table[32],
                1e-15
            };
        }
        else {
            return {
                low_angular_momentum ? constants::lebedev_table[29] : constants::lebedev_table[31],
                low_angular_momentum ? constants::lebedev_table[31] : constants::lebedev_table[32],
                1e-20
            };
        }
    }
}

GridManager::GridManager(const GridConfiguration& config)
    : config_(config) {
}

void GridManager::setupGridsForMolecule(const WFN& wave, const bvec& needs_grid,
    const ivec& atom_list, const cell& unit_cell) {
    if (config_.debug) {
        std::cout << "GridManager: Setting up grids for " << atom_list.size()
            << " atoms with " << config_.getPartitionName() << " partitioning" << std::endl;
    }

    // Clear previous data
    grid_data_.clear();
    prototype_grids_.clear();
    atom_type_list_.clear();

    // Identify unique atom types
    atom_type_list_ = identifyAtomTypes(wave, needs_grid);

    // Setup prototype grids for each atom type
    setupPrototypeGrids(wave, atom_type_list_);
    addTimingPoint("Prototype Grid setup");

    // Generate integration grids for each atom
    generateIntegrationGrids(wave, unit_cell, atom_list);
    addTimingPoint("Atomic Grid setup");

    // Calculate Hirshfeld weights if needed
    if (config_.partition_type == PartitionType::Hirshfeld || config_.debug || config_.all_charges) {
        calculateHirshfeldWeights(wave, unit_cell, atom_list);
        addTimingPoint("Hirshfeld Weights");
    }

    // Prune grid based on cutoff criteria
    pruneGrid();
    addTimingPoint("Grid Pruning");

    // Calculate non-spherical densities
    if (config_.debug) {
        std::cout << std::endl
            << "Using " << wave.get_nmo() << " MOs in temporary wavefunction" << std::endl;
        wave.write_wfn("temp_wavefunction.wfn", false, true);
    }
    calculateNonSphericalDensities(wave, unit_cell);
    addTimingPoint("WFN evaluation on grid");

    if (config_.debug) {
        std::cout << "GridManager: Setup complete. Total grid points: "
            << grid_data_.total_points << std::endl;
    }
}

void GridManager::setupPrototypeGrids(const WFN& wave, const ivec& atom_types) {
    std::cout << "GridManager: Setting up prototype grids for atom types..." << std::endl;
    prototype_grids_.clear();
    prototype_grids_.reserve(atom_types.size());

    for (const int atom_type : atom_types) {
        // Calculate basis set parameters for this atom type
        double alpha_max = 0.0;
        int max_l = 0;
        vec alpha_min(6, 1e8);  // Support up to f functions

        for (int i = 0; i < wave.get_ncen(); ++i) {
            if (wave.get_atom_charge(i) != atom_type) continue;

            for (int b = 0; b < wave.get_nex(); ++b) {
                if (wave.get_center(b) != i + 1) continue;

                alpha_max = std::max(alpha_max, wave.get_exponent(b));

                int l = wave.get_type(b);
                if (l == 1) l = 1;
                else if (l >= 2 && l <= 4) l = 2;
                else if (l >= 5 && l <= 10) l = 3;
                else if (l >= 11 && l <= 20) l = 4;
                else if (l >= 21 && l <= 35) l = 5;
                else if (l >= 36 && l <= 56) l = 6;

                max_l = std::max(max_l, l);
                alpha_min[l - 1] = std::min(alpha_min[l - 1], wave.get_exponent(b));
            }
        }

        int max_l_temp = max_l - 1;
        if (config_.debug) {
            std::cout << "Atom Type: " << atom_type << std::endl;
            std::cout << "max_l: " << max_l_temp << " alpha_max: " << std::scientific << alpha_max << " alpha_min: ";
            for (int l = 0; l <= max_l_temp; ++l) {
                std::cout << std::setw(14) << std::scientific << alpha_min[l];
            }
            std::cout << " accuracy: " << config_.accuracy << std::endl;
        }

        err_checkf(config_.accuracy >= 0, "Negative accuracy is not defined!", std::cout);
        // Get Lebedev grid parameters using the constexpr function
        const auto grid_params = getLebedevGridParams(config_.accuracy, atom_type, max_l_temp);

        // Create the prototype grid with the determined parameters
        prototype_grids_.emplace_back(
            grid_params.radial_accuracy,
            grid_params.lebedev_low,
            grid_params.lebedev_high,
            atom_type,
            alpha_max,
            max_l_temp,
            alpha_min.data(),
            std::cout
        );
        if (config_.debug) {
            std::cout << "Created prototype grid for atom type " << atom_type
                << " with " << prototype_grids_.back().get_num_grid_points()
                << " points." << std::endl;
        }
    }
}

void GridManager::generateIntegrationGrids(const WFN& wave, const cell& unit_cell,
    const ivec& atom_list) {
    std::cout << "GridManager: Generating integration grids for atoms..." << std::endl;
    const int num_atoms_with_grids = atom_list.size();
    grid_data_.resizeForAtoms(num_atoms_with_grids);

    // Setup coordinate arrays for all atoms (including PBC images)
    const int total_atoms = wave.get_ncen() * static_cast<int>(std::pow(config_.pbc * 2 + 1, 3));
    vec x_coords(total_atoms), y_coords(total_atoms), z_coords(total_atoms);
    ivec charges(total_atoms);

    // Fill coordinate arrays
    const int pbc = config_.pbc;
#pragma omp parallel for
    for (int i = 0; i < wave.get_ncen(); ++i) {
        charges[i] = wave.get_atom_charge(i);
        x_coords[i] = wave.get_atom_coordinate(i, 0);
        y_coords[i] = wave.get_atom_coordinate(i, 1);
        z_coords[i] = wave.get_atom_coordinate(i, 2);
        if (pbc != 0) {
            int j = 0;
            for (int pbc_x = -pbc; pbc_x < pbc + 1; pbc_x++)
                for (int pbc_y = -pbc; pbc_y < pbc + 1; pbc_y++)
                    for (int pbc_z = -pbc; pbc_z < pbc + 1; pbc_z++)
                    {
                        if (pbc_x == 0 && pbc_y == 0 && pbc_z == 0)continue;
                        j++;
                        charges[i + j * wave.get_ncen()] = wave.get_atom_charge(i);
                        x_coords[i + j * wave.get_ncen()] = wave.get_atom_coordinate(i, 0) +
                            pbc_x * unit_cell.get_cm(0, 0) +
                            pbc_y * unit_cell.get_cm(0, 1) +
                            pbc_z * unit_cell.get_cm(0, 2);
                        y_coords[i + j * wave.get_ncen()] = wave.get_atom_coordinate(i, 1) +
                            pbc_x * unit_cell.get_cm(1, 0) +
                            pbc_y * unit_cell.get_cm(1, 1) +
                            pbc_z * unit_cell.get_cm(1, 2);
                        z_coords[i + j * wave.get_ncen()] = wave.get_atom_coordinate(i, 2) +
                            pbc_x * unit_cell.get_cm(2, 0) +
                            pbc_y * unit_cell.get_cm(2, 1) +
                            pbc_z * unit_cell.get_cm(2, 2);
                    }

        }
    }

    // Generate grids for each atom
    vec2 chi_matrix;  // For TFVC partitioning
    for (int i = 0; i < num_atoms_with_grids; ++i) {
        const int atom_idx = atom_list[i];
        const int atom_type = wave.get_atom_charge(atom_idx);

        // Find corresponding prototype grid
        int prototype_idx = 0;
        for (int j = 0; j < atom_type_list_.size(); ++j) {
            if (atom_type_list_[j] == atom_type) {
                prototype_idx = j;
                break;
            }
        }

        const int num_points = prototype_grids_[prototype_idx].get_num_grid_points();
        grid_data_.num_points_per_atom[i] = num_points;

        // Resize grid arrays for this atom
        for (int coord = 0; coord < 8; ++coord) {
            grid_data_.atomic_grids[i][coord].resize(num_points);
        }

        // Generate the actual grid
        prototype_grids_[prototype_idx].get_grid(
            total_atoms,
            atom_idx,
            x_coords.data(),
            y_coords.data(),
            z_coords.data(),
            charges.data(),
            grid_data_.atomic_grids[i][GridData::GridIndex::X].data(),
            grid_data_.atomic_grids[i][GridData::GridIndex::Y].data(),
            grid_data_.atomic_grids[i][GridData::GridIndex::Z].data(),
            grid_data_.atomic_grids[i][GridData::GridIndex::WEIGHT].data(),
            grid_data_.atomic_grids[i][GridData::GridIndex::BECKE_WEIGHT].data(),
            grid_data_.atomic_grids[i][GridData::GridIndex::TFVC_WEIGHT].data(),
            wave,
            chi_matrix,
            config_.debug
        );
        if (config_.debug) std::cout << "Generated grid for atom " << i + 1 << "/" << num_atoms_with_grids
            << " (Type " << atom_type << ") with " << num_points << " points." << std::endl;
    }

    grid_data_.total_points = std::accumulate(grid_data_.num_points_per_atom.begin(),
        grid_data_.num_points_per_atom.end(), 0);
    std::cout << "GridManager: Generated total of " << grid_data_.total_points << " grid points." << std::endl;
}

PartitionResults GridManager::calculatePartitionedCharges(const WFN& wave, const cell& unit_cell) {
    if (config_.debug) {
        std::string scheme_name = (config_.debug || config_.all_charges) ? "every" : config_.getPartitionName();
        std::cout << "GridManager: Calculating partitioned charges using "
            << scheme_name << " scheme..." << std::endl;
    }
    enum CHARGE_ORDER { S_BECKE = 0, S_TFVC = 1, S_HIRSH = 2 };

    PartitionResults results;

    // Initialize charge arrays
    const int num_atoms = grid_data_.atomic_grids.size();
    results.atom_charges.resize(3);  // Becke, Hirshfeld, TFVC
    for (int i = 0; i < 3; ++i) {
        results.atom_charges[i].resize(num_atoms, 0.0);
    }

    //Choose grids to consider for charge calculation
    GridData::GridIndex weight_index;
    switch (config_.partition_type) {
    case PartitionType::Becke:     weight_index = GridData::GridIndex::BECKE_WEIGHT;  break;
    case PartitionType::TFVC:      weight_index = GridData::GridIndex::TFVC_WEIGHT;   break;
    case PartitionType::Hirshfeld: weight_index = GridData::GridIndex::HIRSH_WEIGHT;  break;
    default:                       weight_index = GridData::GridIndex::BECKE_WEIGHT;  break;
    }

    const double cutoff = config_.getCutoff();
#pragma omp parallel for schedule(static)
    for (int atom = 0; atom < num_atoms; ++atom) {
        vec2& atomic_grid = grid_data_.atomic_grids[atom];
        const int n_points = grid_data_.num_points_per_atom[atom];

        const double* rho = atomic_grid[GridData::GridIndex::WFN_DENSITY].data();
        const double* wB = atomic_grid[GridData::GridIndex::BECKE_WEIGHT].data();
        const double* wH = atomic_grid[GridData::GridIndex::HIRSH_WEIGHT].data();
        const double* wT = atomic_grid[GridData::GridIndex::TFVC_WEIGHT].data();

        if (config_.debug || config_.all_charges) {
            double accB = 0.0, accH = 0.0, accT = 0.0;
            // Vectorizable inner loop; no shared writes.
#pragma omp simd reduction(+:accB, accH, accT)
            for (int p = 0; p < n_points; ++p) {
                const double r = rho[p];
                accB += r * wB[p];
                accH += r * wH[p];
                accT += r * wT[p];
            }

            results.atom_charges[CHARGE_ORDER::S_BECKE][atom] = accB;
            results.atom_charges[CHARGE_ORDER::S_HIRSH][atom] = accH;
            results.atom_charges[CHARGE_ORDER::S_TFVC][atom] = accT;
        }
        else {
            const double* w = atomic_grid[weight_index].data();
            double acc = 0.0;

#pragma omp simd reduction(+:acc)
            for (int p = 0; p < n_points; ++p) {
                acc += rho[p] * w[p];
            }

            // Store only into the active scheme
            switch (config_.partition_type) {
            case PartitionType::Becke:     results.atom_charges[CHARGE_ORDER::S_BECKE][atom] = acc;  break;
            case PartitionType::TFVC:      results.atom_charges[CHARGE_ORDER::S_TFVC][atom] = acc;  break;
            case PartitionType::Hirshfeld: results.atom_charges[CHARGE_ORDER::S_HIRSH][atom] = acc;  break;
            }
        }

        // Add ECP electrons (only to computed schemes)
        if (wave.get_has_ECPs()) {
            const int ecp_e = wave.get_atom_ECP_electrons(atom); // map atom->basis if needed
            if (config_.debug || config_.all_charges) {
                results.atom_charges[S_BECKE][atom] += ecp_e;
                results.atom_charges[S_HIRSH][atom] += ecp_e;
                results.atom_charges[S_TFVC][atom] += ecp_e;
            }
            else {
                switch (config_.partition_type) {
                case PartitionType::Becke:     results.atom_charges[S_BECKE][atom] += ecp_e; break;
                case PartitionType::TFVC:      results.atom_charges[S_TFVC][atom] += ecp_e; break;
                case PartitionType::Hirshfeld: results.atom_charges[S_HIRSH][atom] += ecp_e; break;
                }
            }
        }
    }

    // Calculate total charges
    results.overall_charges.resize(3, 0.0);
    for (int scheme = 0; scheme < 3; ++scheme) {
        results.overall_charges[scheme] = std::accumulate(
            results.atom_charges[scheme].begin(),
            results.atom_charges[scheme].end(), 0.0);
    }

    return results;
}

void GridManager::getDensityVectors(const WFN& wave, const ivec& atom_list, vec2& d1, vec2& d2, vec2& d3, vec2& dens) {
    if (config_.debug) {
        std::cout << "GridManager: Generating density vectors..." << std::endl;
    }
    const int n_atoms = grid_data_.atomic_grids.size();
    const double cutoff = config_.getCutoff();
    //Choose grids to consider for charge calculation
    GridData::GridIndex idx_single; // default
    switch (config_.partition_type) {
    case PartitionType::Becke:     idx_single = GridData::GridIndex::BECKE_WEIGHT;  break;
    case PartitionType::TFVC:      idx_single = GridData::GridIndex::TFVC_WEIGHT;   break;
    case PartitionType::Hirshfeld: idx_single = GridData::GridIndex::HIRSH_WEIGHT;  break;
    default:                        std::cout << "GridManager: Unknown partition type for density vectors!" << std::endl; exit(1); return;
    }

    d1.resize(n_atoms); d2.resize(n_atoms); d3.resize(n_atoms); dens.resize(n_atoms);
#pragma omp parallel for schedule(dynamic, 1)
    for (int g = 0; g < grid_data_.atomic_grids.size(); ++g) {
        vec2& atomic_grid = grid_data_.atomic_grids[g];
        const int n_points = grid_data_.num_points_per_atom[g];
        dens[g].resize(n_points); d1[g].resize(n_points); d2[g].resize(n_points); d3[g].resize(n_points);


        const double* rho = atomic_grid[GridData::GridIndex::WFN_DENSITY].data();
        const double* w = atomic_grid[idx_single].data();
        double* res = dens[g].data();


        double* d1_ptr = d1[g].data();
        double* d2_ptr = d2[g].data();
        double* d3_ptr = d3[g].data();
        const double* x = atomic_grid[GridData::GridIndex::X].data();
        const double* y = atomic_grid[GridData::GridIndex::Y].data();
        const double* z = atomic_grid[GridData::GridIndex::Z].data();
        const double x0 = wave.get_atom_coordinate(atom_list[g], 0);
        const double y0 = wave.get_atom_coordinate(atom_list[g], 1);
        const double z0 = wave.get_atom_coordinate(atom_list[g], 2);

        size_t accepted_points = 0;
        for (int p = 0; p < n_points; ++p) {
            const double t = rho[p] * w[p];
            const int keep = std::fabs(t) >= cutoff;  // 0/1
            res[accepted_points] = t;
            d1_ptr[accepted_points] = x[p] - x0; d2_ptr[accepted_points] = y[p] - y0; d3_ptr[accepted_points] = z[p] - z0;
            accepted_points += keep;
        }

        //I really like this one... but it seems to be slightly slower than the above version
        //#pragma omp simd
        //        for (int p = 0; p < n_points; ++p) {
        //            res[p] = rho[p] * w[p];
        //        }
        //        for (int p = 0; p < n_points; ++p) {
        //            const int keep = std::fabs(res[p]) > cutoff;
        //            const int idx = accepted_points;
        //            
        //            res[idx] = res[p];
        //            d1_ptr[idx] = x[p] - x0;
        //            d2_ptr[idx] = y[p] - y0;
        //            d3_ptr[idx] = z[p] - z0;
        //            accepted_points += keep;
        //        }

        grid_data_.num_points_per_atom[g] = accepted_points;
        d1[g].resize(accepted_points);
        d2[g].resize(accepted_points);
        d3[g].resize(accepted_points);
        dens[g].resize(accepted_points);
    }
    //Calculate final point count
    grid_data_.total_points = std::accumulate(grid_data_.num_points_per_atom.begin(),
        grid_data_.num_points_per_atom.end(), 0);
}

// Implementation of other methods...

double GridConfiguration::getCutoff() const {
    if (accuracy < 3) return 1e-10;
    else if (accuracy == 3) return 1e-14;
    else return 1e-30;
}

std::string GridConfiguration::getPartitionName() const {
    switch (partition_type) {
    case PartitionType::Hirshfeld: return "Hirshfeld";
    case PartitionType::Becke: return "Becke";
    case PartitionType::TFVC: return "TFVC";
    case PartitionType::RI: return "RI";
    default: return "Unknown";
    }
}

void GridData::clear() {
    atomic_grids.clear();
    num_points_per_atom.clear();
    total_points = 0;
}

void GridData::resizeForAtoms(int num_atoms) {
    atomic_grids.resize(num_atoms);
    num_points_per_atom.resize(num_atoms);

    for (int i = 0; i < num_atoms; ++i) {
        atomic_grids[i].resize(8);  // x, y, z, weight, hirsh_w, becke_w, tfvc_w, wfn_density
    }
}

ivec GridManager::identifyAtomTypes(const WFN& wave, const bvec& needs_grid) {
    std::set<int> unique_types;

    for (int i = 0; i < wave.get_ncen(); ++i) {
        if (i < needs_grid.size() && (needs_grid[i])) {
            const int charge = wave.get_atom_charge(i);
            if (charge != 119) {  // Skip dummy atoms
                unique_types.insert(charge);
            }
        }
    }

    return ivec(unique_types.begin(), unique_types.end());
}

void PartitionResults::printChargeTable(const svec& labels, const WFN& wave, const ivec atom_list, std::ostream& file) const {
    file << "Table of Charges in electrons\n"
        << "    Atom       Becke   Hirshfeld   TFVC\n";

    for (int i = 0; i < atom_list.size(); ++i) {
        const int atom_idx = atom_list[i];
        file << std::setw(10) << labels[i]
            << std::fixed << std::setw(10) << std::setprecision(3)
            << wave.get_atom_charge(atom_idx) - atom_charges[0][i]  // Becke
            << std::fixed << std::setw(10) << std::setprecision(3)
            << wave.get_atom_charge(atom_idx) - atom_charges[2][i]  // Hirshfeld
            << std::fixed << std::setw(10) << std::setprecision(3)
            << wave.get_atom_charge(atom_idx) - atom_charges[1][i]  // TFVC
            << std::endl;
    }

    file << "Total number of electrons:\n"
        << " Becke:     " << std::fixed << std::setw(10) << std::setprecision(6) << overall_charges[0] << "\n"
        << " Hirshfeld: " << std::fixed << std::setw(10) << std::setprecision(6) << overall_charges[2] << "\n"
        << " TVFC:      " << std::fixed << std::setw(10) << std::setprecision(6) << overall_charges[1] << std::endl;
}

void GridManager::calculateHirshfeldWeights(const WFN& wave, const cell& unit_cell, const ivec& atom_list) {
    using namespace std;
    bvec fake_needs_grid(wave.get_ncen(), true);
    ivec complete_type_list = identifyAtomTypes(wave, fake_needs_grid);

    if (config_.debug) {
        std::cout << "GridManager: Calculating spherical densities..." << std::endl;
    }

    // Setup spherical density calculation
    radial_density_.resize(complete_type_list.size());
    radial_distances_.resize(complete_type_list.size());

    for (int i = 0; i < complete_type_list.size(); ++i) {
        radial_density_[i].reserve(1000);
        radial_distances_[i].reserve(1000);
    }

    // Calculate radial densities using Thakkar spherical atoms
    lincr_ = make_sphericals(radial_density_, radial_distances_, complete_type_list,
        std::cout, config_.debug, 1.005, 1.0E-7, config_.accuracy);
    err_checkf(lincr_ != -1000, "error during creations of sphericals", std::cout);

    vec2 single_spherical_density(grid_data_.atomic_grids.size());
    vec2 combined_spherical_density(grid_data_.atomic_grids.size());

    for (int g = 0; g < grid_data_.atomic_grids.size(); ++g) {
        const int num_points = grid_data_.num_points_per_atom[g];
        single_spherical_density[g].resize(num_points, 0.0);
        combined_spherical_density[g].resize(num_points, 0.0);
    }


    // For each atom
    for (int atom_idx = 0; atom_idx < wave.get_ncen(); atom_idx++) {
        // Find the type index for this atom
        int type_idx = -1;
        for (int j = 0; j < complete_type_list.size(); ++j) {
            if (wave.get_atom_charge(atom_idx) == complete_type_list[j]) {
                type_idx = j;
                break;
            }
        }

        if (type_idx == -1) continue; // Skip if atom type not found

        const d3 ax = wave.get_atom_pos(atom_idx);

        // Add this atom's spherical density contribution to all grids
        for (int g = 0; g < grid_data_.atomic_grids.size(); ++g) {
            int g_to_atom_idx = atom_list[g];
            const int num_points = grid_data_.num_points_per_atom[g];

#pragma omp parallel for
            for (int p = 0; p < num_points; ++p) {

                const double dist = array_length(d3{ grid_data_.atomic_grids[g][GridData::GridIndex::X][p],
                grid_data_.atomic_grids[g][GridData::GridIndex::Y][p],
                grid_data_.atomic_grids[g][GridData::GridIndex::Z][p] }, ax);

                const double density = linear_interpolate_spherical_density(
                    radial_density_[type_idx],
                    radial_distances_[type_idx],
                    dist,
                    lincr_,
                    1.0E-7);

                if (atom_idx == g_to_atom_idx) {
                    single_spherical_density[g][p] += density;
                }
                combined_spherical_density[g][p] += density;
            }
        }
    }

    // Add PBC contributions if enabled
    if (config_.pbc != 0) {
        if (config_.debug) {
            std::cout << "GridManager: Adding PBC contributions (pbc=" << config_.pbc << ")" << std::endl;
        }

        for (int pbc_x = -config_.pbc; pbc_x <= config_.pbc; ++pbc_x) {
            for (int pbc_y = -config_.pbc; pbc_y <= config_.pbc; ++pbc_y) {
                for (int pbc_z = -config_.pbc; pbc_z <= config_.pbc; ++pbc_z) {
                    // Skip the original unit cell (already processed above)
                    if (pbc_x == 0 && pbc_y == 0 && pbc_z == 0) continue;

                    // For each atom in the original unit cell
                    for (int atom_idx = 0; atom_idx < wave.get_ncen(); ++atom_idx) {
                        // Find the type index for this atom
                        int type_idx = -1;
                        for (int j = 0; j < complete_type_list.size(); ++j) {
                            if (wave.get_atom_charge(atom_idx) == complete_type_list[j]) {
                                type_idx = j;
                                break;
                            }
                        }

                        if (type_idx == -1) continue;

                        // Calculate PBC image coordinates
                        const d3 a = { wave.get_atom_coordinate(atom_idx, 0) +
                            pbc_x * unit_cell.get_cm(0, 0) +
                            pbc_y * unit_cell.get_cm(0, 1) +
                            pbc_z * unit_cell.get_cm(0, 2),
                          wave.get_atom_coordinate(atom_idx, 1) +
                            pbc_x * unit_cell.get_cm(1, 0) +
                            pbc_y * unit_cell.get_cm(1, 1) +
                            pbc_z * unit_cell.get_cm(1, 2),
                          wave.get_atom_coordinate(atom_idx, 2) +
                            pbc_x * unit_cell.get_cm(2, 0) +
                            pbc_y * unit_cell.get_cm(2, 1) +
                            pbc_z * unit_cell.get_cm(2, 2) };

                        // Add contributions to all grids
                        for (int g = 0; g < grid_data_.atomic_grids.size(); ++g) {
                            const int num_points = grid_data_.num_points_per_atom[g];

#pragma omp parallel for
                            for (int p = 0; p < num_points; ++p) {
                                const d3 g_pos = { grid_data_.atomic_grids[g][0][p],
                                 grid_data_.atomic_grids[g][1][p],
                                 grid_data_.atomic_grids[g][2][p] };

                                const double dist = array_length(g_pos, a);

                                // Use the SAME interpolation function
                                const double density = linear_interpolate_spherical_density(
                                    radial_density_[type_idx],
                                    radial_distances_[type_idx],
                                    dist,
                                    lincr_,
                                    1.0E-7);

                                // Only add to total density for PBC (not to atom-specific density)
                                combined_spherical_density[g][p] += density;
                            }
                        }
                    }
                }
            }
        }
    }

    if (config_.debug) {
        std::cout << "GridManager: Spherical densities calculated, calculating Hirshfeld weights..." << std::endl;
    }
    //Calculate hirshfeld weights
    for (int g = 0; g < grid_data_.atomic_grids.size(); ++g) {
        double* combined_ptr = combined_spherical_density[g].data();
        double* single_ptr = single_spherical_density[g].data();
        vec2& atom_grid = grid_data_.atomic_grids[g];
        for (int p = 0; p < grid_data_.num_points_per_atom[g]; ++p) {
            atom_grid[GridData::GridIndex::HIRSH_WEIGHT][p] =
                (combined_ptr[p] != 0.0) ? (atom_grid[GridData::GridIndex::WEIGHT][p] * (single_ptr[p] / combined_ptr[p])) : 0.0;
        }
    }
}

void GridManager::calculateNonSphericalDensities(const WFN& wave, const cell& unit_cell) {
    using namespace std;

    if (config_.debug) {
        std::cout << "GridManager: Calculating non-spherical densities..." << std::endl;
    }

#pragma omp parallel
    {
        vec2 d_temp(16);
        for (int i = 0; i < 16; i++)
        {
            d_temp[i].resize(wave.get_ncen(), 0.0);
        }
        vec phi_temp(wave.get_nmo(true), 0.0);
        for (int g = 0; g < grid_data_.atomic_grids.size(); ++g) {
            const int num_points = grid_data_.num_points_per_atom[g];
            vec2& atom_grid = grid_data_.atomic_grids[g];
            double* x_ptr = atom_grid[GridData::GridIndex::X].data(); double* y_ptr = atom_grid[GridData::GridIndex::Y].data(); double* z_ptr = atom_grid[GridData::GridIndex::Z].data();
            double* densy_ptr = atom_grid[GridData::GridIndex::WFN_DENSITY].data();
#pragma omp for
            for (int p = 0; p < num_points; ++p) {
                // Calculate WFN density at this point
                densy_ptr[p] = wave.compute_dens({ x_ptr[p], y_ptr[p], z_ptr[p] }, d_temp, phi_temp);
            }
        }
    }
}

void GridManager::pruneGrid() {
    if (config_.debug) {
        std::cout << "GridManager: Pruning grid..." << std::endl;
    }

    const double cutoff = config_.getCutoff();
    const int original_size = grid_data_.total_points;

    //Choose grids to consider for pruning
    GridData::GridIndex weight_index;
    switch (config_.partition_type) {
    case PartitionType::Becke:      weight_index = GridData::GridIndex::BECKE_WEIGHT;  break;
    case PartitionType::TFVC:       weight_index = GridData::GridIndex::TFVC_WEIGHT;   break;
    case PartitionType::Hirshfeld:  weight_index = GridData::GridIndex::HIRSH_WEIGHT;  break;
    default:                        weight_index = GridData::GridIndex::BECKE_WEIGHT;  break;
    }

    bvec2 point_is_kept(grid_data_.atomic_grids.size());
    ivec pruned_num_points(grid_data_.atomic_grids.size());

    // Parallelize safely: each g writes a distinct element.
#pragma omp parallel for
    for (int g = 0; g < static_cast<int>(grid_data_.atomic_grids.size()); ++g) {
        point_is_kept[g].resize(grid_data_.num_points_per_atom[g], false);
        bvec& kept_local = point_is_kept[g]; // For easier access in lambda;

        const int n = grid_data_.num_points_per_atom[g];
        const auto& grid = grid_data_.atomic_grids[g];

        if (config_.debug || config_.all_charges) {
            const double* w_Becke = grid[GridData::GridIndex::BECKE_WEIGHT].data();
            const double* w_TFVC = grid[GridData::GridIndex::TFVC_WEIGHT].data();
            const double* w_Hirsh = grid[GridData::GridIndex::HIRSH_WEIGHT].data();
#pragma omp simd
            for (int p = 0; p < n; ++p) {
                // OR-reduce comparisons; branchless.
                bool keep = fabs(w_Becke[p]) > cutoff;
                keep |= fabs(w_TFVC[p]) > cutoff;
                keep |= fabs(w_Hirsh[p]) > cutoff;

                kept_local[p] = keep;
            }
        }
        else {
            const double* weights = grid[weight_index].data();
#pragma omp simd
            for (int p = 0; p < n; ++p) {
                kept_local[p] = (fabs(weights[p]) > cutoff);//keep;
            }
        }
        pruned_num_points[g] = std::count(kept_local.begin(), kept_local.end(), true);
    }

    // Create pruned grids
    vec3 pruned_atomic_grids(grid_data_.atomic_grids.size());
#pragma omp parallel for
    for (int g = 0; g < grid_data_.atomic_grids.size(); ++g) {
        pruned_atomic_grids[g].resize(8);
        for (int coord = 0; coord < 8; ++coord) {
            pruned_atomic_grids[g][coord].resize(pruned_num_points[g]);
        }
        const bvec& kept_local = point_is_kept[g];

        const vec2& source = grid_data_.atomic_grids[g];
        vec2& dest = pruned_atomic_grids[g];

        int reduced_point = 0;
        for (int p = 0; p < grid_data_.num_points_per_atom[g]; ++p) {
            if (!kept_local[p]) continue;

            dest[GridData::GridIndex::X][reduced_point] = source[GridData::GridIndex::X][p];
            dest[GridData::GridIndex::Y][reduced_point] = source[GridData::GridIndex::Y][p];
            dest[GridData::GridIndex::Z][reduced_point] = source[GridData::GridIndex::Z][p];
            dest[GridData::GridIndex::WEIGHT][reduced_point] = source[GridData::GridIndex::WEIGHT][p];
            dest[GridData::GridIndex::HIRSH_WEIGHT][reduced_point] = source[GridData::GridIndex::HIRSH_WEIGHT][p];
            dest[GridData::GridIndex::BECKE_WEIGHT][reduced_point] = source[GridData::GridIndex::BECKE_WEIGHT][p];
            dest[GridData::GridIndex::TFVC_WEIGHT][reduced_point] = source[GridData::GridIndex::TFVC_WEIGHT][p];

            reduced_point++;
        }
    }

    // Update grid data
    grid_data_.atomic_grids = std::move(pruned_atomic_grids);
    grid_data_.num_points_per_atom = std::move(pruned_num_points);
    grid_data_.total_points = std::accumulate(grid_data_.num_points_per_atom.begin(),
        grid_data_.num_points_per_atom.end(), 0);

    if (config_.debug) {
        std::cout << "GridManager: Grid pruned from " << original_size << " to " << grid_data_.total_points << " points ("
            << std::defaultfloat << (100.0 * grid_data_.total_points / original_size) << "%)" << std::endl;
    }
}

bvec GridManager::determineAtomsNeedingGrids(const WFN& wave, const ivec& asym_atom_list) {
    bvec needs_grid(wave.get_ncen(), false);

    // Mark atoms in the asymmetric unit list as needing grids
    for (int atom_idx : asym_atom_list) {
        if (atom_idx >= 0 && atom_idx < wave.get_ncen()) {
            // Skip dummy atoms (charge 119)
            if (wave.get_atom_charge(atom_idx) != 119) {
                needs_grid[atom_idx] = true;
            }
        }
    }

    return needs_grid;
}
