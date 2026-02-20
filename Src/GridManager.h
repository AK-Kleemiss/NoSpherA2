#pragma once

#include "convenience.h"
#include "wfn_class.h"
#include "AtomGrid.h"
#include "cell.h"
template<typename AtomType>
double make_sphericals(
    vec2& dens,
    vec2& dist,
    const ivec& atom_type_list,
    std::ostream& file,
    const std::vector<std::pair<vec, vec>>& pop_sig,
    bool debug = false,
    double incr_start = 1.005,
    double min_dist = 0.0000001,
    int accuracy = 2);


// Add this structure to hold grid parameters
struct LebedevGridParams {
    int lebedev_low;
    int lebedev_high;
    double radial_accuracy;
};

// Add this constexpr function declaration in the header
constexpr LebedevGridParams getLebedevGridParams(int accuracy, int atom_type, int max_l_temp);


struct GridConfiguration {
    int accuracy = 2;
    int pbc = 0;
    PartitionType partition_type = PartitionType::Hirshfeld;
    bool debug = false;
    bool all_charges = false;

    double getCutoff() const;
    std::string getPartitionName() const;
};

struct GridData {
    const int grid_data_size = 10;
    enum GridIndex { X = 0, Y = 1, Z = 2, WEIGHT = 3, HIRSH_WEIGHT = 4, BECKE_WEIGHT = 5, TFVC_WEIGHT = 6, WFN_DENSITY = 7, MBIS_WEIGHT = 8, EMBIS_WEIGHT = 9 };
    vec3 atomic_grids;           // [atom][coord_type][point]
    ivec num_points_per_atom;    // Number of points for each atom
    int total_points = 0;

    void clear();
    void resizeForAtoms(int num_atoms);
};

struct PartitionResults {
    vec2 atom_charges;           // [atom][partition type]
    vec overall_charges;      // Total charges by partition type
    enum CHARGE_ORDER { S_BECKE = 0, S_TFVC = 1, S_HIRSH = 2, S_MBIS = 3, S_EMBIS = 4 };
    double nrmsd = 0.0;
    double r_value = 0.0;

    void printChargeTable(const svec& labels, const WFN& wave, const ivec atom_list, std::ostream& file) const;
};

class GridManager {
private:
    GridConfiguration config_;
    std::vector<AtomGrid> prototype_grids_;
    GridData grid_data_;
    vec2 radial_density_;
    vec2 radial_distances_;
    ivec atom_type_list_;
    double lincr_ = 0.0;
    double start_dist_ = 1E-7;
    std::vector<std::tuple<std::string, _time_point>> timing_points_;
	bool non_spherical_densities_calculated_ = false;

    // Internal helper methods
    void setupPrototypeGrids(const WFN& wave, const ivec& atom_types);

    void generateIntegrationGrids(const WFN& wave, const cell& unit_cell,
        const ivec& atom_list);
    void getIntegrationGrid1D(const WFN& wave, const int atom_1, const int atom_2, const int num_points, const double padding);

    void calculateSphericalDensities(const WFN& wave, const cell& unit_cell, const ivec& atom_list, vec2& single_spherical_density, vec2& combined_spherical_density, const std::vector<std::pair<vec, vec>> sig_pop = {});
    void calculateHirshfeldWeights(const WFN& wave, const cell& unit_cell, const ivec& atom_list);
    std::vector<std::pair<vec, vec>> calculateMBISWeights(const WFN& wave, const cell& unit_cell, const ivec& atom_list);
    void calculateEMBISWeights(const WFN& wave, const cell& unit_cell, const ivec& atom_list, const std::vector<std::pair<vec, vec>>& MBIS_weights);
    void pruneGrid();

    void addTimingPoint(const std::string& label) {
        timing_points_.emplace_back(label, std::chrono::high_resolution_clock::now());
    };
public:
    GridManager(const GridConfiguration& config = GridConfiguration{});
    ~GridManager() = default;

    // Main interface methods
    void setup3DGridsForMolecule(const WFN& wave, const bvec& needs_grid,
        const ivec& atom_list, const cell& unit_cell = cell());

    void setup1DGridsForMolecule(const WFN& wave, const int atom_1, const int atom_2, const int gridpoints, const double padding);

    void calculateNonSphericalDensities(const WFN& wave, const cell& unit_cell);

    PartitionResults calculatePartitionedCharges(const WFN& wave, const cell& unit_cell = cell());

    void getDensityVectors(const WFN& wave, const ivec& atom_list, vec2& d1, vec2& d2, vec2& d3, vec2& dens);

    // Configuration and data access
    void setConfiguration(const GridConfiguration& config) { config_ = config; }
    const GridConfiguration& getConfiguration() const { return config_; }
    const GridData& getGridData() const { return grid_data_; }
    int getTotalGridPoints() const { return grid_data_.total_points; }
	int getNumPointsForAtom(const int& atom_index) const { return grid_data_.num_points_per_atom[atom_index]; }

    // Utility methods
    static ivec identifyAtomTypes(const WFN& wave, const bvec& needs_grid);
    static bvec determineAtomsNeedingGrids(const WFN& wave, const ivec& asym_atom_list);

    void addTimingInfoToVecs(std::vector<_time_point>& time_points, svec& time_descriptions) {
        for (const auto [label, time] : timing_points_) {
            time_points.push_back(time);
            time_descriptions.push_back(label);
        }
    }
    void writeSimpleGrid(const std::filesystem::path& filename, const vec2& grid_points, std::vector<std::pair<std::string, vec>> data) const;

    vec evaluateFunctionOnGrid(const vec2& grid_points, std::function<double(double, double, double)> func) const;
};