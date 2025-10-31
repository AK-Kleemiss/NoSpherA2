#pragma once

#include "convenience.h"
#include "wfn_class.h"
#include "AtomGrid.h"
#include "cell.h"

double make_sphericals(
    vec2& dens,
    vec2& dist,
    const ivec& atom_type_list,
    std::ostream& file,
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
    enum GridIndex { X = 0, Y = 1, Z = 2, WEIGHT = 3, HIRSH_WEIGHT = 4, BECKE_WEIGHT = 5, TFVC_WEIGHT = 6, WFN_DENSITY = 7};
    vec3 atomic_grids;           // [atom][coord_type][point]
    ivec num_points_per_atom;    // Number of points for each atom
    int total_points = 0;
    
    void clear();
    void resizeForAtoms(int num_atoms);
};

struct PartitionResults {
    vec2 atom_charges;           // [atom][partition type]
    vec overall_charges;      // Total charges by partition type
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
    std::vector<std::tuple<std::string, _time_point>> timing_points_;
    
    // Internal helper methods
    void setupPrototypeGrids(const WFN& wave, const ivec& atom_types);
    void generateIntegrationGrids(const WFN& wave, const cell& unit_cell, 
                                  const bvec& needs_grid, const ivec& atom_list);
    void calculateNonSphericalDensities(const WFN& wave, const cell& unit_cell);
    void calculateHirshfeldWeights(const WFN& wave, const cell& unit_cell, const ivec& atom_list);
    void pruneGrid();
    
    void addTimingPoint(const std::string& label) {
        timing_points_.emplace_back(label, std::chrono::high_resolution_clock::now());
    };
public:
    GridManager(const GridConfiguration& config = GridConfiguration{});
    ~GridManager() = default;
    
    // Main interface methods
    void setupGridsForMolecule(const WFN& wave, const bvec& needs_grid, 
                               const ivec& atom_list, const cell& unit_cell = cell());
    
    PartitionResults calculatePartitionedCharges(const WFN& wave, const cell& unit_cell = cell());
    
    void getDensityVectors(const WFN& wave, vec2& d1, vec2& d2, vec2& d3, vec2& dens);
    
    // Configuration and data access
    void setConfiguration(const GridConfiguration& config) { config_ = config; }
    const GridConfiguration& getConfiguration() const { return config_; }
    const GridData& getGridData() const { return grid_data_; }
    int getTotalGridPoints() const { return grid_data_.total_points; }
    
    // Utility methods
    static ivec identifyAtomTypes(const WFN& wave, const bvec& needs_grid);
    static bvec determineAtomsNeedingGrids(const WFN& wave, const ivec& asym_atom_list);
    
    void add_timing_info_to_vecs(std::vector<_time_point>& time_points, svec& time_descriptions) {
        for (const auto [label, time] : timing_points_) {
            time_points.push_back(time);
            time_descriptions.push_back(label);
        }
    }
};