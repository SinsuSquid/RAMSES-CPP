#ifndef RAMSES_CLUMP_FINDER_HPP
#define RAMSES_CLUMP_FINDER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <vector>
#include <map>

namespace ramses {

struct Clump {
    int id;
    real_t peak_pos[3];
    real_t peak_dens;
    real_t total_mass;
    real_t com[3];
    real_t vel[3];
    int n_cells;
    std::vector<int> cell_indices;
};

/**
 * @brief Identifies clumps (density peaks and bound structures) in the simulation grid.
 */
class ClumpFinder {
public:
    ClumpFinder(AmrGrid& grid, Config& config);
    ~ClumpFinder();

    void find_clumps();
    const std::vector<Clump>& get_clumps() const { return clumps_; }

private:
    AmrGrid& grid_;
    Config& config_;
    std::vector<Clump> clumps_;

    real_t density_threshold_ = -1.0;
    real_t relevance_threshold_ = 2.0;
    real_t saddle_threshold_ = -1.0;

    void identify_peaks();
    void merge_irrelevant_peaks();
    void compute_clump_properties();
};

} // namespace ramses

#endif // RAMSES_CLUMP_FINDER_HPP
