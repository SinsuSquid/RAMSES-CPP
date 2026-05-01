#ifndef RAMSES_INITIALIZER_HPP
#define RAMSES_INITIALIZER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <string>
#include <vector>

namespace ramses {

/**
 * @brief Handles initial condition generation.
 */
class Initializer {
public:
    Initializer(AmrGrid& grid, const Config& config) : grid_(grid), config_(config) {}

    /**
     * @brief Applies initial conditions to the grid based on namelist.
     */
    void apply_all();

private:
    /**
     * @brief Port of region_condinit from condinit.f90
     */
    void region_condinit();

    void initialize_cell(int i, int ilevel, int nreg, real_t gam, int nener,
                         const std::vector<real_t>& drs, const std::vector<real_t>& prs,
                         const std::vector<real_t>& urs, const std::vector<real_t>& vrs,
                         const std::vector<real_t>& wrs, const std::vector<real_t>& Ars,
                         const std::vector<real_t>& Brs, const std::vector<real_t>& Crs,
                         const std::vector<real_t>& xcs, const std::vector<real_t>& lxs,
                         const std::vector<real_t>& ycs, const std::vector<real_t>& lys,
                         const std::vector<real_t>& zcs, const std::vector<real_t>& lzs,
                         const std::vector<std::vector<real_t>>& prads,
                         real_t db, real_t pb, real_t ub, real_t vb, real_t wb,
                         real_t Ab, real_t Bb, real_t Cb, int total_matches[]);

    /**
     * @brief Analytic disk potential initial conditions.
     */
    void ana_disk_potential_condinit();

    /**
     * @brief Orszag-Tang initial conditions.
     */
    void orzag_tang_condinit();

    AmrGrid& grid_;
    const Config& config_;
};

} // namespace ramses

#endif // RAMSES_INITIALIZER_HPP
