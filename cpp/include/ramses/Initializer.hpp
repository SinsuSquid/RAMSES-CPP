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

    AmrGrid& grid_;
    const Config& config_;
};

} // namespace ramses

#endif // RAMSES_INITIALIZER_HPP
