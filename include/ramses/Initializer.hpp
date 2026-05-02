#ifndef RAMSES_INITIALIZER_HPP
#define RAMSES_INITIALIZER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <string>

namespace ramses {

/**
 * @brief Handles initial condition setup.
 * 
 * Ported from amr/init_amr.f90 and hydro/init_hydro.f90.
 */
class Initializer {
public:
    Initializer(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

    void apply_all();
    void region_condinit(int ilevel);

private:
    AmrGrid& grid_;
    Config& config_;
};

} // namespace ramses

#endif // RAMSES_INITIALIZER_HPP
