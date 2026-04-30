#ifndef RAMSES_RT_SOLVER_HPP
#define RAMSES_RT_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <vector>
#include <string>

namespace ramses {

/**
 * @brief Implements Radiative Transfer (RT) using the M1 closure scheme.
 * 
 * Replicates logic from rt/rt_godunov_fine.f90 and rt_flux_module.f90.
 */
class RtSolver {
public:
    RtSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

    /**
     * @brief Initializes RT parameters and loads eigenvalue tables.
     */
    void initialize();

    /**
     * @brief Performs the RT update step on a specific AMR level.
     */
    void godunov_fine(int ilevel, real_t dt, real_t dx);

    /**
     * @brief Synchronizes unew and uold for RT variables.
     */
    void set_uold(int ilevel);
    void set_unew(int ilevel);

private:
    void load_hll_eigenvalues();
    void rt_godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx);

    AmrGrid& grid_;
    Config& config_;

    // HLL Eigenvalue tables [101][101]
    std::vector<std::vector<real_t>> lambda1;
    std::vector<std::vector<real_t>> lambda4;

    int nGroups = 0;
    real_t smallNp = 1e-10;
};

} // namespace ramses

#endif // RAMSES_RT_SOLVER_HPP
