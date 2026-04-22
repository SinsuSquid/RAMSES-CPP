#ifndef RAMSES_HYDRO_SOLVER_HPP
#define RAMSES_HYDRO_SOLVER_HPP

#include "AmrGrid.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Implements the Godunov hydro solver logic.
 * 
 * Replicates hydro/godunov_fine.f90 and related modules.
 */
class HydroSolver {
public:
    HydroSolver(AmrGrid& grid) : grid_(grid) {}

    /**
     * @brief Performs one hydro step for a given level.
     */
    void godunov_fine(int ilevel);

    /**
     * @brief Sets unew to uold before the hydro step.
     */
    void set_unew(int ilevel);

    /**
     * @brief Sets uold to unew after the hydro step.
     */
    void set_uold(int ilevel);

private:
    AmrGrid& grid_;

    /**
     * @brief Inner worker for Godunov solver on a batch of grids.
     */
    void godfine1(const std::vector<int>& ind_grid, int ilevel);
    static void ctoprim(const real_t u[], real_t q[], real_t gamma);
};

} // namespace ramses

#endif // RAMSES_HYDRO_SOLVER_HPP
