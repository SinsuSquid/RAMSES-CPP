#ifndef RAMSES_HYDRO_SOLVER_HPP
#define RAMSES_HYDRO_SOLVER_HPP

#include "AmrGrid.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Implements the Godunov hydro solver logic.
 */
class HydroSolver {
public:
    HydroSolver(AmrGrid& grid) : grid_(grid) {}

    void godunov_fine(int ilevel);
    void set_unew(int ilevel);
    void set_uold(int ilevel);

    static void ctoprim(const real_t u[], real_t q[], real_t gamma);
    static void interpol_hydro(const real_t u1[7][5], real_t u2[8][5]);

    struct LocalStencil {
        real_t uloc[6][6][6][5]; 
        bool refined[6][6][6];
    };

    void godfine1(const std::vector<int>& ind_grid, int ilevel);

    /**
     * @brief Computes the dynamic timestep based on CFL condition.
     */
    real_t compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor);

private:
    AmrGrid& grid_;

    void gather_stencil(int igrid, int ilevel, LocalStencil& stencil);
};

} // namespace ramses

#endif // RAMSES_HYDRO_SOLVER_HPP
