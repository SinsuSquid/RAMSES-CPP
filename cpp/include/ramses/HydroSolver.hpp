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

    /**
     * @brief Performs hydro interpolation (prolongation) from father cells to children.
     * @param u1 Input 7 conservative states (center + 6 neighbors)
     * @param u2 Output 8 conservative states (children)
     */
    static void interpol_hydro(const real_t u1[7][5], real_t u2[8][5]);

    /**
     * @brief A local 6x6x6 (3D) or equivalent 1D/2D stencil for Godunov solver.
     */
    struct LocalStencil {
        // [i][j][k][ivar] - using fixed sizes for performance
        // IU1=-1, IU2=4 => 6 cells in each direction
        real_t uloc[6][6][6][5]; 
        bool refined[6][6][6];
    };

    private:
    AmrGrid& grid_;

    /**
     * @brief Gathers a 6x6x6 stencil for a specific oct.
     */
    void gather_stencil(int igrid, int ilevel, LocalStencil& stencil);


#endif // RAMSES_HYDRO_SOLVER_HPP
