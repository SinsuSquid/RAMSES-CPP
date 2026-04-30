#ifndef RAMSES_HYDRO_SOLVER_HPP
#define RAMSES_HYDRO_SOLVER_HPP

#include "AmrGrid.hpp"
#include <vector>
#include <memory>

namespace ramses {

/**
 * @brief Implements the Godunov hydro solver logic.
 */
class HydroSolver {
public:
    HydroSolver(AmrGrid& grid) : grid_(grid) {}

    void godunov_fine(int ilevel, real_t dt, real_t dx);
    void set_unew(int ilevel);
    void set_uold(int ilevel);
    
    void ctoprim(const real_t u[], real_t q[], real_t gamma);
    void interpol_hydro(const real_t u1[7][20], real_t u2[8][20]);

    struct LocalStencil {
        real_t uloc[6][6][6][20]; 
        bool refined[6][6][6];
    };

    void godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx);

    /**
     * @brief Computes the dynamic timestep based on CFL condition.
     */
    real_t compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor);

    /**
     * @brief Adds gravity source terms to unew.
     */
    void add_gravity_source_terms(int ilevel, real_t dt);

    void get_diagnostics(int ilevel, real_t dx, real_t& mind, real_t& maxv);

private:
    AmrGrid& grid_;

    void gather_stencil(int igrid, int ilevel, LocalStencil& stencil);
    
    // Heap-allocated workspace to avoid stack overflow
    std::unique_ptr<LocalStencil> stencil_ptr_ = std::make_unique<LocalStencil>();
};

} // namespace ramses

#endif // RAMSES_HYDRO_SOLVER_HPP
