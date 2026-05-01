#ifndef RAMSES_HYDRO_SOLVER_HPP
#define RAMSES_HYDRO_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <vector>
#include <memory>

namespace ramses {

/**
 * @brief Implements the Godunov hydro solver logic.
 */
class HydroSolver {
public:
    HydroSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

    void godunov_fine(int ilevel, real_t dt, real_t dx);
    void set_unew(int ilevel);
    void set_uold(int ilevel);
    
    void ctoprim(const real_t u[], real_t q[], real_t gamma);
    void interpol_hydro(const real_t u1[7][20], real_t u2[8][20]);

    void godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx);

    /**
     * @brief Computes the dynamic timestep based on CFL condition.
     */
    real_t compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor);

    /**
     * @brief Adds gravity source terms to unew.
     */
    void add_gravity_source_terms(int ilevel, real_t dt);

    void get_diagnostics(int ilevel, real_t dx, real_t& mind, real_t& maxv, real_t& min_t, real_t& max_t);

    struct LocalStencil {
        real_t uloc[6][6][6][20]; 
        bool refined[6][6][6];
    };

private:
    AmrGrid& grid_;
    Config& config_;

    void trace(const real_t q[], const real_t dq[], real_t dtdx, real_t qm[], real_t qp[], real_t gamma);
    void compute_slopes(int idc, const int icelln[6], int idim, real_t dq[5], int slope_type);

    std::unique_ptr<LocalStencil> stencil_ptr_ = std::make_unique<LocalStencil>();
};

} // namespace ramses

#endif // RAMSES_HYDRO_SOLVER_HPP
