#ifndef RAMSES_MHD_SOLVER_HPP
#define RAMSES_MHD_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Implements the Magnetohydrodynamics (MHD) solver.
 */
class MhdSolver {
public:
    MhdSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

    /**
     * @brief Performs one Godunov step for MHD on a specific level.
     */
    void godunov_fine(int ilevel, real_t dt, real_t dx);

    /**
     * @brief Computes the Courant-limited timestep for MHD.
     */
    real_t compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor);

    /**
     * @brief Sets unew to uold before the step.
     */
    void set_unew(int ilevel);

    /**
     * @brief Sets uold to unew after the step.
     */
    void set_uold(int ilevel);

private:
    AmrGrid& grid_;
    Config& config_;

    // MHD-specific methods
    void hlld(const real_t* qleft, const real_t* qright, real_t* fgdnv, real_t gamma);
    void find_mhd_flux(const real_t* qvar, real_t* cvar, real_t* ff, real_t gamma);
    void find_speed_fast(const real_t* qvar, real_t& vel_info, real_t gamma);
};

} // namespace ramses

#endif
