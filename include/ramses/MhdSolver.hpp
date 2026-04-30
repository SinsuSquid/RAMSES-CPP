#ifndef RAMSES_MHD_SOLVER_HPP
#define RAMSES_MHD_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <vector>
#include <memory>

namespace ramses {

/**
 * @brief Implements the Magnetohydrodynamics (MHD) solver with Constrained Transport (CT).
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

    struct LocalStencil {
        real_t uloc[6][6][6][20]; 
        bool refined[6][6][6];
    };

    /**
     * @brief Computes the maximum divergence of B on a specific level.
     */
    real_t compute_max_div_b(int ilevel, real_t dx);

    void get_diagnostics(int ilevel, real_t dx, real_t& mind, real_t& maxv, real_t& maxdivb);

private:
    AmrGrid& grid_;
    Config& config_;

    void godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx);
    void gather_stencil(int igrid, int ilevel, LocalStencil& stencil);

    // MUSCL-Hancock trace
    void trace(const real_t qloc[6][6][6][20], const real_t bfloc[6][6][6][3][2],
               const real_t dq[6][6][6][3][20], const real_t dbf[6][6][6][3][2],
               real_t dt, real_t dx,
               real_t qm[6][6][6][3][20], real_t qp[6][6][6][3][20]);

    // Riemann solver wrappers
    void cmpflxm(const real_t qm[6][6][6][3][20], const real_t qp[6][6][6][3][20],
                 const real_t bfloc[6][6][6][3][2],
                 int idim, real_t gamma, real_t flux[6][6][6][20]);

    // MHD-specific methods
    void hlld(const real_t* qleft, const real_t* qright, real_t* fgdnv, real_t gamma);
    void llf(const real_t* qleft, const real_t* qright, real_t* fgdnv, real_t gamma);
    void find_mhd_flux(const real_t* qvar, real_t* cvar, real_t* ff, real_t gamma);
    void find_speed_fast(const real_t* qvar, real_t& vel_info, real_t gamma);

    void ctoprim(const real_t u[20], real_t q[20], real_t bf[3][2], real_t gamma);
    
    // Heap-allocated workspace to avoid stack overflow
    std::unique_ptr<LocalStencil> stencil_ptr_ = std::make_unique<LocalStencil>();
};

} // namespace ramses

#endif
