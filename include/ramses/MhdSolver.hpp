#ifndef RAMSES_MHD_SOLVER_HPP
#define RAMSES_MHD_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"

namespace ramses {

/**
 * @brief Godunov-type solver for Magnetohydrodynamics (MHD).
 */
class MhdSolver {
public:
    MhdSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

    void godunov_fine(int ilevel, real_t dt, real_t dx);
    void set_uold(int ilevel);
    real_t compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor);

    void set_nener(int nener) { nener_ = nener; }
    void interpol_mhd(const real_t u1[7][20], real_t u2[8][20]) {
        for(int i=0; i<8; ++i) for(int iv=0; iv<grid_.nvar; ++iv) u2[i][iv] = u1[0][iv];
    }

private:
    AmrGrid& grid_;
    Config& config_;
    int nener_ = 0;
};

} // namespace ramses

#endif // RAMSES_MHD_SOLVER_HPP
