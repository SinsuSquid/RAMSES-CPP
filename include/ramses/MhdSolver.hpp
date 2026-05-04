#ifndef RAMSES_MHD_SOLVER_HPP
#define RAMSES_MHD_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <memory>

namespace ramses {

/**
 * @brief Stencil for 3D MHD calculations.
 */
struct LocalStencil {
    real_t uloc[6][6][6][20];
    real_t bfloc[6][6][6][3][2];
    bool refined[6][6][6];
};

/**
 * @brief Godunov-type solver for Magnetohydrodynamics (MHD).
 */
class MhdSolver {
public:
    MhdSolver(AmrGrid& grid, Config& config);
    ~MhdSolver();

    void godunov_fine(int ilevel, real_t dt, real_t dx);
    void set_uold(int ilevel);
    void set_unew(int ilevel);
    real_t compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor);
    real_t compute_max_div_b(int ilevel, real_t dx);
    void get_diagnostics(int il, real_t dx, real_t& mi, real_t& mv, real_t& mb);

    void gather_stencil(int igrid, int ilevel, LocalStencil& stencil);
    void ctoprim(const real_t u[20], real_t q[20], real_t bf[3][2], real_t gamma);
    void trace(const real_t qloc[6][6][6][20], const real_t bfloc[6][6][6][3][2], const real_t dq[6][6][6][3][20], const real_t dbf[6][6][6][3][2], real_t dt, real_t dx, real_t qm[6][6][6][3][20], real_t qp[6][6][6][3][20]);
    void cmpflxm(const real_t qm[6][6][6][3][20], const real_t qp[6][6][6][3][20], const real_t bfloc[6][6][6][3][2], int idim, real_t gamma, real_t flux[6][6][6][20]);
    void godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx);

    void llf(const real_t* ql, const real_t* qr, real_t* f, real_t gamma);
    void hlld(const real_t* ql_in, const real_t* qr_in, real_t* fgdnv, real_t gamma);
    void find_mhd_flux(const real_t* q, real_t* c, real_t* f, real_t gamma);
    void find_speed_fast(const real_t* q, real_t& v, real_t gamma);

    void set_nener(int nener) { nener_ = nener; }
    void interpol_mhd(const real_t u1[7][20], real_t u2[8][20]) {
        for(int i=0; i<8; ++i) for(int iv=0; iv<grid_.nvar; ++iv) u2[i][iv] = u1[0][iv];
    }

private:
    AmrGrid& grid_;
    Config& config_;
    int nener_ = 0;
    std::unique_ptr<LocalStencil> stencil_ptr_;
};

} // namespace ramses

#endif // RAMSES_MHD_SOLVER_HPP
