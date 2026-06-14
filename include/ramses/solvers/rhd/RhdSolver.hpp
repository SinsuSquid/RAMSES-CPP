#ifndef RAMSES_RHD_SOLVER_HPP
#define RAMSES_RHD_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <string>
#include <vector>

namespace ramses {

/**
 * @brief Godunov-type solver for Relativistic Hydrodynamics.
 */
class RhdSolver {
public:
    RhdSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}
    virtual ~RhdSolver();

    virtual void godunov_fine(int ilevel, real_t dt, real_t dx);
    virtual void set_uold(int ilevel);
    virtual void set_unew(int ilevel);
    virtual real_t compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor);
    virtual void get_diagnostics(int ilevel, real_t dx, real_t& min_d, real_t& max_v, real_t& min_t, real_t& max_t);
    virtual void add_gravity_source_terms(int ilevel, real_t dt);
    virtual void synchro_hydro_fine(int ilevel, real_t dt);
    
    void set_nener(int nener) { nener_ = nener; }
    void set_nvar_hydro(int nvar) { nvar_hydro_ = nvar; }
    virtual void interpol_hydro(const real_t u1[7][64], real_t u2[8][64]);

protected:
    AmrGrid& grid_;
    Config& config_;
    int nener_ = 0;
    int nvar_hydro_ = 5;
    std::vector<real_t> qm_level_, qp_level_; 

private:
    void ctoprim(const real_t u[], real_t q[], real_t gamma);
    void compute_slopes(int idc, const int icelln[6], int idim, real_t dq[20], int slope_type);
    void trace(const real_t q[], const real_t dq[], real_t dtdx, real_t qm[], real_t qp[], real_t gamma);
    
    // RHD specific
    void newton_raphson_mignone(real_t D, real_t M, real_t E, real_t gamma, real_t& R);
};

} // namespace ramses

#endif // RAMSES_RHD_SOLVER_HPP
