#ifndef RAMSES_HYDRO_SOLVER_HPP
#define RAMSES_HYDRO_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <string>
#include <vector>

namespace ramses {

/**
 * @brief Godunov-type solver for Euler equations.
 */
class HydroSolver {
public:
    HydroSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

    void godunov_fine(int ilevel, real_t dt, real_t dx);
    void set_uold(int ilevel);
    void set_unew(int ilevel);
    real_t compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor);
    void get_diagnostics(int ilevel, real_t dx, real_t& min_d, real_t& max_v, real_t& min_t, real_t& max_t);
    void add_gravity_source_terms(int ilevel, real_t dt);
    
    void set_nener(int nener) { nener_ = nener; }
    void interpol_hydro(const real_t u1[7][20], real_t u2[8][20]);

private:
    void ctoprim(const real_t u[], real_t q[], real_t gamma);
    void compute_slopes(int idc, const int icelln[6], int idim, real_t dq[20], int slope_type);
    void trace(const real_t q[], const real_t dq[], real_t dtdx, real_t qm[], real_t qp[], real_t gamma);

    AmrGrid& grid_;
    Config& config_;
    int nener_ = 0;
    std::vector<real_t> qm_level_, qp_level_; 
};

} // namespace ramses

#endif // RAMSES_HYDRO_SOLVER_HPP
