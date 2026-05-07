#ifndef RAMSES_RT_SOLVER_HPP
#define RAMSES_RT_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include "RtChemistry.hpp"
#include <vector>
#include <string>
#include <memory>

namespace ramses {

/**
 * @brief Implements Radiative Transfer (RT) using the M1 closure scheme.
 */
class RtSolver {
public:
    RtSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

    void initialize();
    void godunov_fine(int ilevel, real_t dt, real_t dx);
    void apply_source_terms(int ilevel, real_t dt);
    void set_uold(int ilevel);
    void set_unew(int ilevel);
    real_t get_c_speed() const { return rt_c_speed; }

    void set_nIons(int nIons) { nIons_ = nIons; }
    int get_nIons() const { return nIons_; }
    int get_nGroups() const { return nGroups; }

private:
    void load_hll_eigenvalues();
    void rt_godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx);
    
    // M1 closure specific methods
    void cmp_flux_tensors(int ncache, int iP0, int ilevel);
    void cmp_eigenvals(int ncache, int iP0, int ilevel);
    void inp_eigenvals(real_t ff, real_t omega, real_t& lmin, real_t& lmax);
    void cmp_rt_faces(int ncache, int iP0, int ilevel, real_t dt, real_t dx);

    AmrGrid& grid_;
    Config& config_;
    std::unique_ptr<RtChemistry> chem_;

    // HLL Eigenvalue tables [101][101]
    std::vector<std::vector<real_t>> lambda1;
    std::vector<std::vector<real_t>> lambda4;

    int nGroups = 0;
    int nIons_ = 3;
    real_t smallNp = 1e-10;
    real_t rt_c_speed = 1.0; // Reduced light speed in code units
    bool rt_use_hll = true;

    // Local stencil and flux buffers (equivalent to uloc, cFlx, flux, etc.)
    // Indexed as [ngrid][6][6][6][nvar] or similar
    std::vector<real_t> uloc;
    std::vector<real_t> cFlx; // [ngrid][6][6][6][ndim+1][ndim]
    std::vector<real_t> lmin, lmax; // [ngrid][6][6][6][ndim]
    std::vector<real_t> flux; // [ngrid][4][4][4][nrtvar][ndim]
    std::vector<bool> ok_stencil;
};

} // namespace ramses

#endif // RAMSES_RT_SOLVER_HPP
