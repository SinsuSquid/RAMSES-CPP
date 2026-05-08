#ifndef RAMSES_POISSON_SOLVER_HPP
#define RAMSES_POISSON_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Multigrid Poisson solver for self-gravity.
 */
class PoissonSolver {
public:
    PoissonSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}
    virtual ~PoissonSolver();

    virtual void solve(int ilevel, real_t aexp = 1.0, real_t omega_m = 0.3, real_t rho_tot = 0.0);
    virtual void compute_force(int ilevel);

protected:
    virtual void vcycle(int ilevel, real_t fourpi, real_t rho_tot);
    virtual void smooth(int ilevel, real_t fourpi, real_t rho_tot);
    virtual void restrict(int ilevel);
    virtual void prolong(int ilevel);

    AmrGrid& grid_;
    Config& config_;
    std::vector<real_t> res; 
};

} // namespace ramses

#endif // RAMSES_POISSON_SOLVER_HPP
