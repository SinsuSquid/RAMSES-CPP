#ifndef RAMSES_POISSON_SOLVER_HPP
#define RAMSES_POISSON_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"

namespace ramses {

/**
 * @brief Multigrid-based Poisson solver for gravity.
 */
class PoissonSolver {
public:
    PoissonSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

    void solve(int ilevel);
    void compute_force(int ilevel);

private:
    void smooth(int ilevel);
    void restrict(int ilevel);
    void prolong(int ilevel);
    void vcycle(int ilevel);

    AmrGrid& grid_;
    Config& config_;
    std::vector<real_t> res; // Residual workspace [ncell]
};

} // namespace ramses

#endif // RAMSES_POISSON_SOLVER_HPP
