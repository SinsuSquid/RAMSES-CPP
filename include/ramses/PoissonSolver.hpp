#ifndef RAMSES_POISSON_SOLVER_HPP
#define RAMSES_POISSON_SOLVER_HPP

#include "AmrGrid.hpp"

namespace ramses {

/**
 * @brief Multigrid-based Poisson solver for gravity.
 */
class PoissonSolver {
public:
    PoissonSolver(AmrGrid& grid) : grid_(grid) {}

    void solve(int ilevel);

private:
    void smooth(int ilevel);
    void restrict(int ilevel);
    void prolong(int ilevel);
    void vcycle(int ilevel);

    AmrGrid& grid_;
    std::vector<real_t> res; // Residual workspace [ncell]
};

} // namespace ramses

#endif // RAMSES_POISSON_SOLVER_HPP
