#ifndef RAMSES_RIEMANN_SOLVER_HPP
#define RAMSES_RIEMANN_SOLVER_HPP

#include "Types.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Standalone Riemann solver for hydrodynamics.
 * 
 * Supports LLF and eventually others (HLL, HLLC, Exact).
 */
class RiemannSolver {
public:
    static void solve_llf(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma);
    static void solve_hllc(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma);

private:
    /**
     * @brief Helper to compute conservative variables from primitive ones.
     */
    static void prim_to_cons(const real_t q[], real_t u[], real_t gamma);

    /**
     * @brief Helper to compute flux from primitive state directly.
     */
    static void compute_flux(const real_t q[], real_t f[], real_t gamma);
};

} // namespace ramses

#endif // RAMSES_RIEMANN_SOLVER_HPP
