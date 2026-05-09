#ifndef RAMSES_RELATIVISTIC_RIEMANN_SOLVER_HPP
#define RAMSES_RELATIVISTIC_RIEMANN_SOLVER_HPP

#include "Types.hpp"
#include <vector>
#include <string>

namespace ramses {

/**
 * @brief Standalone Riemann solver for hydrodynamics.
 * 
 * Supports LLF and eventually others (HLL, HLLC, Exact).
 */
class RelativisticRiemannSolver {
public:
    static void solve_llf(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, const std::string& eos);
    static void solve_hll(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, const std::string& eos);
    static void solve_hllc(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, const std::string& eos);

    /**
     * @brief Helper to compute conservative variables from primitive ones.
     */
    static void find_rhd_flux(const real_t q[], real_t u[], real_t f[], real_t gamma, const std::string& eos);

    /**
     * @brief Helper to compute fast wave speed.
     */
    static void find_speed_fast(const real_t q[], real_t& cs, real_t gamma, const std::string& eos);

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

#endif // RAMSES_RELATIVISTIC_RIEMANN_SOLVER_HPP
