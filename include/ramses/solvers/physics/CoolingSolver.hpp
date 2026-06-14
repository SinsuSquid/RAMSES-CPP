#ifndef RAMSES_COOLING_SOLVER_HPP
#define RAMSES_COOLING_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Handles gas cooling and heating processes.
 * 
 * Implements models from hydro/cooling_module.f90 and cooling_module_ism.f90.
 */
class CoolingSolver {
public:
    CoolingSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}
    virtual ~CoolingSolver();

    /**
     * @brief Applies cooling to a specific AMR level.
     */
    virtual void apply_cooling(int ilevel, real_t dt);

protected:
    /**
     * @brief The analytic ISM cooling model (Hennebelle 2005).
     */
    void solve_cooling_ism(real_t& d, real_t& e_tot, real_t dt_sec, real_t dx, real_t mu);

    /**
     * @brief Low temperature cooling rates (T < 1e4 K).
     */
    real_t cooling_low(real_t T, real_t n);

    /**
     * @brief High temperature cooling rates (T > 1e4 K).
     */
    real_t cooling_high(real_t T, real_t n);

    AmrGrid& grid_;
    Config& config_;
};

} // namespace ramses

#endif // RAMSES_COOLING_SOLVER_HPP
