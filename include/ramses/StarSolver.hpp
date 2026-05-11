#ifndef RAMSES_STAR_SOLVER_HPP
#define RAMSES_STAR_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <random>

namespace ramses {

/**
 * @brief Solver for star formation processes.
 * 
 * Handles gas-to-star conversion based on Poisson statistics and local gas properties.
 */
class StarSolver {
public:
    StarSolver(AmrGrid& grid, Config& config);
    ~StarSolver();

    void init();
    void form_stars(int ilevel, real_t dt);

private:
    AmrGrid& grid_;
    Config& config_;

    // Parameters
    real_t n_star_;      ///< SF density threshold [H/cc]
    real_t eps_star_;    ///< SF efficiency per free-fall time
    real_t m_star_;      ///< Base star particle mass [Solar Masses]
    real_t T2_star_;     ///< Polytropic temperature threshold
    real_t g_star_;      ///< Polytropic index
    
    // Internal state
    int nstar_tot_ = 0;
};

} // namespace ramses

#endif // RAMSES_STAR_SOLVER_HPP
