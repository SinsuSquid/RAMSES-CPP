#ifndef RAMSES_TURBULENCE_SOLVER_HPP
#define RAMSES_TURBULENCE_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <vector>
#include <string>

namespace ramses {

/**
 * @brief Module for turbulent forcing (stochastic driving).
 * 
 * Supports Ornstein-Uhlenbeck process in Fourier space and 
 * spatial interpolation to AMR cells.
 */
class TurbulenceSolver {
public:
    TurbulenceSolver(AmrGrid& grid, Config& config);
    ~TurbulenceSolver();

    /**
     * @brief Initialize the turbulent fields.
     */
    void init();

    /**
     * @brief Update the turbulent fields (OU process).
     * @param dt Current timestep.
     */
    void update_fields(real_t dt);

    /**
     * @brief Calculate and apply turbulent forcing to a level.
     * @param ilevel Level to apply forcing.
     * @param dt Timestep for this level.
     */
    void apply_forcing(int ilevel, real_t dt);

private:
    AmrGrid& grid_;
    Config& config_;

    int turb_gs_ = 64; // Grid size for uniform forcing field
    real_t turb_min_rho_ = 0.0;
    real_t sol_frac_ = 1.0; // Solenoidal fraction

    // Uniform forcing field (spatial)
    std::vector<real_t> afield_now_; // Size: NDIM * turb_gs_^NDIM

    // Internal methods
    void interpolate_force(const real_t x[3], real_t force[3]);
    void generate_forcing_field(); // If FFTW available, or fallback
};

} // namespace ramses

#endif // RAMSES_TURBULENCE_SOLVER_HPP
