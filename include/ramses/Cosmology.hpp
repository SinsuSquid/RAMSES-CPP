#ifndef RAMSES_COSMOLOGY_HPP
#define RAMSES_COSMOLOGY_HPP

#include "Types.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Handles cosmological expansion and Friedmann equations.
 * 
 * Replicates amr/init_time.f90 (Friedman model).
 */
class Cosmology {
public:
    Cosmology() = default;

    /**
     * @brief Computes look-up tables for expansion factors and Hubble parameters.
     */
    void solve_friedman(real_t omega_m, real_t omega_l, real_t omega_k, real_t aexp_ini);

    /**
     * @brief Interpolates cosmological parameters at a given conformal time tau.
     */
    void get_cosmo_params(real_t tau, real_t& aexp, real_t& hexp, real_t& texp);

    /**
     * @brief Interpolates conformal time tau at a given expansion factor aexp.
     */
    real_t get_tau(real_t aexp);

    /**
     * @brief Computes the growth factor f = d log D1 / d log a.
     */
    real_t fpeebl(real_t a);

private:
    real_t d1a(real_t a);
    real_t rombint(real_t a, real_t b, real_t tol);
    real_t fy(real_t a, real_t omega_m, real_t omega_l, real_t omega_k);
    real_t dadtau(real_t a, real_t omega_m, real_t omega_l, real_t omega_k);
    real_t dadt(real_t a, real_t omega_m, real_t omega_l, real_t omega_k);

    std::vector<real_t> aexp_frw;
    std::vector<real_t> hexp_frw;
    std::vector<real_t> tau_frw;
    std::vector<real_t> t_frw;

    real_t omega_m_ = 0.3;
    real_t omega_l_ = 0.7;
    real_t omega_k_ = 0.0;
};

} // namespace ramses

#endif // RAMSES_COSMOLOGY_HPP
