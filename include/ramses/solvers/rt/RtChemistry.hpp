#ifndef RAMSES_RT_CHEMISTRY_HPP
#define RAMSES_RT_CHEMISTRY_HPP

#include "Types.hpp"
#include <vector>
#include <cmath>

namespace ramses {

/**
 * @brief Handles RT-gas chemistry and ionization source terms.
 * 
 * Replicates logic from rt/rt_cooling_module.f90.
 */
class RtChemistry {
public:
    RtChemistry(int nGroups) : nGroups_(nGroups) {
        sig_n.assign(nGroups, std::vector<real_t>(3, 0.0)); // HI, HeI, HeII
        sig_e.assign(nGroups, std::vector<real_t>(3, 0.0));
        group_egy.assign(nGroups, 0.0);
    }

    /**
     * @brief Evolves ionization states and gas energy in a single cell.
     */
    void solve_chemistry(real_t& T2, real_t xion[3], real_t Np[], real_t Fp[][3], 
                         real_t nH, real_t dt, real_t aexp);

private:
    // Recombination and Ionization rates (Hydrogen)
    real_t get_alphaA_HII(real_t T) const;
    real_t get_alphaB_HII(real_t T) const;
    real_t get_beta_HI(real_t T) const;

    int nGroups_;
    std::vector<std::vector<real_t>> sig_n; // Cross-sections [cm2]
    std::vector<std::vector<real_t>> sig_e; // Energy-weighted cross-sections [cm2]
    std::vector<real_t> group_egy;         // Group mean energy [eV]
    
    const real_t X = 0.76; // Hydrogen mass fraction
    const real_t mH = 1.6733e-24;
    const real_t kB = 1.3806e-16;
};

} // namespace ramses

#endif
