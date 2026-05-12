#include "ramses/RtChemistry.hpp"
#include "ramses/Parameters.hpp"
#include <algorithm>

namespace ramses {

void RtChemistry::solve_chemistry(real_t& T2, real_t xion[3], real_t Np[], real_t Fp[][3], 
                                real_t nH, real_t dt, real_t aexp) {
    // T2 is T/mu [K]
    // xion[0] = xHII
    
    const real_t smallx = 1e-10;
    real_t xHII = std::max(smallx, std::min(1.0 - smallx, xion[0]));
    real_t T = T2 * (1.0 + xHII); // mu approx 1/(1+xHII) for H-only
    
    // Convert dt from code units to physical seconds for CGS rates
    real_t dt_sec = dt * params::units_time;

    // 1. Photo-ionization and heating rates
    real_t photo_ion_rate = 0.0;
    real_t photo_heat_rate = 0.0;
    for (int ig = 0; ig < nGroups_; ++ig) {
        photo_ion_rate += sig_n[ig][0] * Np[ig] * 3e10; // c approx 3e10
        photo_heat_rate += sig_e[ig][0] * Np[ig] * 3e10 * group_egy[ig] * 1.6022e-12; // erg
    }

    // 2. Cooling and Recombination rates
    real_t alphaB = get_alphaB_HII(T);
    real_t betaHI = get_beta_HI(T);

    // 3. Fully-implicit update for xHII
    // dx/dt = C_coeff * (1-x) - D_coeff * x
    real_t C_coeff = nH * betaHI + photo_ion_rate;
    real_t D_coeff = nH * alphaB;
    
    xHII = (xHII + C_coeff * dt_sec) / (1.0 + (C_coeff + D_coeff) * dt_sec);
    xion[0] = xHII;

    // 4. Energy update (simplified)
    // dE/dt = photo_heat_rate - cooling_rate
    // For now just update T2 assuming constant cooling or simple model
    T2 += (photo_heat_rate / (1.5 * nH * kB)) * dt_sec;
}

real_t RtChemistry::get_alphaB_HII(real_t T) const {
    real_t lambda = 315614.0 / std::max(T, 1.0);
    return 2.753e-14 * std::pow(lambda, 1.5) / std::pow(1.0 + std::pow(lambda / 2.74, 0.407), 2.242);
}

real_t RtChemistry::get_beta_HI(real_t T) const {
    real_t lambda = 315614.0 / std::max(T, 1.0);
    return 5.85e-11 * std::sqrt(T) * std::exp(-lambda / 2.0) / (1.0 + std::pow(lambda / 0.522, 0.47));
}

real_t RtChemistry::get_alphaA_HII(real_t T) const {
    real_t lambda = 315614.0 / std::max(T, 1.0);
    return 1.269e-13 * std::pow(lambda, 1.503) / std::pow(1.0 + std::pow(lambda / 0.522, 0.47), 1.923);
}

} // namespace ramses
