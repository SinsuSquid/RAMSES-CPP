#include "ramses/solvers/hydro/SlopeLimiter.hpp"
#include <cmath>
#include <algorithm>

namespace ramses {

real_t SlopeLimiter::compute_slope(real_t ql, real_t qc, real_t qr, int slope_type, real_t theta, real_t nu) {
    if (slope_type == 0) return 0.0;

    real_t dlft = qc - ql;
    real_t drgt = qr - qc;

    // Superbee (4) and Ultrabee (5) limiters
    if (slope_type == 4 || slope_type == 5) {
        real_t dlft_s, drgt_s;
        if (slope_type == 4) {
            dlft_s = 2.0 / (1.0 + nu) * dlft;
            drgt_s = 2.0 / (1.0 - nu) * drgt;
        } else {
            if (nu >= 0.0) {
                dlft_s = 2.0 / (nu + 1e-10) * dlft;
                drgt_s = 2.0 / (1.0 - nu) * drgt;
            } else {
                dlft_s = 2.0 / (1.0 + nu) * dlft;
                drgt_s = 2.0 / (-nu + 1e-10) * drgt;
            }
        }
        real_t dcen = 0.5 * (dlft + drgt);
        real_t dsgn = (dcen >= 0.0) ? 1.0 : -1.0;
        real_t dlim = std::min(std::abs(dlft_s), std::abs(drgt_s));
        if (dlft * drgt <= 0.0) dlim = 0.0;
        return dsgn * dlim;
    }

    // Centered slope (6)
    if (slope_type == 6) {
        return 0.5 * (dlft + drgt);
    }

    // Standard TVD limiters require same sign for left/right differences
    if (dlft * drgt <= 0.0) return 0.0;

    real_t dcen = 0.5 * (dlft + drgt);
    real_t dsgn = (dcen >= 0.0) ? 1.0 : -1.0;
    real_t dlim = 0.0;

    if (slope_type == 1) { // MinMod
        dlim = std::min(std::abs(dlft), std::abs(drgt));
    } else if (slope_type == 2 || slope_type == 3) { // MonCen
        dlim = std::min({2.0 * std::abs(dlft), 2.0 * std::abs(drgt), std::abs(dcen)});
    } else if (slope_type == 7) { // van Leer
        dlim = 2.0 * std::abs(dlft) * std::abs(drgt) / (std::abs(dlft) + std::abs(drgt));
    } else if (slope_type == 8) { // Generalized MonCen
        dlim = std::min({theta * std::abs(dlft), theta * std::abs(drgt), std::abs(dcen)});
    } else { // Fallback to MinMod
        dlim = std::min(std::abs(dlft), std::abs(drgt));
    }

    return dsgn * dlim;
}

} // namespace ramses
