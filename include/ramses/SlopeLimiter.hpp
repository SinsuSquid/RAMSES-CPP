#ifndef RAMSES_SLOPE_LIMITER_HPP
#define RAMSES_SLOPE_LIMITER_HPP

#include "Types.hpp"

namespace ramses {

/**
 * @brief Implements TVD slope limiters for MUSCL reconstruction.
 * 
 * Replicates logic from hydro/umuscl.f90.
 */
class SlopeLimiter {
public:
    /**
     * @brief Computes a limited slope given cell values.
     * @param ql Value in the left cell
     * @param qc Value in the central cell
     * @param qr Value in the right cell
     * @param slope_type 1=MinMod, 2=MonCen, 3=Average, etc.
     * @param theta Parameter for generalized MonCen (default 1.5 - 2.0)
     * @param nu Courant number for velocity-dependent limiters (default 0.0)
     */
    static real_t compute_slope(real_t ql, real_t qc, real_t qr, int slope_type, real_t theta = 2.0, real_t nu = 0.0);
};

} // namespace ramses

#endif // RAMSES_SLOPE_LIMITER_HPP
