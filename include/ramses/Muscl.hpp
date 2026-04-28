#ifndef RAMSES_MUSCL_HPP
#define RAMSES_MUSCL_HPP

#include "Types.hpp"
#include <vector>

namespace ramses {

/**
 * @brief MUSCL slope limiting and reconstruction.
 */
class Muscl {
public:
    // Compute slopes using a specific limiter
    static void compute_slopes(const real_t* q, real_t* dq, int nx, int nvar, int slope_type);

    // Reconstruct left/right interface states
    static void reconstruct(const real_t* q, const real_t* dq, real_t* qL, real_t* qR, int nx, int nvar);

    // MUSCL prediction
    static void predict(const real_t q[], const real_t dq[], const real_t s0[], 
                       real_t dt_dx, real_t qm[], real_t qp[], int nvar);
};

} // namespace ramses

#endif
