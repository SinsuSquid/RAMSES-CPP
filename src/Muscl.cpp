#include "ramses/Muscl.hpp"
#include <algorithm>
#include <cmath>

namespace ramses {

void Muscl::compute_slopes(const real_t* q, real_t* dq, int nx, int nvar, int slope_type) {
    // 1D MinMod slope limiter
    for (int i = 1; i < nx - 1; ++i) {
        for (int iv = 0; iv < nvar; ++iv) {
            real_t q_prev = q[(i - 1) * nvar + iv];
            real_t q_curr = q[i * nvar + iv];
            real_t q_next = q[(i + 1) * nvar + iv];

            real_t dl = q_curr - q_prev;
            real_t dr = q_next - q_curr;

            if (dl * dr <= 0.0) {
                dq[i * nvar + iv] = 0.0;
            } else {
                real_t sgn = (dl >= 0.0) ? 1.0 : -1.0;
                dq[i * nvar + iv] = sgn * std::min(std::abs(dl), std::abs(dr));
            }
        }
    }
    // Boundary slopes (zero-gradient)
    for (int iv = 0; iv < nvar; ++iv) {
        dq[0 * nvar + iv] = 0.0;
        dq[(nx - 1) * nvar + iv] = 0.0;
    }
}

void Muscl::reconstruct(const real_t* q, const real_t* dq, real_t* qL, real_t* qR, int nx, int nvar) {
    for (int i = 0; i < nx; ++i) {
        for (int iv = 0; iv < nvar; ++iv) {
            real_t q_val = q[i * nvar + iv];
            real_t dq_val = dq[i * nvar + iv];
            
            // Interfaces: q_i-1/2, q_i+1/2
            qL[i * nvar + iv] = q_val - 0.5 * dq_val;
            qR[i * nvar + iv] = q_val + 0.5 * dq_val;
        }
    }
}

void Muscl::predict(const real_t q[], const real_t dq[], const real_t s0[], 
                   real_t dt_dx, real_t qm[], real_t qp[], int nvar) {
    const real_t smallr = 1e-10;
    for (int iv = 0; iv < nvar; ++iv) {
        real_t src_term = s0[iv] * dt_dx * 0.5;
        qp[iv] = q[iv] - 0.5 * dq[iv] + src_term;
        qm[iv] = q[iv] + 0.5 * dq[iv] + src_term;
    }
    if (qp[0] < smallr) qp[0] = q[0];
    if (qm[0] < smallr) qm[0] = q[0];
}

} // namespace ramses
