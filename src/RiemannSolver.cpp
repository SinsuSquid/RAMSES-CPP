#include "ramses/RiemannSolver.hpp"
#include <cmath>
#include <algorithm>

namespace ramses {

void RiemannSolver::solve_llf(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma) {
    real_t fl[20], fr[20];
    compute_flux(ql, fl, gamma);
    compute_flux(qr, fr, gamma);
    real_t ul[20], ur[20];
    prim_to_cons(ql, ul, gamma);
    prim_to_cons(qr, ur, gamma);
    
    int ipress = NDIM + 1;
    real_t al = std::sqrt(gamma * ql[ipress] / std::max(ql[0], 1e-10));
    real_t ar = std::sqrt(gamma * qr[ipress] / std::max(qr[0], 1e-10));
    real_t a_max = std::max(std::abs(ql[1]) + al, std::abs(qr[1]) + ar);
    for (int i = 0; i < NDIM + 2; ++i) {
        flux[i] = 0.5 * (fl[i] + fr[i]) - 0.5 * a_max * (ur[i] - ul[i]);
    }
}

void RiemannSolver::solve_hll(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma) {
    real_t fl[20], fr[20];
    compute_flux(ql, fl, gamma);
    compute_flux(qr, fr, gamma);
    real_t ul_vec[20], ur_vec[20];
    prim_to_cons(ql, ul_vec, gamma);
    prim_to_cons(qr, ur_vec, gamma);
    
    int ipress = NDIM + 1;
    real_t al = std::sqrt(gamma * std::max(ql[ipress], 0.0) / std::max(ql[0], 1e-10));
    real_t ar = std::sqrt(gamma * std::max(qr[ipress], 0.0) / std::max(qr[0], 1e-10));
    
    real_t sl = std::min(0.0, std::min(ql[1] - al, qr[1] - ar));
    real_t sr = std::max(0.0, std::max(ql[1] + al, qr[1] + ar));
    
    if (sr - sl < 1e-20) {
        for (int i = 0; i < NDIM + 2; ++i) flux[i] = 0.5 * (fl[i] + fr[i]);
    } else {
        for (int i = 0; i < NDIM + 2; ++i) {
            flux[i] = (sr * fl[i] - sl * fr[i] + sr * sl * (ur_vec[i] - ul_vec[i])) / (sr - sl);
        }
    }
}

void RiemannSolver::solve_hllc(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma) {
    const real_t smallr = 1e-10;
    const real_t smallp = 1e-10;

    int ipress = NDIM + 1;
    real_t rl = std::max(ql[0], smallr);
    real_t ul = ql[1];
    real_t pl = std::max(ql[ipress], rl * smallp);

    real_t rr = std::max(qr[0], smallr);
    real_t ur = qr[1];
    real_t pr = std::max(qr[ipress], rr * smallp);

    // Sound speeds
    real_t cl = std::sqrt(gamma * pl / rl);
    real_t cr = std::sqrt(gamma * pr / rr);

    // Einfeldt wave speed estimates for robustness
    real_t sl = std::min(ul - cl, ur - cr);
    real_t sr = std::max(ul + cl, ur + cr);
    
    // Ensure the wave speeds are symmetric enough for zero velocity
    if (std::abs(ul) < 1e-10 && std::abs(ur) < 1e-10) {
        real_t cmax = std::max(cl, cr);
        sl = -cmax; sr = cmax;
    }

    if (sl > 0.0) {
        compute_flux(ql, flux, gamma);
    } else if (sr < 0.0) {
        compute_flux(qr, flux, gamma);
    } else {
        // Star state
        real_t rcl = rl * (ul - sl);
        real_t rcr = rr * (sr - ur);
        real_t ustar = (rcr * ur + rcl * ul + (pl - pr)) / (rcr + rcl);
        real_t pstar = (rcr * pl + rcl * pr + rcl * rcr * (ul - ur)) / (rcr + rcl);

        if (ustar > 0.0) {
            real_t rstarl = rl * (sl - ul) / (std::min(sl - ustar, -1e-10));
            real_t qstarl[20] = {0};
            qstarl[0] = rstarl; qstarl[1] = ustar;
            for(int i=2; i<=NDIM; ++i) qstarl[i] = ql[i];
            qstarl[ipress] = pstar;
            
            real_t fl_arr[20], ul_vec[20], ustarl_vec[20];
            compute_flux(ql, fl_arr, gamma);
            prim_to_cons(ql, ul_vec, gamma);
            prim_to_cons(qstarl, ustarl_vec, gamma);
            for (int i = 0; i < NDIM + 2; ++i) flux[i] = fl_arr[i] + sl * (ustarl_vec[i] - ul_vec[i]);
        } else {
            real_t rstarr = rr * (sr - ur) / (std::max(sr - ustar, 1e-10));
            real_t qstarr[20] = {0};
            qstarr[0] = rstarr; qstarr[1] = ustar;
            for(int i=2; i<=NDIM; ++i) qstarr[i] = qr[i];
            qstarr[ipress] = pstar;

            real_t fr_arr[20], ur_vec[20], ustarr_vec[20];
            compute_flux(qr, fr_arr, gamma);
            prim_to_cons(qr, ur_vec, gamma);
            prim_to_cons(qstarr, ustarr_vec, gamma);
            for (int i = 0; i < NDIM + 2; ++i) flux[i] = fr_arr[i] + sr * (ustarr_vec[i] - ur_vec[i]);
        }
    }
}

void RiemannSolver::prim_to_cons(const real_t q[], real_t u[], real_t gamma) {
    u[0] = q[0];
    real_t v2 = 0;
    for (int i = 1; i <= NDIM; ++i) {
        u[i] = q[0] * q[i];
        v2 += q[i] * q[i];
    }
    int ipress = NDIM + 1;
    real_t e_kin = 0.5 * q[0] * v2;
    real_t e_int = q[ipress] / (gamma - 1.0);
    u[ipress] = e_kin + e_int;
}

void RiemannSolver::compute_flux(const real_t q[], real_t f[], real_t gamma) {
    real_t rho = q[0]; real_t u = q[1]; 
    int ipress = NDIM + 1;
    real_t p = q[ipress];
    
    real_t v2 = 0;
    for(int i=1; i<=NDIM; ++i) v2 += q[i]*q[i];
    
    real_t e_kin = 0.5 * rho * v2;
    real_t e_int = p / (gamma - 1.0);
    real_t E = e_kin + e_int;
    
    f[0] = rho * u;
    f[1] = rho * u * u + p;
    for(int i=2; i<=NDIM; ++i) f[i] = rho * u * q[i];
    f[ipress] = (E + p) * u;
}

} // namespace ramses
