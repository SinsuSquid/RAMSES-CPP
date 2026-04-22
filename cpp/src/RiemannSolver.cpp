#include "ramses/RiemannSolver.hpp"
#include <cmath>
#include <algorithm>

namespace ramses {

void RiemannSolver::solve_llf(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma) {
    real_t fl[5], fr[5];
    compute_flux(ql, fl, gamma);
    compute_flux(qr, fr, gamma);
    real_t ul[5], ur[5];
    prim_to_cons(ql, ul, gamma);
    prim_to_cons(qr, ur, gamma);
    real_t al = std::sqrt(gamma * ql[4] / ql[0]);
    real_t ar = std::sqrt(gamma * qr[4] / qr[0]);
    real_t a_max = std::max(std::abs(ql[1]) + al, std::abs(qr[1]) + ar);
    for (int i = 0; i < 5; ++i) {
        flux[i] = 0.5 * (fl[i] + fr[i]) - 0.5 * a_max * (ur[i] - ul[i]);
    }
}

void RiemannSolver::solve_hllc(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma) {
    const real_t smallr = 1e-10;
    const real_t smallp = 1e-10;

    real_t rl = std::max(ql[0], smallr);
    real_t ul = ql[1];
    real_t pl = std::max(ql[4], rl * smallp);

    real_t rr = std::max(qr[0], smallr);
    real_t ur = qr[1];
    real_t pr = std::max(qr[4], rr * smallp);

    // Sound speeds
    real_t cl = std::sqrt(gamma * pl / rl);
    real_t cr = std::sqrt(gamma * pr / rr);

    // Wave speed estimates (Davis)
    real_t sl = std::min(ul, ur) - std::max(cl, cr);
    real_t sr = std::max(ul, ur) + std::max(cl, cr);

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
            real_t rstarl = rl * (sl - ul) / (sl - ustar);
            real_t qstarl[5] = {rstarl, ustar, ql[2], ql[3], pstar};
            real_t fl[5], ul_vec[5], ustarl_vec[5];
            compute_flux(ql, fl, gamma);
            prim_to_cons(ql, ul_vec, gamma);
            prim_to_cons(qstarl, ustarl_vec, gamma);
            for (int i = 0; i < 5; ++i) flux[i] = fl[i] + sl * (ustarl_vec[i] - ul_vec[i]);
        } else {
            real_t rstarr = rr * (sr - ur) / (sr - ustar);
            real_t qstarr[5] = {rstarr, ustar, qr[2], qr[3], pstar};
            real_t fr[5], ur_vec[5], ustarr_vec[5];
            compute_flux(qr, fr, gamma);
            prim_to_cons(qr, ur_vec, gamma);
            prim_to_cons(qstarr, ustarr_vec, gamma);
            for (int i = 0; i < 5; ++i) flux[i] = fr[i] + sr * (ustarr_vec[i] - ur_vec[i]);
        }
    }
}

void RiemannSolver::prim_to_cons(const real_t q[], real_t u[], real_t gamma) {
    u[0] = q[0];
    u[1] = q[0] * q[1];
    u[2] = q[0] * q[2];
    u[3] = q[0] * q[3];
    real_t e_kin = 0.5 * q[0] * (q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    real_t e_int = q[4] / (gamma - 1.0);
    u[4] = e_kin + e_int;
}

void RiemannSolver::compute_flux(const real_t q[], real_t f[], real_t gamma) {
    real_t rho = q[0]; real_t u = q[1]; real_t v = q[2]; real_t w = q[3]; real_t p = q[4];
    real_t e_kin = 0.5 * rho * (u*u + v*v + w*w);
    real_t e_int = p / (gamma - 1.0);
    real_t E = e_kin + e_int;
    f[0] = rho * u;
    f[1] = rho * u * u + p;
    f[2] = rho * u * v;
    f[3] = rho * u * w;
    f[4] = (E + p) * u;
}

} // namespace ramses
