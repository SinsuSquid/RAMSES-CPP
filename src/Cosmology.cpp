#include "ramses/Cosmology.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace ramses {

void Cosmology::solve_friedman(real_t omega_m, real_t omega_l, real_t omega_k, real_t aexp_ini) {
    omega_m_ = omega_m;
    omega_l_ = omega_l;
    omega_k_ = omega_k;

    const int n_frw = 1000;
    const real_t alpha = 1e-6;
    const real_t axp_min = std::min(aexp_ini, static_cast<real_t>(1e-6));

    aexp_frw.assign(n_frw + 1, 0.0);
    hexp_frw.assign(n_frw + 1, 0.0);
    tau_frw.assign(n_frw + 1, 0.0);
    t_frw.assign(n_frw + 1, 0.0);

    real_t axp_tau = 1.0;
    real_t axp_t_val = 1.0;
    real_t tau = 0.0;
    real_t t = 0.0;
    int nstep = 0;

    // First pass to count steps
    while (axp_tau >= axp_min || axp_t_val >= axp_min) {
        nstep++;
        real_t dtau = alpha * axp_tau / dadtau(axp_tau, omega_m, omega_l, omega_k);
        real_t axp_tau_pre = axp_tau - dadtau(axp_tau, omega_m, omega_l, omega_k) * dtau * 0.5;
        axp_tau -= dadtau(axp_tau_pre, omega_m, omega_l, omega_k) * dtau;
        tau -= dtau;

        real_t dt = alpha * axp_t_val / dadt(axp_t_val, omega_m, omega_l, omega_k);
        real_t axp_t_pre = axp_t_val - dadt(axp_t_val, omega_m, omega_l, omega_k) * dt * 0.5;
        axp_t_val -= dadt(axp_t_pre, omega_m, omega_l, omega_k) * dt;
        t -= dt;
    }

    int nskip = std::max(1, nstep / n_frw);

    // Second pass to fill tables
    axp_tau = 1.0;
    axp_t_val = 1.0;
    tau = 0.0;
    t = 0.0;
    int curr_step = 0;
    int nout = 0;

    aexp_frw[0] = axp_tau;
    hexp_frw[0] = dadtau(axp_tau, omega_m, omega_l, omega_k) / axp_tau;
    tau_frw[0] = tau;
    t_frw[0] = t;

    while ((axp_tau >= axp_min || axp_t_val >= axp_min) && nout < n_frw) {
        curr_step++;
        real_t dtau = alpha * axp_tau / dadtau(axp_tau, omega_m, omega_l, omega_k);
        real_t axp_tau_pre = axp_tau - dadtau(axp_tau, omega_m, omega_l, omega_k) * dtau * 0.5;
        axp_tau -= dadtau(axp_tau_pre, omega_m, omega_l, omega_k) * dtau;
        tau -= dtau;

        real_t dt = alpha * axp_t_val / dadt(axp_t_val, omega_m, omega_l, omega_k);
        real_t axp_t_pre = axp_t_val - dadt(axp_t_val, omega_m, omega_l, omega_k) * dt * 0.5;
        axp_t_val -= dadt(axp_t_pre, omega_m, omega_l, omega_k) * dt;
        t -= dt;

        if (curr_step % nskip == 0) {
            nout++;
            aexp_frw[nout] = axp_tau;
            hexp_frw[nout] = dadtau(axp_tau, omega_m, omega_l, omega_k) / axp_tau;
            tau_frw[nout] = tau;
            t_frw[nout] = t;
        }
    }
    // Fill remaining slots if any
    for (int i = nout + 1; i <= n_frw; ++i) {
        aexp_frw[i] = axp_tau;
        hexp_frw[i] = dadtau(axp_tau, omega_m, omega_l, omega_k) / axp_tau;
        tau_frw[i] = tau;
        t_frw[i] = t;
    }
}

void Cosmology::get_cosmo_params(real_t tau, real_t& aexp, real_t& hexp, real_t& texp) {
    if (tau_frw.empty()) return;

    // Binary search for tau
    auto it = std::lower_bound(tau_frw.rbegin(), tau_frw.rend(), tau);
    int i = std::distance(tau_frw.rbegin(), it);
    i = tau_frw.size() - 1 - i;

    if (i <= 0) {
        aexp = aexp_frw[0]; hexp = hexp_frw[0]; texp = t_frw[0];
    } else {
        real_t frac = (tau - tau_frw[i - 1]) / (tau_frw[i] - tau_frw[i - 1]);
        aexp = aexp_frw[i - 1] + frac * (aexp_frw[i] - aexp_frw[i - 1]);
        hexp = hexp_frw[i - 1] + frac * (hexp_frw[i] - hexp_frw[i - 1]);
        texp = t_frw[i - 1] + frac * (t_frw[i] - t_frw[i - 1]);
    }
}

real_t Cosmology::get_tau(real_t aexp_target) {
    if (aexp_frw.empty()) return 0.0;

    auto it = std::lower_bound(aexp_frw.rbegin(), aexp_frw.rend(), aexp_target);
    int i = std::distance(aexp_frw.rbegin(), it);
    i = aexp_frw.size() - 1 - i;

    if (i <= 0) return tau_frw[0];
    real_t frac = (aexp_target - aexp_frw[i - 1]) / (aexp_frw[i] - aexp_frw[i - 1]);
    return tau_frw[i - 1] + frac * (tau_frw[i] - tau_frw[i - 1]);
}

real_t Cosmology::fpeebl(real_t a) {
    real_t y = omega_m_ * (1.0 / a - 1.0) + omega_l_ * (a * a - 1.0) + 1.0;
    real_t fact = rombint(1e-6, a, 1e-6);
    return (omega_l_ * a * a - 0.5 * omega_m_ / a) / y - 1.0 + a * fy(a, omega_m_, omega_l_, omega_k_) / fact;
}

real_t Cosmology::d1a(real_t a) {
    real_t y = omega_m_ * (1.0 / a - 1.0) + omega_l_ * (a * a - 1.0) + 1.0;
    real_t y12 = std::sqrt(std::max(static_cast<real_t>(0.0), y));
    return y12 / a * rombint(1e-6, a, 1e-6);
}

real_t Cosmology::rombint(real_t a, real_t b, real_t tol) {
    const int maxiter = 16;
    const int maxj = 5;
    std::vector<real_t> g(100, 0.0);
    
    real_t h = 0.5 * (b - a);
    real_t gmax = h * (fy(a, omega_m_, omega_l_, omega_k_) + fy(b, omega_m_, omega_l_, omega_k_));
    g[0] = gmax;
    int nint = 1;
    real_t error = 1e20;
    
    for (int i = 1; i <= maxiter; ++i) {
        real_t g0 = 0.0;
        for (int k = 1; k <= nint; ++k) {
            g0 += fy(a + (k + k - 1) * h, omega_m_, omega_l_, omega_k_);
        }
        g0 = 0.5 * g[0] + h * g0;
        h = 0.5 * h;
        nint *= 2;
        int jmax = std::min(i, maxj);
        real_t fourj = 1.0;
        
        for (int j = 1; j <= jmax; ++j) {
            fourj *= 4.0;
            real_t g1 = g0 + (g0 - g[j - 1]) / (fourj - 1.0);
            g[j - 1] = g0;
            g0 = g1;
        }
        
        if (std::abs(g0) > tol) {
            error = 1.0 - gmax / g0;
        } else {
            error = gmax;
        }
        gmax = g0;
        g[jmax] = g0;
        if (i > 5 && std::abs(error) < tol) return g0;
    }
    return g[0];
}

real_t Cosmology::fy(real_t a, real_t omega_m, real_t omega_l, real_t omega_k) {
    real_t y = omega_m * (1.0 / a - 1.0) + omega_l * (a * a - 1.0) + 1.0;
    return 1.0 / std::pow(y, 1.5);
}

real_t Cosmology::dadtau(real_t a, real_t omega_m, real_t omega_l, real_t omega_k) {
    real_t d = a * a * a * (omega_m + omega_l * a * a * a + omega_k * a);
    return std::sqrt(std::max(static_cast<real_t>(0.0), d));
}

real_t Cosmology::dadt(real_t a, real_t omega_m, real_t omega_l, real_t omega_k) {
    real_t d = (1.0 / a) * (omega_m + omega_l * a * a * a + omega_k * a);
    return std::sqrt(std::max(static_cast<real_t>(0.0), d));
}

} // namespace ramses
