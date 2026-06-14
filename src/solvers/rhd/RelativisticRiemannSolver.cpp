#include "ramses/solvers/rhd/RelativisticRiemannSolver.hpp"
#include <cmath>
#include <algorithm>

namespace ramses {

// ============================================================================
// 1. Local Lax-Friedrichs (LLF) Riemann Solver
// ============================================================================
void RelativisticRiemannSolver::solve_llf(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, const std::string& eos) {
    // Left and right flux vectors
    real_t fl[20], fr[20];
    // Left and right conserved variable vectors
    real_t ul[20], ur[20];
    // Left and right fast wave speeds (sound speeds)
    real_t cfl, cfr;

    // Compute conserved variables and fluxes for both left and right states
    find_rhd_flux(ql, ul, fl, gamma, eos);
    find_rhd_flux(qr, ur, fr, gamma, eos);

    // Compute fast magnetosonic/acoustic speeds
    find_speed_fast(ql, cfl, gamma, eos);
    find_speed_fast(qr, cfr, gamma, eos);

    // Estimate the maximum wave propagation speed:
    // s_max = max( |v_n,L| + c_L, |v_n,R| + c_R )
    // Under C++ layout:
    // q[0]: density, q[1]: normal velocity, q[2]/q[3]: transverse velocities, q[4]: pressure
    real_t s_max = std::max(std::abs(ql[1]) + cfl, std::abs(qr[1]) + cfr);

    // Local Lax-Friedrichs flux formula:
    // F_LLF = 0.5 * (F_L + F_R) - 0.5 * s_max * (U_R - U_L)
    for (int i = 0; i < 20; ++i) {
        flux[i] = 0.5 * (fl[i] + fr[i]) - 0.5 * s_max * (ur[i] - ul[i]);
    }
}

// ============================================================================
// 2. HLL (Harten-Lax-van Leer) Riemann Solver
// ============================================================================
void RelativisticRiemannSolver::solve_hll(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, const std::string& eos) {
    // Left and right flux vectors
    real_t fl[20], fr[20];
    // Left and right conserved variable vectors
    real_t ul[20], ur[20];
    // Left and right fast wave speeds
    real_t cfl, cfr;

    // Compute conserved variables and fluxes for both left and right states
    find_rhd_flux(ql, ul, fl, gamma, eos);
    find_rhd_flux(qr, ur, fr, gamma, eos);

    // Compute sound speeds
    find_speed_fast(ql, cfl, gamma, eos);
    find_speed_fast(qr, cfr, gamma, eos);

    // Under C++ layout:
    // q[0]: density, q[1]: normal velocity, q[2]: transverse velocity 1, q[3]: transverse velocity 2, q[4]: pressure
    real_t vnl = ql[1];
    real_t vnr = qr[1];

    // Compute velocity squares (v^2)
    real_t vlsq = ql[1]*ql[1] + ql[2]*ql[2] + ql[3]*ql[3];
    real_t vrsq = qr[1]*qr[1] + qr[2]*qr[2] + qr[3]*qr[3];

    // Compute transverse velocity squares (vt^2)
    real_t vplsq = ql[2]*ql[2] + ql[3]*ql[3];
    real_t vprsq = qr[2]*qr[2] + qr[3]*qr[3];

    // Lorentz factors (lor = 1 / sqrt(1 - v^2))
    real_t lorl = 1.0 / std::sqrt(1.0 - vlsq);
    real_t lorr = 1.0 / std::sqrt(1.0 - vrsq);

    // Relativistic wave speed correction factors:
    // factor = c_s * sqrt(1 - vn^2 - c_s^2 * vt^2) / lor
    // (Note: / lor is equivalent to multiplying by sqrt(1 - v^2))
    real_t factorl = cfl * std::sqrt(1.0 - vnl*vnl - cfl*cfl * vplsq) / lorl;
    real_t factorr = cfr * std::sqrt(1.0 - vnr*vnr - cfr*cfr * vprsq) / lorr;

    // Left and right propagation speeds of acoustic waves (lambda_L/R^\pm):
    // lambda^\pm = (vn * (1 - c_s^2) \pm factor) / (1 - c_s^2 * v^2)
    real_t lleftp = (vnl * (1.0 - cfl*cfl) + factorl) / (1.0 - cfl*cfl * vlsq);
    real_t lleftm = (vnl * (1.0 - cfl*cfl) - factorl) / (1.0 - cfl*cfl * vlsq);
    real_t lrightp = (vnr * (1.0 - cfr*cfr) + factorr) / (1.0 - cfr*cfr * vrsq);
    real_t lrightm = (vnr * (1.0 - cfr*cfr) - factorr) / (1.0 - cfr*cfr * vrsq);

    // HLL wave speed bounds:
    // S_R = max(0, lambda_L^+, lambda_R^+)
    // S_L = min(0, lambda_L^-, lambda_R^-)
    real_t sr = std::max({0.0, lleftp, lrightp});
    real_t sl = std::min({0.0, lleftm, lrightm});

    if (sl >= 0.0) {
        // Entirely supersonic flow to the right: return left flux
        for (int i = 0; i < 20; ++i) {
            flux[i] = fl[i];
        }
    } else if (sr <= 0.0) {
        // Entirely supersonic flow to the left: return right flux
        for (int i = 0; i < 20; ++i) {
            flux[i] = fr[i];
        }
    } else {
        // Subsonic region: HLL weighted average flux
        // F_HLL = (S_R * F_L - S_L * F_R + S_L * S_R * (U_R - U_L)) / (S_R - S_L)
        real_t inv_srs = 1.0 / (sr - sl);
        for (int i = 0; i < 20; ++i) {
            flux[i] = (sr * fl[i] - sl * fr[i] + sr * sl * (ur[i] - ul[i])) * inv_srs;
        }
    }
}

// ============================================================================
// 3. HLLC (Harten-Lax-van Leer-Contact) Riemann Solver
// ============================================================================
void RelativisticRiemannSolver::solve_hllc(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, const std::string& eos) {
    // Left and right flux vectors
    real_t fl[20], fr[20];
    // Left and right conserved variable vectors
    real_t ul[20], ur[20];
    // Left and right fast wave speeds
    real_t cfl, cfr;

    // Compute conserved variables and fluxes for both left and right states
    find_rhd_flux(ql, ul, fl, gamma, eos);
    find_rhd_flux(qr, ur, fr, gamma, eos);

    // Compute sound speeds
    find_speed_fast(ql, cfl, gamma, eos);
    find_speed_fast(qr, cfr, gamma, eos);

    // Under C++ layout:
    // q[0]: density, q[1]: normal velocity, q[2]: transverse velocity 1, q[3]: transverse velocity 2, q[4]: pressure
    real_t vnl = ql[1];
    real_t vnr = qr[1];

    // Compute velocity squares
    real_t vlsq = ql[1]*ql[1] + ql[2]*ql[2] + ql[3]*ql[3];
    real_t vrsq = qr[1]*qr[1] + qr[2]*qr[2] + qr[3]*qr[3];

    // Compute transverse velocity squares
    real_t vplsq = ql[2]*ql[2] + ql[3]*ql[3];
    real_t vprsq = qr[2]*qr[2] + qr[3]*qr[3];

    // Lorentz factors
    real_t lorl = 1.0 / std::sqrt(1.0 - vlsq);
    real_t lorr = 1.0 / std::sqrt(1.0 - vrsq);

    // Wave speed correction factors
    real_t factorl = cfl * std::sqrt(1.0 - vnl*vnl - cfl*cfl * vplsq) / lorl;
    real_t factorr = cfr * std::sqrt(1.0 - vnr*vnr - cfr*cfr * vprsq) / lorr;

    // Acoustic wave propagation speeds
    real_t lleftp = (vnl * (1.0 - cfl*cfl) + factorl) / (1.0 - cfl*cfl * vlsq);
    real_t lleftm = (vnl * (1.0 - cfl*cfl) - factorl) / (1.0 - cfl*cfl * vlsq);
    real_t lrightp = (vnr * (1.0 - cfr*cfr) + factorr) / (1.0 - cfr*cfr * vrsq);
    real_t lrightm = (vnr * (1.0 - cfr*cfr) - factorr) / (1.0 - cfr*cfr * vrsq);

    // HLL wave speed bounds
    real_t sr = std::max({lleftp, lrightp, 0.0});
    real_t sl = std::min({lleftm, lrightm, 0.0});

    if (sl >= 0.0) {
        // Entirely supersonic flow to the right
        for (int i = 0; i < 20; ++i) {
            flux[i] = fl[i];
        }
    } else if (sr <= 0.0) {
        // Entirely supersonic flow to the left
        for (int i = 0; i < 20; ++i) {
            flux[i] = fr[i];
        }
    } else {
        // HLL flux and conserved state in the star region
        real_t inv_srs = 1.0 / (sr - sl);
        real_t fhll[20], uhll[20];
        for (int i = 0; i < 20; ++i) {
            fhll[i] = (sr * fl[i] - sl * fr[i] + sr * sl * (ur[i] - ul[i])) * inv_srs;
            uhll[i] = (sr * ur[i] - sl * ul[i] + fl[i] - fr[i]) * inv_srs;
        }

        // Quadratic equation for contact/star wave speed sstar (Mignone Eq 18):
        // a * sstar^2 + b * sstar + c = 0
        // where:
        // a = fhll[E]        (fhll[1] is energy flux)
        // b = -(uhll[E] + fhll[mn])  (uhll[1] is energy density, fhll[2] is normal momentum flux)
        // c = uhll[mn]       (uhll[2] is normal momentum density)
        real_t a = fhll[1];
        real_t b = -(uhll[1] + fhll[2]);
        real_t c = uhll[2];
        real_t sstar;

        // Solve using the physical root (with minus sign)
        if (b > 0) {
            real_t quad = -0.5 * (b + std::sqrt(std::max(0.0, b*b - 4.0*a*c)));
            sstar = c / quad;
        } else {
            real_t quad = -0.5 * (b - std::sqrt(std::max(0.0, b*b - 4.0*a*c)));
            sstar = c / quad;
        }

        // Compute fluxes based on contact speed sstar
        if (sstar >= 0.0) {
            // Star region pressure: ps = -Fhll[E] * sstar + Fhll[mn]
            real_t ps = -fhll[1] * sstar + fhll[2];
            real_t den = 1.0 / (sl - sstar);

            // Compute star conserved variables on the left side:
            // U_L^* = ( (S_L - v_n) * U_L + (0, p^* s^* - p v_n, p^* - p, 0, 0)^T ) / (S_L - S^*)
            // Note indexing: 0 = density, 1 = energy, 2 = normal momentum, 3 = transverse 1, 4 = transverse 2
            real_t usl[20];
            usl[0] = ul[0] * (sl - ql[1]) * den;
            usl[2] = (ul[2] * (sl - ql[1]) + ps - ql[4]) * den;
            usl[3] = ul[3] * (sl - ql[1]) * den;
            usl[4] = ul[4] * (sl - ql[1]) * den;
            usl[1] = (ul[1] * (sl - ql[1]) + ps * sstar - ql[4] * ql[1]) * den;

            // Compute intercell HLLC flux: F_L^* = S_L * (U_L^* - U_L) + F_L
            for (int i = 0; i < 5; ++i) {
                flux[i] = sl * (usl[i] - ul[i]) + fl[i];
            }
            // Passive scalars flux
            for (int i = 5; i < 20; ++i) {
                flux[i] = (flux[0] > 0) ? flux[0] * ql[i] : flux[0] * qr[i];
            }
        } else {
            // Star region pressure on the right side
            real_t ps = -fhll[1] * sstar + fhll[2];
            real_t den = 1.0 / (sr - sstar);

            // Compute star conserved variables on the right side
            real_t usr[20];
            usr[0] = ur[0] * (sr - qr[1]) * den;
            usr[2] = (ur[2] * (sr - qr[1]) + ps - qr[4]) * den;
            usr[3] = ur[3] * (sr - qr[1]) * den;
            usr[4] = ur[4] * (sr - qr[1]) * den;
            usr[1] = (ur[1] * (sr - qr[1]) + ps * sstar - qr[4] * qr[1]) * den;

            // Compute intercell HLLC flux: F_R^* = S_R * (U_R^* - U_R) + F_R
            for (int i = 0; i < 5; ++i) {
                flux[i] = sr * (usr[i] - ur[i]) + fr[i];
            }
            // Passive scalars flux
            for (int i = 5; i < 20; ++i) {
                flux[i] = (flux[0] > 0) ? flux[0] * ql[i] : flux[0] * qr[i];
            }
        }
    }
}

// ============================================================================
// 4. Flux and Conserved Variable Calculations
// ============================================================================
void RelativisticRiemannSolver::find_rhd_flux(const real_t q[], real_t u[], real_t f[], real_t gamma, const std::string& eos) {
    // Under C++ layout:
    // q[0]: density (d), q[1]: normal velocity (vx), q[2]: transverse velocity 1 (vy)
    // q[3]: transverse velocity 2 (vz), q[4]: pressure (p)
    real_t d = q[0];
    real_t p = q[4];
    real_t vx = q[1];
    real_t vy = q[2];
    real_t vz = q[3];

    // Compute velocity squared and Lorentz factor:
    // v^2 = vx^2 + vy^2 + vz^2
    // lor = 1 / sqrt(1 - v^2)
    real_t v2 = vx*vx + vy*vy + vz*vz;
    real_t lor = 1.0 / std::sqrt(1.0 - v2);

    // Compute specific enthalpy (enth):
    // For Ideal gas EoS:  h = d + p * [gamma / (gamma - 1)]
    // For Taub-Mathews:   h = d * (2.5 * tau + 1.5 * sqrt(tau^2 + 4/9))  with tau = p/d
    real_t enth;
    if (eos == "TM") {
        real_t tau = p / d;
        enth = d * (2.5 * tau + 1.5 * std::sqrt(tau*tau + 4.0/9.0));
    } else {
        enth = d + p / (gamma - 1.0) + p;
    }

    // Conserved variable vector:
    // u[0] = D = d * lor                     (conserved density)
    // u[1] = E = enth * lor^2 - p            (conserved total energy density)
    // u[2] = Mx = enth * lor^2 * vx          (conserved normal momentum density)
    // u[3] = My = enth * lor^2 * vy          (conserved transverse momentum 1)
    // u[4] = Mz = enth * lor^2 * vz          (conserved transverse momentum 2)
    u[0] = d * lor;
    u[1] = enth * lor*lor - p;
    u[2] = enth * lor*lor * vx;
    u[3] = enth * lor*lor * vy;
    u[4] = enth * lor*lor * vz;

    // Passive scalars: u_s = D * q_s
    for (int i = 5; i < 20; ++i) {
        u[i] = d * q[i] * lor;
    }

    // Flux vector:
    // f[0] = D * vx = d * vx * lor
    // f[1] = E * vx + p * vx = enth * lor^2 * vx
    // f[2] = Mx * vx + p = enth * lor^2 * vx^2 + p
    // f[3] = My * vx = enth * lor^2 * vx * vy
    // f[4] = Mz * vx = enth * lor^2 * vx * vz
    f[0] = d * vx * lor;
    f[1] = enth * lor*lor * vx;
    f[2] = enth * lor*lor * vx*vx + p;
    f[3] = enth * lor*lor * vx * vy;
    f[4] = enth * lor*lor * vx * vz;

    // Passive scalars fluxes: f_s = u_s * vx = D * q_s * vx
    for (int i = 5; i < 20; ++i) {
        f[i] = d * vx * q[i] * lor;
    }
}

// ============================================================================
// 5. Sound Speed Calculation
// ============================================================================
void RelativisticRiemannSolver::find_speed_fast(const real_t q[], real_t& cs, real_t gamma, const std::string& eos) {
    // Under C++ layout:
    // q[0]: density (d), q[4]: pressure (p)
    real_t d = q[0];
    real_t p = q[4];

    // Compute sound speed cs:
    // cs^2 = gamma * p / (rho * h)
    if (eos == "TM") {
        real_t tau = p / d;
        real_t enth = 2.5 * tau + 1.5 * std::sqrt(tau*tau + 4.0/9.0);
        cs = std::sqrt(tau / (3.0 * enth) * (5.0 * enth - 8.0 * tau) / (enth - tau));
    } else {
        real_t enth = d + p / (gamma - 1.0) + p;
        cs = std::sqrt(gamma * p / enth);
    }
}

} // namespace ramses