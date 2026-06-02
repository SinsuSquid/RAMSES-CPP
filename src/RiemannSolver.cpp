#include "ramses/RiemannSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace ramses {

real_t RiemannSolver::get_cs2(real_t d, real_t p, real_t gamma, const real_t q[], int nener, const std::vector<real_t>& gamma_rad) {
    if (params::barotropic_eos) {
        real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
        real_t v2_unit = params::units_velocity * params::units_velocity;
        if (params::barotropic_eos_form == "isothermal") return T2 / v2_unit;
        if (params::barotropic_eos_form == "polytrope") return params::polytrope_index * T2 * std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0) / v2_unit;
        if (params::barotropic_eos_form == "double_polytrope") {
            real_t fac = std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
            return T2 * (1.0 + params::polytrope_index * fac) / v2_unit;
        }
    }
    real_t cs2 = gamma * p;
    if (q != nullptr) {
        for (int ie = 0; ie < nener; ++ie) {
            cs2 += gamma_rad[ie] * q[NDIM + 2 + ie];
        }
    }
    return std::clamp(cs2 / std::max(d, 1e-10), 1e-2, 1e12);
}

void RiemannSolver::solve_llf(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, int nener, const std::vector<real_t>& gamma_rad) {
    real_t fl[20], fr[20];
    compute_flux(ql, fl, gamma, nener, gamma_rad);
    compute_flux(qr, fr, gamma, nener, gamma_rad);
    real_t ul[20], ur[20];
    prim_to_cons(ql, ul, gamma, nener, gamma_rad);
    prim_to_cons(qr, ur, gamma, nener, gamma_rad);

    int ipress = NDIM + 1;
    real_t al = std::sqrt(get_cs2(ql[0], ql[ipress], gamma, ql, nener, gamma_rad));
    real_t ar = std::sqrt(get_cs2(qr[0], qr[ipress], gamma, qr, nener, gamma_rad));
    real_t a_max = std::max(std::abs(ql[1]) + al, std::abs(qr[1]) + ar);
    for (int i = 0; i < NDIM + 2; ++i) {
        flux[i] = 0.5 * (fl[i] + fr[i]) - 0.5 * a_max * (ur[i] - ul[i]);
    }
}

void RiemannSolver::solve_hll(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, int nener, const std::vector<real_t>& gamma_rad) {
    real_t fl[20], fr[20];
    compute_flux(ql, fl, gamma, nener, gamma_rad);
    compute_flux(qr, fr, gamma, nener, gamma_rad);
    real_t ul_vec[20], ur_vec[20];
    prim_to_cons(ql, ul_vec, gamma, nener, gamma_rad);
    prim_to_cons(qr, ur_vec, gamma, nener, gamma_rad);

    int ipress = NDIM + 1;
    real_t al = std::sqrt(get_cs2(ql[0], ql[ipress], gamma, ql, nener, gamma_rad));
    real_t ar = std::sqrt(get_cs2(qr[0], qr[ipress], gamma, qr, nener, gamma_rad));

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

void RiemannSolver::solve_hllc(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, int nener, const std::vector<real_t>& gamma_rad) {
    const real_t smallc = 1e-10;
    const real_t smallr = smallc;
    const real_t smallp = smallc * smallc / gamma;

    int ipress = NDIM + 1;
    real_t rl = std::max(ql[0], smallr);
    real_t ul = ql[1];
    real_t pl = std::max(ql[ipress], rl * smallp);
    for (int ie = 0; ie < nener; ++ie) {
        pl += ql[ipress + 1 + ie];
    }

    real_t rr = std::max(qr[0], smallr);
    real_t ur = qr[1];
    real_t pr = std::max(qr[ipress], rr * smallp);
    for (int ie = 0; ie < nener; ++ie) {
        pr += qr[ipress + 1 + ie];
    }

    real_t cl = std::sqrt(get_cs2(rl, ql[ipress], gamma, ql, nener, gamma_rad));
    real_t cr = std::sqrt(get_cs2(rr, qr[ipress], gamma, qr, nener, gamma_rad));

    real_t sl = std::min(ul, ur) - std::max(cl, cr);
    real_t sr = std::max(ul, ur) + std::max(cl, cr);

    if (sl > 0.0) {
        compute_flux(ql, flux, gamma, nener, gamma_rad);
    } else if (sr < 0.0) {
        compute_flux(qr, flux, gamma, nener, gamma_rad);
    } else {
        real_t rcl = rl * (ul - sl);
        real_t rcr = rr * (sr - ur);
        real_t div = rcr + rcl;
        if (div < 1e-5) {
            real_t fl_fallback[20], fr_fallback[20];
            compute_flux(ql, fl_fallback, gamma, nener, gamma_rad);
            compute_flux(qr, fr_fallback, gamma, nener, gamma_rad);
            for (int i = 0; i < NDIM + 2; ++i) flux[i] = 0.5 * (fl_fallback[i] + fr_fallback[i]);
            return;
        }
        real_t ustar = (rcr * ur + rcl * ul + (pl - pr)) / div;
        ustar = std::max(-1e6, std::min(1e6, ustar));
        real_t pstar = (rcr * pl + rcl * pr + rcl * rcr * (ul - ur)) / div;
        pstar = std::max(pstar, 0.0);

        if (ustar > 0.0) {
            real_t rstarl = rl * (sl - ul) / (std::min(sl - ustar, -1e-10));
            real_t qstarl[20] = {0};
            qstarl[0] = rstarl; qstarl[1] = ustar;
            for(int i=2; i<=NDIM; ++i) qstarl[i] = ql[i];
            
            real_t sum_prad_star = 0.0;
            for (int ie = 0; ie < nener; ++ie) {
                real_t prad_star = ql[ipress + 1 + ie] * (sl - ul) / (std::min(sl - ustar, -1e-10));
                qstarl[ipress + 1 + ie] = prad_star;
                sum_prad_star += prad_star;
            }
            qstarl[ipress] = std::max(pstar - sum_prad_star, rstarl * smallp);

            real_t fl_arr[20], ul_vec[20], ustarl_vec[20];
            compute_flux(ql, fl_arr, gamma, nener, gamma_rad);
            prim_to_cons(ql, ul_vec, gamma, nener, gamma_rad);
            prim_to_cons(qstarl, ustarl_vec, gamma, nener, gamma_rad);
            ustarl_vec[ipress] = ((sl - ul) * ul_vec[ipress] - pl * ul + pstar * ustar) / (std::min(sl - ustar, -1e-10));
            for (int i = 0; i < NDIM + 2; ++i) flux[i] = fl_arr[i] + sl * (ustarl_vec[i] - ul_vec[i]);
        } else {
            real_t rstarr = rr * (sr - ur) / (std::max(sr - ustar, 1e-10));
            real_t qstarr[20] = {0};
            qstarr[0] = rstarr; qstarr[1] = ustar;
            for(int i=2; i<=NDIM; ++i) qstarr[i] = qr[i];

            real_t sum_prad_star = 0.0;
            for (int ie = 0; ie < nener; ++ie) {
                real_t prad_star = qr[ipress + 1 + ie] * (sr - ur) / (std::max(sr - ustar, 1e-10));
                qstarr[ipress + 1 + ie] = prad_star;
                sum_prad_star += prad_star;
            }
            qstarr[ipress] = std::max(pstar - sum_prad_star, rstarr * smallp);

            real_t fr_arr[20], ur_vec[20], ustarr_vec[20];
            compute_flux(qr, fr_arr, gamma, nener, gamma_rad);
            prim_to_cons(qr, ur_vec, gamma, nener, gamma_rad);
            prim_to_cons(qstarr, ustarr_vec, gamma, nener, gamma_rad);
            ustarr_vec[ipress] = ((sr - ur) * ur_vec[ipress] - pr * ur + pstar * ustar) / (std::max(sr - ustar, 1e-10));
            for (int i = 0; i < NDIM + 2; ++i) flux[i] = fr_arr[i] + sr * (ustarr_vec[i] - ur_vec[i]);
        }
    }
}

void RiemannSolver::solve_godunov_nr(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma) {
    const real_t smallc = 1e-10;
    const real_t smallr = smallc;
    const real_t smallp = smallc * smallc / gamma;
    const int niter = 10;
    const real_t gamma6 = (gamma + 1.0) / (2.0 * gamma);

    int ipress = NDIM + 1;
    real_t rl = std::max(ql[0], smallr);
    real_t ul = ql[1];
    real_t pl = std::max(ql[ipress], rl * smallp);

    real_t rr = std::max(qr[0], smallr);
    real_t ur = qr[1];
    real_t pr = std::max(qr[ipress], rr * smallp);

    real_t cl_sq = gamma * pl * rl;
    real_t cr_sq = gamma * pr * rr;
    real_t wl = std::sqrt(std::max(cl_sq, 0.0));
    real_t wr = std::sqrt(std::max(cr_sq, 0.0));

    real_t pstar = (wr * pl + wl * pr + wl * wr * (ul - ur)) / (wl + wr + 1e-20);
    pstar = std::max(pstar, 0.0);
    real_t pold = pstar;

    // Newton-Raphson iterations to find pstar
    for (int iter = 0; iter < niter; ++iter) {
        real_t wwl_sq = cl_sq * (1.0 + gamma6 * (pold - pl) / (pl + 1e-20));
        real_t wwr_sq = cr_sq * (1.0 + gamma6 * (pold - pr) / (pr + 1e-20));
        real_t wwl = std::sqrt(std::max(wwl_sq, 0.0));
        real_t wwr = std::sqrt(std::max(wwr_sq, 0.0));

        real_t ql_nr = (2.0 * std::pow(wwl, 3)) / (wwl * wwl + cl_sq + 1e-20);
        real_t qr_nr = (2.0 * std::pow(wwr, 3)) / (wwr * wwr + cr_sq + 1e-20);
        real_t usl = ul - (pold - pl) / (wwl + 1e-20);
        real_t usr = ur + (pold - pr) / (wwr + 1e-20);

        real_t delp = (qr_nr * ql_nr / (qr_nr + ql_nr + 1e-20)) * (usl - usr);
        delp = std::max(delp, -pold);
        pold = pold + delp;
        if (std::abs(delp / (pold + 1e-10)) < 1e-6) break;
    }
    pstar = pold;


    real_t wwl_f = std::sqrt(std::max(cl_sq * (1.0 + gamma6 * (pstar - pl) / (pl + 1e-20)), 0.0));
    real_t wwr_f = std::sqrt(std::max(cr_sq * (1.0 + gamma6 * (pstar - pr) / (pr + 1e-20)), 0.0));
    real_t ustar = 0.5 * (ul + (pl - pstar) / (wwl_f + 1e-20) + ur - (pr - pstar) / (wwr_f + 1e-20));

    real_t ro, uo, po;
    real_t sgnm = (ustar > 0) ? 1.0 : -1.0;

    if (sgnm == 1.0) {
        ro = rl; uo = ul; po = pl;
    } else {
        ro = rr; uo = ur; po = pr;
    }

    real_t co = std::sqrt(std::max(gamma * po / (ro + 1e-20), 1e-20));
    real_t rstar_val;
    if (pstar >= po) {
        real_t cl_unperturbed = (sgnm == 1.0) ? cl_sq : cr_sq;
        rstar_val = ro * ((gamma + 1.0) * pstar + (gamma - 1.0) * po) / ((gamma - 1.0) * pstar + (gamma + 1.0) * po);
    } else {
        rstar_val = ro * std::pow(pstar / (po + 1e-20), 1.0 / gamma);
    }
    rstar_val = std::max(rstar_val, smallr);

    real_t cstar = std::sqrt(std::max(gamma * pstar / (rstar_val + 1e-20), 1e-20));
    real_t spout = co - sgnm * uo;
    real_t spin = cstar - sgnm * ustar;

    if (pstar >= po) {
        real_t wo = (sgnm == 1.0) ? wwl_f : wwr_f;
        real_t ushock = wo / (ro + 1e-20) - sgnm * uo;
        spout = ushock;
        spin = ushock;
    }

    // Use HLLC-style wave speed estimates for robust sampling
    real_t cl_hllc = std::sqrt(get_cs2(rl, pl, gamma));
    real_t cr_hllc = std::sqrt(get_cs2(rr, pr, gamma));
    real_t sl_hllc = std::min(ul, ur) - std::max(cl_hllc, cr_hllc);
    real_t sr_hllc = std::max(ul, ur) + std::max(cl_hllc, cr_hllc);

    if (sl_hllc > 0.0) {
        ro = rl; uo = ul; po = pl;
    } else if (sr_hllc < 0.0) {
        ro = rr; uo = ur; po = pr;
    } else {
        ro = rstar_val; uo = ustar; po = pstar;
    }

    if (std::abs(pl - pr) > 1e-3) {
        static int count_nr = 0; if(count_nr++ < 5) {
            std::cout << "[NR DEBUG] pstar=" << pstar << " ro=" << ro << " uo=" << uo << " po=" << po << " spout=" << spout << " spin=" << spin << std::endl;
        }
    }

    real_t e_int = po / (gamma - 1.0);
    real_t e_kin = 0.5 * ro * uo * uo;
    real_t E = e_int + e_kin;

    flux[0] = ro * uo;
    flux[1] = ro * uo * uo + po;
    for(int i=2; i<=NDIM; ++i) flux[i] = 0.0;
    flux[ipress] = (E + po) * uo;
}

void RiemannSolver::prim_to_cons(const real_t q[], real_t u[], real_t gamma, int nener, const std::vector<real_t>& gamma_rad) {
    u[0] = q[0];
    real_t v2 = 0;
    for (int i = 1; i <= NDIM; ++i) {
        u[i] = q[0] * q[i];
        v2 += q[i] * q[i];
    }
    int ipress = NDIM + 1;
    real_t e_kin = 0.5 * q[0] * v2;
    real_t e_int = q[ipress] / (gamma - 1.0);
    real_t e_nonthermal = 0.0;
    for (int ie = 0; ie < nener; ++ie) {
        e_nonthermal += q[ipress + 1 + ie] / (gamma_rad[ie] - 1.0);
    }
    if (params::barotropic_eos && params::barotropic_eos_form == "isothermal") {
    }
    u[ipress] = e_kin + e_int + e_nonthermal;
}

void RiemannSolver::compute_flux(const real_t q[], real_t f[], real_t gamma, int nener, const std::vector<real_t>& gamma_rad) {
    real_t rho = q[0]; real_t u = q[1];
    int ipress = NDIM + 1;
    real_t p = q[ipress];

    real_t v2 = 0;
    for(int i=1; i<=NDIM; ++i) v2 += q[i]*q[i];

    real_t e_kin = 0.5 * rho * v2;
    real_t e_int = p / (gamma - 1.0);
    real_t e_nonthermal = 0.0;
    real_t p_tot = p;
    for (int ie = 0; ie < nener; ++ie) {
        real_t prad = q[ipress + 1 + ie];
        p_tot += prad;
        e_nonthermal += prad / (gamma_rad[ie] - 1.0);
    }
    real_t E = e_kin + e_int + e_nonthermal;

    f[0] = rho * u;
    f[1] = rho * u * u + p_tot;
    for(int i=2; i<=NDIM; ++i) f[i] = rho * u * q[i];
    f[ipress] = (E + p_tot) * u;
}

} // namespace ramses
