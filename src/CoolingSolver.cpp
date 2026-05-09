#include "ramses/CoolingSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/MpiManager.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace ramses {

CoolingSolver::~CoolingSolver() {}

void CoolingSolver::apply_cooling(int ilevel, real_t dt) {
    if (!config_.get_bool("cooling_params", "cooling", false)) return;

    real_t gamma = config_.get_double("hydro_params", "gamma", 1.6666667);
    real_t mu = config_.get_double("cooling_params", "mu_gas", 1.4);
    real_t kB = 1.3806e-16;
    real_t mH = 1.67e-24;
    real_t dt_sec = dt * params::units_time;
    
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son[idc - 1] != 0) continue;

            real_t d = grid_.uold(idc, 1);
            real_t u = grid_.uold(idc, 2) / d;
            real_t v = grid_.uold(idc, 3) / d;
            real_t w = grid_.uold(idc, 4) / d;
            real_t etot = grid_.uold(idc, NDIM + 2);
            
            real_t ekin = 0.5 * d * (u*u + v*v + w*w);
            real_t emag = 0.0;
#ifdef MHD
            real_t A = 0.5 * (grid_.uold(idc, 6) + grid_.uold(idc, 9));
            real_t B = 0.5 * (grid_.uold(idc, 7) + grid_.uold(idc, 10));
            real_t C = 0.5 * (grid_.uold(idc, 8) + grid_.uold(idc, 11));
            emag = 0.5 * (A*A + B*B + C*C);
#endif
            real_t eint = etot - ekin - emag;
            if (eint < 0) eint = d * 1e-10; // Floor

            // Convert to physical units
            real_t nH_phys = (d * params::units_density) / (mu * mH);
            real_t T_phys = eint * params::units_pressure / (d * params::units_density) * (mu * mH / kB) * (gamma - 1.0);
            
            // Solve cooling (modifies T_phys)
            real_t T_new = T_phys;
            real_t time = 0.0;
            int iter = 0;
            real_t alpha_ct = nH_phys * kB / (gamma - 1.0);
            real_t dt_sub = 0.1 * dt_sec;

            while (time < dt_sec && iter < 100) {
                real_t T_old = T_new;
                real_t rate, rate2;
                if (T_new < 10035.0) {
                    rate = cooling_low(T_new, nH_phys);
                    rate2 = cooling_low(T_new * 1.00001, nH_phys);
                } else {
                    rate = cooling_high(T_new, nH_phys);
                    rate2 = cooling_high(T_new * 1.00001, nH_phys);
                }
                
                real_t dratedT = (rate2 - rate) / (T_new * 0.00001);
                real_t dTemp = rate / (alpha_ct / dt_sub - dratedT);
                
                // Limit temperature change to 20%
                if (std::abs(dTemp / T_new) > 0.2) dTemp = 0.2 * T_new * (dTemp > 0 ? 1.0 : -1.0);
                
                T_new = std::max(10.0, T_old + dTemp);
                time += dt_sub;
                dt_sub = std::min(dt_sec - time, dt_sub * 1.2);
                iter++;
            }

            // Convert back to code units
            real_t eint_new = T_new / (params::units_pressure / (d * params::units_density) * (mu * mH / kB) * (gamma - 1.0));
            real_t etot_new = eint_new + ekin + emag;
            
            grid_.uold(idc, NDIM + 2) = etot_new;
            grid_.unew(idc, NDIM + 2) = etot_new;
        }
        igrid = grid_.next[igrid - 1];
    }
}

real_t CoolingSolver::cooling_high(real_t T, real_t n) {
    real_t logT = std::log10(std::max(T, 1.0));
    real_t cold;
    if (logT < 4.0) cold = 0.1343*std::pow(logT,3) - 1.3906*std::pow(logT,2) + 5.1554*logT - 31.967;
    else if (logT < 4.25) cold = 12.64*logT - 75.56;
    else if (logT < 4.35) cold = -0.3*logT - 20.565;
    else if (logT < 4.9)  cold = 1.745*logT - 29.463;
    else if (logT < 5.4)  cold = -20.9125;
    else if (logT < 5.9)  cold = -1.795*logT - 11.219;
    else if (logT < 6.2)  cold = -21.8095;
    else if (logT < 6.7)  cold = -1.261*logT - 13.991;
    else cold = -22.44;
    
    return -std::pow(10.0, cold) * n * n;
}

real_t CoolingSolver::cooling_low(real_t T, real_t n) {
    real_t ne = 2.4e-3 * std::pow(T/100.0, 0.25) / 0.5;
    real_t x = std::clamp(ne / n, 3.5e-4 * 0.4, 0.1);
    
    real_t cold_cII = 92.0 * 1.38e-16 * 2.0 * (2.8e-7 * std::pow(T/100.0, -0.5) * x + 8e-10 * std::pow(T/100.0, 0.07)) * 3.5e-4 * 0.4 * std::exp(-92.0/T);
    real_t cold_o = 1e-26 * std::sqrt(T) * (24.0 * std::exp(-228.0/T) + 7.0 * std::exp(-326.0/T)) * 4.5e-4;
    
    real_t cold_cII_m = 6.2e4 * 1.38e-16 * (2.3e-8 * std::pow(T/10000.0, -0.5) * x + 1e-12) * std::exp(-6.2e4/T) * 3.5e-4 * 0.4;
    
    real_t cold_o_m;
    if (T <= 1e4) {
        cold_o_m = 2.3e4 * 1.38e-16 / 3.0 * (5.1e-9 * std::pow(T/10000.0, 0.57) * x + 1e-12) * std::exp(-2.3e4/T);
        cold_o_m += 4.9e4 * 1.38e-16 / 3.0 * (2.5e-9 * std::pow(T/10000.0, 0.57) * x + 1e-12) * std::exp(-4.9e4/T);
        cold_o_m += 2.6e4 * 1.38e-16 * (5.2e-9 * std::pow(T/10000.0, 0.57) * x + 1e-12) * std::exp(-2.6e4/T);
    } else {
        cold_o_m = 2.3e4 * 1.38e-16 / 3.0 * (5.1e-9 * std::pow(T/10000.0, 0.17) * x + 1e-12) * std::exp(-2.3e4/T);
        cold_o_m += 4.9e4 * 1.38e-16 / 3.0 * (2.5e-9 * std::pow(T/10000.0, 0.13) * x + 1e-12) * std::exp(-4.9e4/T);
        cold_o_m += 2.6e4 * 1.38e-16 * (5.2e-9 * std::pow(T/10000.0, 0.15) * x + 1e-12) * std::exp(-2.6e4/T);
    }
    cold_o_m *= 4.5e-4;
    
    real_t cold = (cold_cII + cold_o + cold_o_m + cold_cII_m);
    
    real_t G0 = 1.0 / 1.7;
    real_t param = G0 * std::sqrt(T) / (n * x);
    real_t epsilon = 4.9e-2 / (1.0 + std::pow(param/1925.0, 0.73)) + 3.7e-2 * std::pow(T/1e4, 0.7) / (1.0 + (param/5e3));
    real_t hot = 1e-24 * epsilon * G0;
    
    real_t bet = 0.74 / std::pow(T, 0.068);
    real_t cold_rec = 4.65e-30 * std::pow(T, 0.94) * std::pow(param, bet) * x;
    
    return hot * n - n * n * (cold + cold_rec);
}

void CoolingSolver::solve_cooling_ism(real_t& d, real_t& e_tot, real_t dt_sec, real_t dx, real_t mu) {
    // This is now integrated into apply_cooling for convenience
}

} // namespace ramses
