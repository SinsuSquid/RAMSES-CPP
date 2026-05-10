#include "ramses/HydroSolver.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/RiemannSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/SolverFactory.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

namespace ramses {

HydroSolver::~HydroSolver() {}

void HydroSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    int myid = MpiManager::instance().rank() + 1;
    std::vector<int> octs;
    if (ilevel > 1) {
        int igrid = grid_.get_headl(myid, ilevel);
        while (igrid > 0) { octs.push_back(igrid); igrid = grid_.next[igrid - 1]; }
    }
    
    // Level 1 doesn't have grids, it has coarse cells.
    bool do_level_1 = (ilevel == 1);
    if (octs.empty() && !do_level_1) return;

    real_t gamma = grid_.gamma;
    real_t dtdx = dt / dx;
    int slope_type = config_.get_int("hydro_params", "slope_type", 1);
    std::string riemann = config_.get("hydro_params", "riemann", "hllc");

    if (qm_level_.size() < (size_t)grid_.ncell * 3 * grid_.nvar) {
        qm_level_.assign((size_t)grid_.ncell * 3 * grid_.nvar, 0.0);
        qp_level_.assign((size_t)grid_.ncell * 3 * grid_.nvar, 0.0);
    }
    
    auto get_q = [&](int idc_0, int idim, int iv) -> real_t& {
        return qm_level_[(idc_0 * 3 + idim) * grid_.nvar + (iv - 1)];
    };
    auto get_qp = [&](int idc_0, int idim, int iv) -> real_t& {
        return qp_level_[(idc_0 * 3 + idim) * grid_.nvar + (iv - 1)];
    };

    // 1. Trace step
    if (do_level_1) {
        for (int idc_0 = 0; idc_0 < grid_.ncoarse; ++idc_0) {
            real_t u_c[20], q_c[20];
            for(int iv=1; iv<=grid_.nvar; ++iv) u_c[iv-1] = grid_.uold(idc_0 + 1, iv);
            ctoprim(u_c, q_c, gamma);
            for (int idim = 0; idim < NDIM; ++idim) {
                real_t dq[20] = {0}, qm_tmp[20], qp_tmp[20];
                trace(q_c, dq, dtdx, qm_tmp, qp_tmp, gamma);
                for(int iv=1; iv<=grid_.nvar; ++iv) { get_q(idc_0, idim, iv) = qm_tmp[iv-1]; get_qp(idc_0, idim, iv) = qp_tmp[iv-1]; }
            }
        }
    }

    for (int ig : octs) {
        int ign[7]; grid_.get_nbor_grids(ig, ign);
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc_0 = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
            int icn[6]; grid_.get_nbor_cells(ign, ic, icn, ig);
            real_t u_c[20], q_c[20];
            for(int iv=1; iv<=grid_.nvar; ++iv) u_c[iv-1] = grid_.uold(idc_0 + 1, iv);
            ctoprim(u_c, q_c, gamma);
            for (int idim = 0; idim < NDIM; ++idim) {
                real_t dq[20], q_rot[20];
                for(int iv=0; iv<grid_.nvar; ++iv) q_rot[iv] = q_c[iv];
                compute_slopes(idc_0 + 1, icn, idim, dq, slope_type);
                if (idim > 0) { std::swap(q_rot[1], q_rot[1+idim]); std::swap(dq[1], dq[1+idim]); }
                real_t qm_tmp[20], qp_tmp[20];
                trace(q_rot, dq, dtdx, qm_tmp, qp_tmp, gamma);
                if (idim > 0) { std::swap(qm_tmp[1], qm_tmp[1+idim]); std::swap(qp_tmp[1], qp_tmp[1+idim]); }
                for(int iv=1; iv<=grid_.nvar; ++iv) { get_q(idc_0, idim, iv) = qm_tmp[iv-1]; get_qp(idc_0, idim, iv) = qp_tmp[iv-1]; }
            }
        }
    }

    // 2. Flux step helper
    auto compute_fluxes = [&](int idc_0, const int icn[6], int ig_info) {
        if (grid_.son.at(idc_0) != 0) return;
        real_t flux_sum[20] = {0};
        for (int idim = 0; idim < NDIM; ++idim) {
            for (int side = 0; side < 2; ++side) {
                int id_n = icn[idim * 2 + side];
                real_t ql_f[20], qr_f[20], flux[20];
                if (id_n <= 0) {
                    int ibound = -id_n; int btype = 1;
                    if (ibound > 0 && ibound <= (int)grid_.bound_type.size()) btype = grid_.bound_type.at(ibound - 1);
                    real_t u_nb[20], q_nb[20], qm_nb[20], qp_nb[20], dq_null[20] = {0};
                    for(int iv=1; iv<=grid_.nvar; ++iv) u_nb[iv-1] = grid_.uold(idc_0 + 1, iv);
                    if (btype == 1) u_nb[1 + idim] *= -1.0;
                    ctoprim(u_nb, q_nb, gamma);
                    if (idim > 0) std::swap(q_nb[1], q_nb[1+idim]);
                    trace(q_nb, dq_null, dtdx, qm_nb, qp_nb, gamma);
                    if (idim > 0) { std::swap(qm_nb[1], qm_nb[1+idim]); std::swap(qp_nb[1], qp_nb[1+idim]); }
                    if (side == 0) { for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = qp_nb[iv-1]; qr_f[iv-1] = get_q(idc_0, idim, iv); } }
                    else { for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qp(idc_0, idim, iv); qr_f[iv-1] = qm_nb[iv-1]; } }
                } else {
                    real_t u_nb[20], q_nb[20], dq_nb[20], qm_nb[20], qp_nb[20];
                    for(int iv=1; iv<=grid_.nvar; ++iv) u_nb[iv-1] = grid_.uold(id_n, iv);
                    ctoprim(u_nb, q_nb, gamma);
                    int icn_nb[6], ign_nb[7] = {0}, ic_nb = 1, ig_nb = 0;
                    if (id_n > grid_.ncoarse) {
                        ig_nb = ((id_n - grid_.ncoarse - 1) % grid_.ngridmax) + 1;
                        ic_nb = ((id_n - grid_.ncoarse - 1) / grid_.ngridmax) + 1;
                        grid_.get_nbor_grids(ig_nb, ign_nb);
                    }
                    grid_.get_nbor_cells(ign_nb, ic_nb, icn_nb, ig_nb);
                    compute_slopes(id_n, icn_nb, idim, dq_nb, slope_type);
                    if (idim > 0) { std::swap(q_nb[1], q_nb[1+idim]); std::swap(dq_nb[1], dq_nb[1+idim]); }
                    trace(q_nb, dq_nb, dtdx, qm_nb, qp_nb, gamma);
                    if (idim > 0) { std::swap(qm_nb[1], qm_nb[1+idim]); std::swap(qp_nb[1], qp_nb[1+idim]); }
                    if (side == 0) { for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = qp_nb[iv-1]; qr_f[iv-1] = get_q(idc_0, idim, iv); } }
                    else { for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qp(idc_0, idim, iv); qr_f[iv-1] = qm_nb[iv-1]; } }
                }
                if (idim > 0) { std::swap(ql_f[1], ql_f[1+idim]); std::swap(qr_f[1], qr_f[1+idim]); }
                if (riemann == "hllc") RiemannSolver::solve_hllc(ql_f, qr_f, flux, gamma);
                else if (riemann == "llf") RiemannSolver::solve_llf(ql_f, qr_f, flux, gamma);
                else RiemannSolver::solve_hll(ql_f, qr_f, flux, gamma);
                real_t mass_flux = flux[0];
                for (int iv = NDIM + 2; iv < grid_.nvar; ++iv) flux[iv] = mass_flux * ((mass_flux > 0) ? ql_f[iv] : qr_f[iv]);
                if (idim > 0) std::swap(flux[1], flux[1+idim]);
                real_t sign = (side == 0) ? 1.0 : -1.0;
                for(int iv=0; iv<grid_.nvar; ++iv) flux_sum[iv] += sign * flux[iv];
            }
        }
        for (int iv = 1; iv <= nvar_hydro_; ++iv) grid_.unew(idc_0 + 1, iv) = grid_.uold(idc_0 + 1, iv) + dtdx * flux_sum[iv-1];
        real_t d_curr = std::max(grid_.unew(idc_0 + 1, 1), 1e-10); grid_.unew(idc_0 + 1, 1) = d_curr;
        real_t v2_curr = 0.0; for(int i=1; i<=NDIM; ++i) { real_t v = grid_.unew(idc_0 + 1, 1+i)/d_curr; v2_curr += v*v; }
        int iener = NDIM + 2;
        real_t ei_curr = grid_.unew(idc_0 + 1, iener) - 0.5*d_curr*v2_curr;
        for(int ie=0; ie<nener_; ++ie) ei_curr -= grid_.unew(idc_0 + 1, iener+1+ie);
        if (ei_curr < d_curr*1e-10/(gamma-1.0)) {
            real_t e_non = 0; for(int ie=0; ie<nener_; ++ie) e_non += grid_.unew(idc_0 + 1, iener+1+ie);
            grid_.unew(idc_0 + 1, iener) = d_curr*1e-10/(gamma-1.0) + 0.5*d_curr*v2_curr + e_non;
        }
    };

    if (do_level_1) {
        for (int idc_0 = 0; idc_0 < grid_.ncoarse; ++idc_0) {
            int icn[6];
            for (int idim = 0; idim < 3; ++idim) {
                for (int side = 0; side < 2; ++side) {
                    // This is a simplified boundary logic for coarse cells
                    int nx = grid_.nx, ny = grid_.ny, nz = grid_.nz;
                    int idx = idc_0;
                    int iz = idx / (nx * ny); int rem = idx % (nx * ny);
                    int iy = rem / nx; int ix = rem % nx;
                    int ixyz[3] = {ix, iy, iz};
                    int n[3] = {nx, ny, nz};
                    if (side == 0) ixyz[idim] = (ixyz[idim] - 1 + n[idim]) % n[idim];
                    else ixyz[idim] = (ixyz[idim] + 1) % n[idim];
                    icn[idim * 2 + side] = ixyz[2] * nx * ny + ixyz[1] * nx + ixyz[0] + 1;
                }
            }
            compute_fluxes(idc_0, icn, 0);
        }
    }

    for (int ig : octs) {
        int ign[7]; grid_.get_nbor_grids(ig, ign);
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc_0 = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
            int icn[6]; grid_.get_nbor_cells(ign, ic, icn, ig);
            compute_fluxes(idc_0, icn, ig);
        }
    }
}

void HydroSolver::set_unew(int ilevel) {
    int myid = MpiManager::instance().rank() + 1, n2d_val = (1 << NDIM);
    if (ilevel == 1) {
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.unew(i, iv) = grid_.uold(i, iv);
        }
    }
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= n2d_val; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.unew(idc, iv) = grid_.uold(idc, iv);
        }
        igrid = grid_.next[igrid - 1];
    }
}

void HydroSolver::set_uold(int ilevel) {
    int myid = MpiManager::instance().rank() + 1, n2d_val = (1 << NDIM);
    if (ilevel == 1) {
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            if (grid_.son.at(i - 1) == 0) {
                for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.uold(i, iv) = grid_.unew(i, iv);
            }
        }
    }
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= n2d_val; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son.at(idc - 1) == 0) {
                for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.uold(idc, iv) = grid_.unew(idc, iv);
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

void HydroSolver::ctoprim(const real_t u[], real_t q[], real_t gamma) {
    real_t d = std::max(u[0], 1e-10); q[0] = d;
    real_t v2 = 0.0; for (int i = 1; i <= NDIM; ++i) { q[i] = u[i] / d; v2 += q[i] * q[i]; }
    int iener = NDIM + 1; // 0-based index in local array
    
    if (params::barotropic_eos) {
        // Calculate pressure from density only
        real_t nH = d * params::units_density / (params::mu_gas * constants::mH);
        real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
        real_t temp_mu = T2;
        if (params::barotropic_eos_form == "polytrope") {
            temp_mu = T2 * std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
        } else if (params::barotropic_eos_form == "double_polytrope") {
            temp_mu = T2 * (1.0 + std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0));
        }
        q[iener] = d * temp_mu / (params::units_velocity * params::units_velocity);
    } else {
        real_t e_thermal_dens = u[iener] - 0.5 * d * v2;
        for (int ie = 0; ie < nener_; ++ie) e_thermal_dens -= u[iener + 1 + ie];
        q[iener] = std::max(e_thermal_dens * (gamma - 1.0), d * 1e-10);
    }
    
    for (int iv = iener + 1; iv < grid_.nvar; ++iv) { if (iv < iener + 1 + nener_) q[iv] = u[iv] * (gamma - 1.0); else q[iv] = u[iv] / d; }
}

void HydroSolver::compute_slopes(int idc, const int icelln[6], int idim, real_t dq[20], int slope_type) {
    real_t ql[20], qc[20], qr[20], u_c[20];
    for(int iv=1; iv<=grid_.nvar; ++iv) u_c[iv-1]=grid_.uold(idc, iv);
    ctoprim(u_c, qc, grid_.gamma);
    auto get_nb_q = [&](int id_n, int side, real_t q_nb[20]) {
        if (id_n > 0) { real_t u_nb[20]; for(int iv=1; iv<=grid_.nvar; ++iv) u_nb[iv-1]=grid_.uold(id_n, iv); ctoprim(u_nb, q_nb, grid_.gamma); }
        else { for(int iv=0; iv<grid_.nvar; ++iv) q_nb[iv] = qc[iv]; int ib = -id_n; if(ib > 0 && ib <= (int)grid_.bound_type.size() && grid_.bound_type.at(ib-1) == 1) q_nb[1 + idim] *= -1.0; }
    };
    get_nb_q(icelln[idim*2], 0, ql); get_nb_q(icelln[idim*2+1], 1, qr);
    for (int iv = 0; iv < grid_.nvar; ++iv) {
        real_t dlft = qc[iv] - ql[iv], drgt = qr[iv] - qc[iv];
        if (dlft * drgt <= 0.0) dq[iv] = 0.0;
        else { real_t sgn = (dlft >= 0.0) ? 1.0 : -1.0; dq[iv] = sgn * std::min(std::abs(dlft), std::abs(drgt)); }
    }
}

void HydroSolver::trace(const real_t q[], const real_t dq[], real_t dtdx, real_t qm[], real_t qp[], real_t gamma) {
    real_t r = std::max(q[0], 1e-10), u = q[1], p = std::max(q[NDIM+1], 1e-10);
    real_t sr0 = -u * dq[0] - dq[1] * r;
    real_t su0 = -u * dq[1] - dq[NDIM+1] / r;
    int iener = NDIM + 1;
    for(int ie=0; ie<nener_; ++ie) su0 -= dq[iener+1+ie] / r;
    
    real_t cs2 = gamma * p / r;
    if (params::barotropic_eos) {
        real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
        real_t v2_unit = params::units_velocity * params::units_velocity;
        if (params::barotropic_eos_form == "isothermal") cs2 = T2 / v2_unit;
        else if (params::barotropic_eos_form == "polytrope") cs2 = params::polytrope_index * T2 * std::pow(r * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0) / v2_unit;
        else if (params::barotropic_eos_form == "double_polytrope") {
            real_t fac = std::pow(r * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
            cs2 = T2 * (1.0 + params::polytrope_index * fac) / v2_unit;
        }
    }
    
    real_t sp0 = -u * dq[iener] - dq[1] * r * cs2;
    for (int iv = 0; iv < grid_.nvar; ++iv) {
        real_t dqi = dq[iv], src = -u * dqi;
        if (iv == 0) src = sr0; else if (iv == 1) src = su0; else if (iv == iener) src = sp0;
        else if (iv > iener && iv <= iener + nener_) src = -u * dqi - dq[1] * gamma * q[iv];
        qp[iv] = q[iv] - 0.5 * dqi + 0.5 * dtdx * src;
        qm[iv] = q[iv] + 0.5 * dqi + 0.5 * dtdx * src;
    }
}

void HydroSolver::interpol_hydro(const real_t u1[7][64], real_t u2[8][64]) {
    int nvar = grid_.nvar; real_t gam = grid_.gamma;
    real_t q1[7][64], q2[8][64];
    for (int i = 0; i < 7; ++i) ctoprim(u1[i], q1[i], gam);
    for (int iv = 0; iv < nvar; ++iv) {
        real_t slopes[3] = {0,0,0};
        for (int idim = 0; idim < NDIM; ++idim) {
            real_t dlft = q1[0][iv] - q1[2*idim+1][iv], drgt = q1[2*idim+2][iv] - q1[0][iv];
            if (dlft * drgt > 0.0) { real_t sgn = (dlft >= 0.0) ? 1.0 : -1.0; slopes[idim] = sgn * std::min({2.0*std::abs(dlft), 2.0*std::abs(drgt), 0.5*std::abs(dlft+drgt)}); }
        }
        for (int i = 0; i < 8; ++i) {
            int ix = i & 1, iy = (i & 2) >> 1, iz = (i & 4) >> 2;
            q2[i][iv] = q1[0][iv] + (ix-0.5)*slopes[0] + (iy-0.5)*slopes[1] + (iz-0.5)*slopes[2];
        }
    }
    int iener = NDIM + 1;
    for (int i = 0; i < 8; ++i) {
        real_t d = std::max(q2[i][0], 1e-10);
        real_t p = 0;
        if (params::barotropic_eos) {
            real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
            real_t temp_mu = T2;
            if (params::barotropic_eos_form == "polytrope") temp_mu = T2 * std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
            else if (params::barotropic_eos_form == "double_polytrope") temp_mu = T2 * (1.0 + std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0));
            p = d * temp_mu / (params::units_velocity * params::units_velocity);
        } else {
            p = std::max(q2[i][iener], 1e-20);
        }
        
        u2[i][0] = d;
        real_t v2 = 0;
        for(int d_idx=1; d_idx<=NDIM; ++d_idx) { u2[i][d_idx] = d * q2[i][d_idx]; v2 += q2[i][d_idx] * q2[i][d_idx]; }
        real_t e_thermal = p / (gam - 1.0), e_kinetic = 0.5 * d * v2, e_nonthermal = 0.0;
        for (int ie = 0; ie < nener_; ++ie) { real_t e_rad = q2[i][iener+1+ie] / (gam - 1.0); u2[i][iener+1+ie] = e_rad; e_nonthermal += e_rad; }
        u2[i][iener] = e_thermal + e_kinetic + e_nonthermal;
        for (int iv = iener + 1 + nener_; iv < nvar; ++iv) u2[i][iv] = d * q2[i][iv];
    }
}

real_t HydroSolver::compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor) {
    int myid = MpiManager::instance().rank() + 1;
    real_t dt_max = 1e30;
    
    auto get_cs = [&](real_t d, real_t p) {
        if (params::barotropic_eos) {
            // Sound speed for polytrope: cs^2 = dP/drho
            real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
            real_t v2_unit = params::units_velocity * params::units_velocity;
            if (params::barotropic_eos_form == "isothermal") return std::sqrt(T2 / v2_unit);
            if (params::barotropic_eos_form == "polytrope") {
                return std::sqrt(params::polytrope_index * T2 * std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0) / v2_unit);
            }
            if (params::barotropic_eos_form == "double_polytrope") {
                real_t fac = std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
                real_t cs2 = T2 * (1.0 + params::polytrope_index * fac) / v2_unit;
                return std::sqrt(cs2);
            }
        }
        return std::sqrt(gamma * p / d);
    };

    if (ilevel == 1) {
        for (int id = 1; id <= grid_.ncoarse; ++id) {
            if (grid_.son.at(id - 1) > 0) continue;
            real_t d = std::max(grid_.uold(id, 1), 1e-10), v2 = 0.0;
            for (int i = 1; i <= NDIM; ++i) { real_t v = grid_.uold(id, 1 + i) / d; v2 += v * v; }
            int iener = NDIM + 2;
            
            real_t p = 0;
            if (params::barotropic_eos) {
                real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
                real_t temp_mu = T2;
                if (params::barotropic_eos_form == "polytrope") temp_mu = T2 * std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
                else if (params::barotropic_eos_form == "double_polytrope") temp_mu = T2 * (1.0 + std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0));
                p = d * temp_mu / (params::units_velocity * params::units_velocity);
            } else {
                p = std::max((grid_.uold(id, iener) - 0.5 * d * v2) * (gamma - 1.0), d * 1e-10);
            }
            
            real_t cs = get_cs(d, p);
            real_t v_mag = std::sqrt(v2);
            dt_max = std::min(dt_max, courant_factor * dx / (v_mag + cs));
        }
    } else {
        int igrid = grid_.get_headl(myid, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int id = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                if (grid_.son.at(id - 1) > 0) continue;
                real_t d = std::max(grid_.uold(id, 1), 1e-10), v2 = 0.0;
                for (int i = 1; i <= NDIM; ++i) { real_t v = grid_.uold(id, 1 + i) / d; v2 += v * v; }
                int iener = NDIM + 2;
                
                real_t p = 0;
                if (params::barotropic_eos) {
                    real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
                    real_t temp_mu = T2;
                    if (params::barotropic_eos_form == "polytrope") temp_mu = T2 * std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
                    else if (params::barotropic_eos_form == "double_polytrope") temp_mu = T2 * (1.0 + std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0));
                    p = d * temp_mu / (params::units_velocity * params::units_velocity);
                } else {
                    p = std::max((grid_.uold(id, iener) - 0.5 * d * v2) * (gamma - 1.0), d * 1e-10);
                }
                
                real_t cs = get_cs(d, p);
                real_t v_mag = std::sqrt(v2);
                dt_max = std::min(dt_max, courant_factor * dx / (v_mag + cs));
                for(int idim=1; idim<=NDIM; ++idim) {
                    real_t acc = std::abs(grid_.f(id, idim));
                    if (acc > 0) dt_max = std::min(dt_max, courant_factor * std::sqrt(dx / acc));
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    return dt_max;
}

void HydroSolver::get_diagnostics(int ilevel, real_t dx, real_t& min_d, real_t& max_v, real_t& min_t, real_t& max_t) {}

void HydroSolver::add_gravity_source_terms(int ilevel, real_t dt) {
    int myid = MpiManager::instance().rank() + 1, n2d_val = (1 << NDIM);
    int iener = NDIM + 2;
    if (ilevel == 1) {
        for (int idc = 1; idc <= grid_.ncoarse; ++idc) {
            if (grid_.son.at(idc - 1) > 0) continue;
            real_t d = std::max(grid_.uold(idc, 1), 1e-10);
            for (int idim = 1; idim <= NDIM; ++idim) {
                real_t f = grid_.f(idc, idim);
                grid_.unew(idc, 1 + idim) += d * f * dt;
                grid_.unew(idc, iener) += grid_.uold(idc, 1 + idim) * f * dt + 0.5 * d * f * f * dt * dt;
            }
        }
    }
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= n2d_val; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son.at(idc - 1) > 0) continue; 
            real_t d = std::max(grid_.uold(idc, 1), 1e-10);
            for (int idim = 1; idim <= NDIM; ++idim) {
                real_t f = grid_.f(idc, idim);
                grid_.unew(idc, 1 + idim) += d * f * dt;
                grid_.unew(idc, iener) += grid_.uold(idc, 1 + idim) * f * dt + 0.5 * d * f * f * dt * dt;
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

void HydroSolver::synchro_hydro_fine(int ilevel, real_t dt) {
    int myid = MpiManager::instance().rank() + 1, n2d_val = (1 << NDIM);
    int iener = NDIM + 2;
    real_t gam = grid_.gamma;

    auto apply_synchro = [&](int idc) {
        if (grid_.son.at(idc - 1) > 0) return;
        real_t d = std::max(grid_.uold(idc, 1), 1e-10);
        real_t v2_old = 0;
        for(int idim=1; idim<=NDIM; ++idim) { real_t v = grid_.uold(idc, 1+idim)/d; v2_old += v*v; }
        
        real_t e_int = 0;
        if (params::barotropic_eos) {
            real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
            real_t temp_mu = T2;
            if (params::barotropic_eos_form == "polytrope") temp_mu = T2 * std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
            else if (params::barotropic_eos_form == "double_polytrope") temp_mu = T2 * (1.0 + std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0));
            real_t p = d * temp_mu / (params::units_velocity * params::units_velocity);
            e_int = p / (gam - 1.0);
        } else {
            e_int = std::max(grid_.uold(idc, iener) - 0.5*d*v2_old, d * 1e-10 / (gam - 1.0));
        }

        for (int idim = 1; idim <= NDIM; ++idim) grid_.uold(idc, 1 + idim) += d * grid_.f(idc, idim) * dt;
        
        real_t v2_new = 0;
        for(int idim=1; idim<=NDIM; ++idim) { real_t v = grid_.uold(idc, 1+idim)/d; v2_new += v*v; }
        grid_.uold(idc, iener) = e_int + 0.5*d*v2_new;
    };

    if (ilevel == 1) {
        for (int idc = 1; idc <= grid_.ncoarse; ++idc) apply_synchro(idc);
    } else {
        int igrid = grid_.get_headl(myid, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= n2d_val; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                apply_synchro(idc);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

} // namespace ramses
