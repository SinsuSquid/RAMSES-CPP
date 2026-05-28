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
    if (ilevel > 0) {
        int igrid = grid_.get_headl(myid, ilevel);
        while (igrid > 0) { octs.push_back(igrid); igrid = grid_.next[igrid - 1]; }
    }
    
    // Level 0 has coarse cells.
    bool do_level_0 = (ilevel == 0);
    if (octs.empty() && !do_level_0) return;

    real_t gamma = grid_.gamma;
    real_t dtdx = dt / dx;
    int slope_type = config_.get_int("hydro_params", "slope_type", 1);
    bool verbose = config_.get_bool("run_params", "verbose", false);
    static bool first_print = true;
    if (first_print && ilevel == 0 && MpiManager::instance().rank() == 0 && verbose) {
        std::cout << "[HydroSolver] Using slope_type=" << slope_type << " Riemann=" << config_.get("hydro_params", "riemann", "hllc") << std::endl;
        first_print = false;
    }
    std::string riemann = config_.get("hydro_params", "riemann", "hllc");

    if (qm_level_.size() < (size_t)grid_.ncell * 3 * grid_.nvar) {
        qm_level_.assign((size_t)grid_.ncell * 3 * grid_.nvar, 0.0);
        qp_level_.assign((size_t)grid_.ncell * 3 * grid_.nvar, 0.0);
    }
    
    auto get_qr = [&](int idc_0, int idim, int iv) -> real_t& {
        return qm_level_[(idc_0 * 3 + idim) * grid_.nvar + (iv - 1)];
    };
    auto get_ql = [&](int idc_0, int idim, int iv) -> real_t& {
        return qp_level_[(idc_0 * 3 + idim) * grid_.nvar + (iv - 1)];
    };

    // 1. Trace step
    if (do_level_0) {
        for (int idc_0 = 0; idc_0 < grid_.ncoarse; ++idc_0) {
            real_t u_c[20], q_c[20];
            for(int iv=1; iv<=grid_.nvar; ++iv) u_c[iv-1] = grid_.uold(idc_0 + 1, iv);
            ctoprim(u_c, q_c, gamma);
            for (int idim = 0; idim < NDIM; ++idim) {
                real_t dq[20] = {0}, qm_tmp[20], qp_tmp[20];
                trace(q_c, dq, dt, dx, qm_tmp, qp_tmp, gamma);
                for(int iv=1; iv<=grid_.nvar; ++iv) {
                    get_qr(idc_0, idim, iv) = qm_tmp[iv-1];
                    get_ql(idc_0, idim, iv) = qp_tmp[iv-1];
                }
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
                compute_slopes(idc_0 + 1, icn, idim, dq, slope_type, dt, dx);
                if (idim > 0) { std::swap(q_rot[1], q_rot[1+idim]); std::swap(dq[1], dq[1+idim]); }
                real_t qm_tmp[20], qp_tmp[20];
                trace(q_rot, dq, dt, dx, qm_tmp, qp_tmp, gamma);
                if (idim > 0) { std::swap(qm_tmp[1], qm_tmp[1+idim]); std::swap(qp_tmp[1], qp_tmp[1+idim]); }
                for(int iv=1; iv<=grid_.nvar; ++iv) { 
                    get_qr(idc_0, idim, iv) = qm_tmp[iv-1]; 
                    get_ql(idc_0, idim, iv) = qp_tmp[iv-1]; 
                }
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
                    trace(q_nb, dq_null, dt, dx, qm_nb, qp_nb, gamma);
                    if (idim > 0) { std::swap(qm_nb[1], qm_nb[1+idim]); std::swap(qp_nb[1], qp_nb[1+idim]); }
                    
                    if (side == 0) { // Cell i, interface i-1/2. Neighbor is to the left.
                        for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = qm_nb[iv-1]; qr_f[iv-1] = get_ql(idc_0, idim, iv); }
                    } else { // Cell i, interface i+1/2. Neighbor is to the right.
                        for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qr(idc_0, idim, iv); qr_f[iv-1] = qp_nb[iv-1]; }
                    }
                } else {
                    int id_n0 = id_n - 1;
                    if (id_n > grid_.ncoarse && !grid_.son[id_n0] && grid_.get_cell_level(id_n) == ilevel) {
                        // Neighbor is at the same level; use pre-computed traces
                        if (side == 0) { // Cell i, interface i-1/2. Neighbor is i-1.
                            for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qr(id_n0, idim, iv); qr_f[iv-1] = get_ql(idc_0, idim, iv); }
                        } else { // Cell i, interface i+1/2. Neighbor is i+1.
                            for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qr(idc_0, idim, iv); qr_f[iv-1] = get_ql(id_n0, idim, iv); }
                        }
                    } else {
                        // Fallback for interface cells: use raw cell primitives
                        real_t u_nb[20], q_nb[20];
                        for(int iv=1; iv<=grid_.nvar; ++iv) u_nb[iv-1] = grid_.uold(id_n, iv);
                        ctoprim(u_nb, q_nb, gamma);
                        if (side == 0) { for(int iv=0; iv<grid_.nvar; ++iv) { ql_f[iv] = q_nb[iv]; qr_f[iv] = get_ql(idc_0, idim, iv+1); } }
                        else { for(int iv=0; iv<grid_.nvar; ++iv) { ql_f[iv] = get_qr(idc_0, idim, iv+1); qr_f[iv] = q_nb[iv]; } }
                    }
                }
                if (idim > 0) { std::swap(ql_f[1], ql_f[1+idim]); std::swap(qr_f[1], qr_f[1+idim]); }
                std::string riemann_type = config_.get("hydro_params", "riemann", "llf");
                if (riemann_type == "hllc") RiemannSolver::solve_hllc(ql_f, qr_f, flux, gamma);
                else if (riemann_type == "hll") RiemannSolver::solve_hll(ql_f, qr_f, flux, gamma);
                else RiemannSolver::solve_llf(ql_f, qr_f, flux, gamma);

                real_t mass_flux = flux[0];
                for (int iv = NDIM + 2; iv < grid_.nvar; ++iv) flux[iv] = mass_flux * ((mass_flux > 0) ? ql_f[iv] : qr_f[iv]);
                if (idim > 0) std::swap(flux[1], flux[1+idim]);
                real_t sign = (side == 0) ? 1.0 : -1.0;
                for(int iv=0; iv<grid_.nvar; ++iv) flux_sum[iv] += sign * flux[iv];
            }
        }
        for (int iv = 1; iv <= grid_.nvar; ++iv) {
            grid_.unew(idc_0 + 1, iv) = grid_.uold(idc_0 + 1, iv) + dtdx * flux_sum[iv-1];
            if (std::abs(grid_.unew(idc_0 + 1, iv)) > 1e10) {
                std::cout << "[DEBUG] Huge energy! idc=" << idc_0 + 1 << " val=" << grid_.unew(idc_0 + 1, iv)
                          << " uold=" << grid_.uold(idc_0 + 1, iv) << " flux=" << flux_sum[iv-1]
                          << " dtdx=" << dtdx << std::endl;
            }
        }
        real_t d_curr = std::clamp(grid_.unew(idc_0 + 1, 1), 1e-10, 1e6);
        grid_.unew(idc_0 + 1, 1) = d_curr;
        for (int i = 1; i <= NDIM; ++i) {
            grid_.unew(idc_0 + 1, 1 + i) = std::clamp(grid_.unew(idc_0 + 1, 1 + i), -d_curr * 10.0, d_curr * 10.0);
        }
        real_t v2_curr = 0.0; for(int i=1; i<=NDIM; ++i) { real_t v = grid_.unew(idc_0 + 1, 1+i)/d_curr; v2_curr += v*v; }
        int iener = NDIM + 2;
        real_t ei_curr = grid_.unew(idc_0 + 1, iener) - 0.5*d_curr*v2_curr;
        for(int ie=0; ie<nener_; ++ie) ei_curr -= grid_.unew(idc_0 + 1, iener+1+ie);
        if (ei_curr < d_curr*1e-10/(gamma-1.0)) {
            real_t e_non = 0; for(int ie=0; ie<nener_; ++ie) e_non += grid_.unew(idc_0 + 1, iener+1+ie);
            grid_.unew(idc_0 + 1, iener) = d_curr*1e-10/(gamma-1.0) + 0.5*d_curr*v2_curr + e_non;
        }
        grid_.unew(idc_0 + 1, iener) = std::min(grid_.unew(idc_0 + 1, iener), d_curr * 1e8);
    };

    if (do_level_0) {
        for (int idc_0 = 0; idc_0 < grid_.ncoarse; ++idc_0) {
            int icn[6];
            for (int idim = 0; idim < NDIM; ++idim) {
                for (int side = 0; side < 2; ++side) {
                    // This is a simplified boundary logic for coarse cells
                    int nx = grid_.nx, ny = grid_.ny, nz = grid_.nz;
                    int idx = idc_0;
                    int iz = idx / (nx * ny); int rem = idx % (nx * ny);
                    int iy = rem / nx; int ix = rem % nx;
                    int ixyz[3] = {ix, iy, iz};
                    int n[3] = {nx, ny, nz};
                    if (side == 0) {
                        if (ixyz[idim] == 0) icn[idim * 2 + side] = -1;
                        else {
                            int nxt_ixyz[3] = {ixyz[0], ixyz[1], ixyz[2]};
                            nxt_ixyz[idim]--;
                            icn[idim * 2 + side] = nxt_ixyz[2] * nx * ny + nxt_ixyz[1] * nx + nxt_ixyz[0] + 1;
                        }
                    } else {
                        if (ixyz[idim] == n[idim] - 1) icn[idim * 2 + side] = -2;
                        else {
                            int nxt_ixyz[3] = {ixyz[0], ixyz[1], ixyz[2]};
                            nxt_ixyz[idim]++;
                            icn[idim * 2 + side] = nxt_ixyz[2] * nx * ny + nxt_ixyz[1] * nx + nxt_ixyz[0] + 1;
                        }
                    }
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
    if (ilevel == 0) {
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
    if (ilevel == 0) {
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
    real_t d = std::max(u[0], 1e-10);
    d = std::min(d, 1e6);
    q[0] = d;
    real_t v2 = 0.0; for (int i = 1; i <= NDIM; ++i) {
        q[i] = std::clamp(u[i] / d, -10.0, 10.0);
        v2 += q[i] * q[i];
    }
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
        // Cap pressure to prevent sound speed explosion in vacuum
        q[iener] = std::min(q[iener], d * 1e6);
    }
    
    for (int iv = iener + 1; iv < grid_.nvar; ++iv) { if (iv < iener + 1 + nener_) q[iv] = u[iv] * (gamma - 1.0); else q[iv] = u[iv] / d; }
}

void HydroSolver::compute_slopes(int idc, const int icelln[6], int idim, real_t dq[20], int slope_type, real_t dt, real_t dx) {
    real_t ql[20], qc[20], qr[20], u_c[20];
    for(int iv=1; iv<=grid_.nvar; ++iv) u_c[iv-1]=grid_.uold(idc, iv);
    ctoprim(u_c, qc, grid_.gamma);
    auto get_nb_q = [&](int id_n, int side, real_t q_nb[20]) {
        if (id_n > 0) { real_t u_nb[20]; for(int iv=1; iv<=grid_.nvar; ++iv) u_nb[iv-1]=grid_.uold(id_n, iv); ctoprim(u_nb, q_nb, grid_.gamma); }
        else { for(int iv=0; iv<grid_.nvar; ++iv) q_nb[iv] = qc[iv]; int ib = -id_n; if(ib > 0 && ib <= (int)grid_.bound_type.size() && grid_.bound_type.at(ib-1) == 1) q_nb[1 + idim] *= -1.0; }
    };
    get_nb_q(icelln[idim*2], 0, ql); get_nb_q(icelln[idim*2+1], 1, qr);

    // slope_type 4 (Superbee) and 5 (Ultrabee): velocity-dependent limiters, 1D only
    // Type 5 applies ultrabee to density only; zero slope for all other vars
    if (slope_type == 4 || slope_type == 5) {
        real_t nu = qc[1 + idim] * dt / dx;  // Courant number: vel_component * dt/dx
        for (int iv = 0; iv < grid_.nvar; ++iv) {
            if (slope_type == 5 && iv != 0) { dq[iv] = 0.0; continue; }
            real_t dlft_raw = qc[iv] - ql[iv], drgt_raw = qr[iv] - qc[iv];
            real_t dlft_s, drgt_s;
            if (slope_type == 4) {
                // Superbee: symmetric scaling
                dlft_s = 2.0 / (1.0 + nu) * dlft_raw;
                drgt_s = 2.0 / (1.0 - nu) * drgt_raw;
            } else {
                // Ultrabee: one-sided scaling based on flow direction
                if (nu >= 0.0) {
                    dlft_s = 2.0 / (nu + 1e-10) * dlft_raw;
                    drgt_s = 2.0 / (1.0 - nu) * drgt_raw;
                } else {
                    dlft_s = 2.0 / (1.0 + nu) * dlft_raw;
                    drgt_s = 2.0 / (-nu + 1e-10) * drgt_raw;
                }
            }
            real_t dcen = 0.5 * (dlft_raw + drgt_raw);
            real_t dsgn = (dcen >= 0.0) ? 1.0 : -1.0;
            real_t dlim = std::min(std::abs(dlft_s), std::abs(drgt_s));
            if (dlft_s * drgt_s <= 0.0) dlim = 0.0;
            dq[iv] = dsgn * std::min(dlim, std::abs(dcen)) / dx;
        }
        return;
    }

    for (int iv = 0; iv < grid_.nvar; ++iv) {
        real_t dlft = qc[iv] - ql[iv], drgt = qr[iv] - qc[iv];
        if (dlft * drgt <= 0.0) {
            dq[iv] = 0.0;
        } else {
            real_t dcen = 0.5 * (dlft + drgt);
            real_t dsgn = (dcen >= 0.0) ? 1.0 : -1.0;
            real_t dlim = 0.0;
            if (slope_type == 1) {
                dlim = std::min(std::abs(dlft), std::abs(drgt));
            } else if (slope_type == 2 || slope_type == 3) {
                // MC (MonCen)
                dlim = std::min({2.0 * std::abs(dlft), 2.0 * std::abs(drgt), std::abs(dcen)});
            } else {
                dlim = std::min(std::abs(dlft), std::abs(drgt));
            }
            dq[iv] = dsgn * dlim / dx;
        }
    }
}

void HydroSolver::trace(const real_t q[], const real_t dq[], real_t dt, real_t dx, real_t qm[], real_t qp[], real_t gamma) {
    real_t r = std::max(q[0], 1e-10), u = q[1], p = std::max(q[NDIM+1], r * 1e-20);
    real_t sr0 = -u * dq[0] - dq[1] * r;
    real_t su0 = -u * dq[1] - dq[NDIM+1] / r;
    int iener = NDIM + 1;
    for(int ie=0; ie<nener_; ++ie) su0 -= dq[iener+1+ie] / r;
    
    real_t cs2 = RiemannSolver::get_cs2(r, p, gamma);
    real_t sp0 = -u * dq[NDIM+1] - dq[1] * r * cs2;

    // Stabilize source terms in vacuum regions to prevent timestep collapse
    su0 = std::max(-r * 1e3, std::min(r * 1e3, su0));
    sp0 = std::max(-r * 1e3, std::min(r * 1e3, sp0));

    for (int iv = 0; iv < grid_.nvar; ++iv) {
        real_t dqi = dq[iv];
        real_t dq_dt = 0;
        if (iv == 0) dq_dt = sr0;
        else if (iv == 1) dq_dt = su0;
        else if (iv == NDIM + 1) dq_dt = sp0;

        qp[iv] = q[iv] - 0.5 * dx * dqi + 0.5 * dt * dq_dt;
        qm[iv] = q[iv] + 0.5 * dx * dqi + 0.5 * dt * dq_dt;
        if (iv == 0) {
            if (qp[iv] < 1e-10) qp[iv] = r;
            if (qp[iv] > 1e6) qp[iv] = 1e6;
            if (qm[iv] < 1e-10) qm[iv] = r;
            if (qm[iv] > 1e6) qm[iv] = 1e6;
        } else if (iv == iener) {
            if (qp[iv] < r * 1e-20) qp[iv] = r * 1e-20;
            if (qm[iv] < r * 1e-20) qm[iv] = r * 1e-20;
            qp[iv] = std::min(qp[iv], r * 1e2);
            qm[iv] = std::min(qm[iv], r * 1e2);
        } else if (iv > 0 && iv <= NDIM) {
            qp[iv] = std::clamp(qp[iv], -10.0, 10.0);
            qm[iv] = std::clamp(qm[iv], -10.0, 10.0);
        }
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
    std::cout << "[DEBUG] compute_courant_step called for level " << ilevel << std::endl;
    int myid = MpiManager::instance().rank() + 1;
    real_t dt_max = 1e30;
    real_t max_dtdx = 0.1;
    real_t max_dt = max_dtdx * dx;

    auto get_cs = [&](real_t d, real_t p) {
        return std::sqrt(RiemannSolver::get_cs2(d, p, gamma));
    };

    if (ilevel == 0) {
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
            p = std::min(p, d * 1e2);
            }
            
            real_t cs = get_cs(d, p);
            real_t v_mag = std::sqrt(v2);
            if (id == 1 && ilevel == 0) {
                std::cout << "[DEBUG] Courant check: d=" << d << " p=" << p << " cs=" << cs << " v_mag=" << v_mag << " dt=" << courant_factor * dx / (v_mag + cs) << std::endl;
            }
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
            p = std::min(p, d * 1e2);
                }
                
                real_t cs = get_cs(d, p);
                real_t v_mag = std::sqrt(v2);
                static bool first_print = true;
                if (first_print) {
                    std::cout << "[DEBUG] Courant check (L" << ilevel << "): d=" << d << " p=" << p << " cs=" << cs << " v_mag=" << v_mag << " dt=" << courant_factor * dx / (v_mag + cs) << std::endl;
                    first_print = false;
                }
                real_t dt_cfl = courant_factor * dx / (v_mag + cs);
                real_t dt_grav = 1e30;
                real_t max_acc = 0;
                for(int idim=1; idim<=NDIM; ++idim) {
                    real_t acc = std::abs(grid_.f(id, idim));
                    if (acc > max_acc) max_acc = acc;
                    if (acc > 0) dt_grav = std::min(dt_grav, courant_factor * std::sqrt(dx / acc));
                }
                real_t dt_cand = std::min(dt_cfl, dt_grav);
                if (dt_cand < 1e-8) {
                    std::cout << "[DEBUG] Timestep collapse! idc=" << id << " max_acc=" << max_acc << " dt_cfl=" << dt_cfl << " dt_grav=" << dt_grav << " v_mag=" << v_mag << " cs=" << cs << std::endl;
                }
                if (id == 1 && ilevel == 0) {
                    // we can check flux later
                }
                dt_max = std::min(dt_max, dt_cand);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    real_t final_dt = std::min(dt_max, max_dt);
    if (MpiManager::instance().rank() == 0 && ilevel == 0) {
        std::cout << "[DEBUG] compute_courant_step L0 returning dt=" << final_dt << " (max_dt=" << max_dt << ")" << std::endl;
    }
    return final_dt;
}

void HydroSolver::get_diagnostics(int ilevel, real_t dx, real_t& min_d, real_t& max_v, real_t& min_t, real_t& max_t) {}

void HydroSolver::add_gravity_source_terms(int ilevel, real_t dt) {
    int myid = MpiManager::instance().rank() + 1, n2d_val = (1 << NDIM);
    int iener = NDIM + 2;
    if (ilevel == 1) {
        for (int idc = 1; idc <= grid_.ncoarse; ++idc) {
            if (grid_.son.at(idc - 1) > 0) continue;
            real_t d = std::clamp(grid_.uold(idc, 1), 1e-10, 1e6);
            for (int idim = 1; idim <= NDIM; ++idim) {
                real_t f = grid_.f(idc, idim);
                grid_.unew(idc, 1 + idim) += d * f * dt;
                grid_.unew(idc, 1 + idim) = std::clamp(grid_.unew(idc, 1 + idim), -d * 10.0, d * 10.0);
                grid_.unew(idc, iener) += grid_.uold(idc, 1 + idim) * f * dt + 0.5 * d * f * f * dt * dt;
            }
        }
    }
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= n2d_val; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son.at(idc - 1) > 0) continue; 
            real_t d = std::clamp(grid_.uold(idc, 1), 1e-10, 1e6);
            for (int idim = 1; idim <= NDIM; ++idim) {
                real_t f = grid_.f(idc, idim);
                grid_.unew(idc, 1 + idim) += d * f * dt;
                grid_.unew(idc, 1 + idim) = std::clamp(grid_.unew(idc, 1 + idim), -d * 10.0, d * 10.0);
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
        real_t d = std::clamp(grid_.uold(idc, 1), 1e-10, 1e6);
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

        for (int idim = 1; idim <= NDIM; ++idim) {
            grid_.uold(idc, 1 + idim) += d * grid_.f(idc, idim) * dt;
            grid_.uold(idc, 1 + idim) = std::clamp(grid_.uold(idc, 1 + idim), -d * 10.0, d * 10.0);
        }
        
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
