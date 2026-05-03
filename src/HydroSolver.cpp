#include "ramses/HydroSolver.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/RiemannSolver.hpp"
#include "ramses/Parameters.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace ramses {

void HydroSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    int myid = MpiManager::instance().rank() + 1;
    std::vector<int> octs;
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) { octs.push_back(igrid); igrid = grid_.next[igrid - 1]; }
    if (octs.empty()) return;

    real_t gamma = grid_.gamma;
    real_t dtdx = dt / dx;
    int slope_type = config_.get_int("hydro_params", "slope_type", 1);
    std::string riemann = config_.get("hydro_params", "riemann", "hllc");

    if (qm_level_.size() < (size_t)grid_.ncell * 3 * grid_.nvar) {
        qm_level_.assign(grid_.ncell * 3 * grid_.nvar, 0.0);
        qp_level_.assign(grid_.ncell * 3 * grid_.nvar, 0.0);
    }
    
    auto get_q = [&](int icell, int idim, int iv) -> real_t& {
        return qm_level_[((icell - 1) * 3 + idim) * grid_.nvar + (iv - 1)];
    };
    auto get_qp = [&](int icell, int idim, int iv) -> real_t& {
        return qp_level_[((icell - 1) * 3 + idim) * grid_.nvar + (iv - 1)];
    };

    for (int ig : octs) {
        int ign[7]; grid_.get_nbor_grids(ig, ign);
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
            int icn[6]; grid_.get_nbor_cells(ign, ic, icn, ig);
            real_t u_c[20], q_c[20];
            for(int iv=1; iv<=grid_.nvar; ++iv) u_c[iv-1] = grid_.uold(idc, iv);
            ctoprim(u_c, q_c, gamma);
            for (int idim = 0; idim < NDIM; ++idim) {
                real_t dq[20], q_rot[20];
                for(int iv=0; iv<grid_.nvar; ++iv) q_rot[iv] = q_c[iv];
                compute_slopes(idc, icn, idim, dq, slope_type);
                if (idim > 0) { std::swap(q_rot[1], q_rot[1+idim]); std::swap(dq[1], dq[1+idim]); }
                real_t qm_tmp[20], qp_tmp[20];
                trace(q_rot, dq, dtdx, qm_tmp, qp_tmp, gamma);
                if (idim > 0) { std::swap(qm_tmp[1], qm_tmp[1+idim]); std::swap(qp_tmp[1], qp_tmp[1+idim]); }
                for(int iv=1; iv<=grid_.nvar; ++iv) { get_q(idc, idim, iv) = qm_tmp[iv-1]; get_qp(idc, idim, iv) = qp_tmp[iv-1]; }
            }
        }
    }

    for (int ig : octs) {
        int ign[7]; grid_.get_nbor_grids(ig, ign);
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
            if (grid_.son[idc] != 0) continue;
            int icn[6]; grid_.get_nbor_cells(ign, ic, icn, ig);
            real_t flux_sum[20] = {0};
            for (int idim = 0; idim < NDIM; ++idim) {
                for (int side = 0; side < 2; ++side) {
                    int id_n = icn[idim * 2 + side];
                    real_t ql_f[20], qr_f[20], flux[20];
                    
                    if (id_n <= 0) {
                        int ibound = (idim == 0) ? ((side == 0) ? 1 : 2) : (idim == 1) ? ((side == 0) ? 3 : 4) : ((side == 0) ? 5 : 6);
                        int btype = 1;
                        if (grid_.nboundary > 0 && ibound <= (int)grid_.bound_type.size()) btype = grid_.bound_type[ibound - 1];
                        real_t u_nb[20], q_nb[20], qm_nb[20], qp_nb[20];
                        for(int iv=1; iv<=grid_.nvar; ++iv) u_nb[iv-1] = grid_.uold(idc, iv);
                        if (btype == 1) u_nb[1 + idim] *= -1.0;
                        ctoprim(u_nb, q_nb, gamma);
                        real_t dq_null[20] = {0};
                        if (idim > 0) std::swap(q_nb[1], q_nb[1+idim]);
                        trace(q_nb, dq_null, dtdx, qm_nb, qp_nb, gamma);
                        if (idim > 0) { std::swap(qm_nb[1], qm_nb[1+idim]); std::swap(qp_nb[1], qp_nb[1+idim]); }
                        if (side == 0) { for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = qp_nb[iv-1]; qr_f[iv-1] = get_q(idc, idim, iv); } }
                        else { for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qp(idc, idim, iv); qr_f[iv-1] = qm_nb[iv-1]; } }
                    } else if (id_n <= grid_.ncoarse) {
                        real_t u_nb[20], q_nb[20], dq_nb[20], qm_nb[20], qp_nb[20];
                        for(int iv=1; iv<=grid_.nvar; ++iv) u_nb[iv-1] = grid_.uold(id_n, iv);
                        ctoprim(u_nb, q_nb, gamma);
                        int icn_nb[6]; int ign_nb[7] = {0};
                        grid_.get_nbor_cells(ign_nb, id_n, icn_nb, 0);
                        compute_slopes(id_n, icn_nb, idim, dq_nb, slope_type);
                        if (idim > 0) { std::swap(q_nb[1], q_nb[1+idim]); std::swap(dq_nb[1], dq_nb[1+idim]); }
                        trace(q_nb, dq_nb, dtdx, qm_nb, qp_nb, gamma);
                        if (idim > 0) { std::swap(qm_nb[1], qm_nb[1+idim]); std::swap(qp_nb[1], qp_nb[1+idim]); }
                        if (side == 0) { for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = qp_nb[iv-1]; qr_f[iv-1] = get_q(idc, idim, iv); } }
                        else { for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qp(idc, idim, iv); qr_f[iv-1] = qm_nb[iv-1]; } }
                    } else {
                        if (side == 0) { for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qp(id_n, idim, iv); qr_f[iv-1] = get_q(idc, idim, iv); } }
                        else { for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qp(idc, idim, iv); qr_f[iv-1] = get_q(id_n, idim, iv); } }
                    }
                    if (idim > 0) { std::swap(ql_f[1], ql_f[1+idim]); std::swap(qr_f[1], qr_f[1+idim]); }
                    
                    if (riemann == "hllc") RiemannSolver::solve_hllc(ql_f, qr_f, flux, gamma);
                    else if (riemann == "llf") RiemannSolver::solve_llf(ql_f, qr_f, flux, gamma);
                    else RiemannSolver::solve_hll(ql_f, qr_f, flux, gamma);

                    real_t mass_flux = flux[0];
                    for (int iv = 5; iv < grid_.nvar; ++iv) {
                        real_t q_scalar = (mass_flux > 0) ? ql_f[iv] : qr_f[iv];
                        flux[iv] = mass_flux * q_scalar;
                    }
                    
                    if (idim > 0) std::swap(flux[1], flux[1+idim]);
                    real_t sign = (side == 0) ? 1.0 : -1.0;
                    for(int iv=0; iv<grid_.nvar; ++iv) flux_sum[iv] += sign * flux[iv];
                }
            }
            for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.unew(idc, iv) = grid_.uold(idc, iv) + dtdx * flux_sum[iv-1];
            grid_.unew(idc, 1) = std::max(grid_.unew(idc, 1), 1e-10);
            real_t d_curr = grid_.unew(idc, 1), v2_curr = 0.0;
            for(int i=1; i<=3; ++i) { real_t v = grid_.unew(idc, 1+i)/d_curr; v2_curr += v*v; }
            real_t ei_curr = grid_.unew(idc, 5) - 0.5*d_curr*v2_curr;
            for(int ie=0; ie<nener_; ++ie) ei_curr -= grid_.unew(idc, 6+ie);
            if (ei_curr < d_curr*1e-10/(gamma-1.0)) {
                real_t e_non = 0; for(int ie=0; ie<nener_; ++ie) e_non += grid_.unew(idc, 6+ie);
                grid_.unew(idc, 5) = d_curr*1e-10/(gamma-1.0) + 0.5*d_curr*v2_curr + e_non;
            }
        }
    }
}

void HydroSolver::set_unew(int ilevel) {
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.unew(idc, iv) = grid_.uold(idc, iv);
        }
        igrid = grid_.next[igrid - 1];
    }
}

void HydroSolver::set_uold(int ilevel) {
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son[idc] == 0) {
                for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.uold(idc, iv) = grid_.unew(idc, iv);
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

void HydroSolver::ctoprim(const real_t u[], real_t q[], real_t gamma) {
    real_t d = std::max(u[0], 1e-10); q[0] = d;
    real_t v2 = 0.0; for (int i = 1; i <= 3; ++i) { q[i] = u[i] / d; v2 += q[i] * q[i]; }
    real_t e_thermal_dens = u[4] - 0.5 * d * v2;
    for (int ie = 0; ie < nener_; ++ie) e_thermal_dens -= u[5 + ie];
    q[4] = std::max(e_thermal_dens * (gamma - 1.0), d * 1e-10);
    for (int iv = 5; iv < grid_.nvar; ++iv) q[iv] = u[iv] * (gamma - 1.0);
}

void HydroSolver::compute_slopes(int idc, const int icelln[6], int idim, real_t dq[20], int slope_type) {
    int id_l = icelln[idim * 2], id_r = icelln[idim * 2 + 1];
    if (id_l <= 0) id_l = idc; if (id_r <= 0) id_r = idc;
    real_t ql[20], qc[20], qr[20], u_l[20], u_c[20], u_r[20];
    for(int iv=1; iv<=grid_.nvar; ++iv) { u_l[iv-1]=grid_.uold(id_l, iv); u_c[iv-1]=grid_.uold(idc, iv); u_r[iv-1]=grid_.uold(id_r, iv); }
    ctoprim(u_l, ql, grid_.gamma); ctoprim(u_c, qc, grid_.gamma); ctoprim(u_r, qr, grid_.gamma);
    for (int iv = 0; iv < grid_.nvar; ++iv) {
        real_t dlft = qc[iv] - ql[iv], drgt = qr[iv] - qc[iv];
        if (dlft * drgt <= 0.0) dq[iv] = 0.0;
        else { real_t sgn = (dlft >= 0.0) ? 1.0 : -1.0; dq[iv] = sgn * std::min(std::abs(dlft), std::abs(drgt)); }
    }
}

void HydroSolver::trace(const real_t q[], const real_t dq[], real_t dtdx, real_t qm[], real_t qp[], real_t gamma) {
    real_t r = std::max(q[0], 1e-10), u = q[1], p = std::max(q[4], 1e-10);
    real_t sr0 = -u * dq[0] - dq[1] * r;
    real_t su0 = -u * dq[1] - dq[4] / r;
    for(int ie=0; ie<nener_; ++ie) su0 -= dq[5+ie] / r;
    real_t sp0 = -u * dq[4] - dq[1] * gamma * p;

    for (int iv = 0; iv < grid_.nvar; ++iv) {
        real_t dqi = dq[iv], src = -u * dqi;
        if (iv == 0) src = sr0; else if (iv == 1) src = su0; else if (iv == 4) src = sp0;
        else if (iv >= 5 && iv < 5+nener_) src = -u * dq[iv] - dq[1] * gamma * q[iv];
        qp[iv] = q[iv] - 0.5 * dqi + 0.5 * dtdx * src;
        qm[iv] = q[iv] + 0.5 * dqi + 0.5 * dtdx * src;
    }
    if (qp[0] < 1e-10) qp[0] = q[0]; if (qm[0] < 1e-10) qm[0] = q[0];
}

void HydroSolver::interpol_hydro(const real_t u1[7][20], real_t u2[8][20]) {
    int interpol_type = config_.get_int("refine_params", "interpol_type", 2);
    int nvar = grid_.nvar;
    real_t gam = grid_.gamma;

    if (interpol_type == 0) {
        for (int i = 0; i < 8; ++i) for (int iv = 1; iv <= nvar; ++iv) u2[i][iv-1] = u1[0][iv-1];
        return;
    }

    real_t q1[7][20], q2[8][20];
    for (int i = 0; i < 7; ++i) ctoprim(u1[i], q1[i], gam);

    for (int iv = 0; iv < nvar; ++iv) {
        real_t slopes[3] = {0,0,0};
        for (int idim = 0; idim < NDIM; ++idim) {
            real_t dlft = q1[0][iv] - q1[2*idim+1][iv];
            real_t drgt = q1[2*idim+2][iv] - q1[0][iv];
            if (dlft * drgt <= 0.0) slopes[idim] = 0.0;
            else {
                real_t sgn = (dlft >= 0.0) ? 1.0 : -1.0;
                slopes[idim] = sgn * std::min(std::abs(dlft), std::abs(drgt));
            }
        }
        for (int i = 0; i < 8; ++i) {
            int ix = i & 1, iy = (i & 2) >> 1, iz = (i & 4) >> 2;
            q2[i][iv] = q1[0][iv] + (ix-0.5)*slopes[0] + (iy-0.5)*slopes[1] + (iz-0.5)*slopes[2];
        }
    }

    for (int i = 0; i < 8; ++i) {
        real_t d = q2[i][0], u = q2[i][1], v = q2[i][2], w = q2[i][3], p = q2[i][4];
        u2[i][0] = d; u2[i][1] = d*u; u2[i][2] = d*v; u2[i][3] = d*w;
        real_t e_thermal = p / (gam - 1.0);
        real_t e_kinetic = 0.5 * d * (u*u + v*v + w*w);
        real_t e_nonthermal = 0.0;
        for (int ie = 0; ie < nener_; ++ie) {
            real_t e_rad = q2[i][5+ie] / (gam - 1.0);
            u2[i][5+ie] = e_rad;
            e_nonthermal += e_rad;
        }
        u2[i][4] = e_thermal + e_kinetic + e_nonthermal;
    }
}

real_t HydroSolver::compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor) {
    real_t dt_max = 1e30;
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int id = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            real_t d = std::max(grid_.uold(id, 1), 1e-10);
            real_t v2 = 0.0, v_max = 0.0;
            for (int i = 1; i <= 3; ++i) { real_t v = grid_.uold(id, 1 + i) / d; v2 += v * v; v_max = std::max(v_max, std::abs(v)); }
            real_t p = std::max((grid_.uold(id, 5) - 0.5 * d * v2) * (gamma - 1.0), d * 1e-10);
            real_t cs = std::sqrt(gamma * p / d);
            if (v_max + cs > 1e-20) dt_max = std::min(dt_max, courant_factor * dx / (v_max + cs));
        }
        igrid = grid_.next[igrid - 1];
    }
    return dt_max;
}

void HydroSolver::get_diagnostics(int ilevel, real_t dx, real_t& min_d, real_t& max_v, real_t& min_t, real_t& max_t) {
    min_d = 1e30; max_v = -1e30; min_t = 1e30; max_t = -1e30;
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int id = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            real_t d = grid_.uold(id, 1); min_d = std::min(min_d, d);
            real_t v2 = 0.0; for (int i = 1; i <= 3; ++i) { real_t v = grid_.uold(id, 1 + i) / std::max(d, 1e-10); v2 += v * v; }
            max_v = std::max(max_v, std::sqrt(v2));
            real_t p = (grid_.uold(id, 5) - 0.5 * std::max(d, 1e-10) * v2) * (grid_.gamma - 1.0);
            real_t T = p / std::max(d, 1e-10); min_t = std::min(min_t, T); max_t = std::max(max_t, T);
        }
        igrid = grid_.next[igrid - 1];
    }
}

void HydroSolver::add_gravity_source_terms(int ilevel, real_t dt) {}

} // namespace ramses
