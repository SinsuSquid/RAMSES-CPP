#include "ramses/HydroSolver.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/RiemannSolver.hpp"
#include "ramses/SlopeLimiter.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace ramses {

void HydroSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    int myid = MpiManager::instance().rank() + 1;
    std::vector<int> octs;
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) { octs.push_back(igrid); igrid = grid_.next[igrid - 1]; }
    if (!octs.empty()) godfine1(octs, ilevel, dt, dx);
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
    real_t d = std::max(u[0], 1e-10);
    q[0] = d;
    real_t v2 = 0.0;
    for (int i = 1; i <= 3; ++i) { q[i] = u[i] / d; v2 += q[i] * q[i]; }
    real_t e_int = u[4] - 0.5 * d * v2;
    q[4] = std::max(e_int * (gamma - 1.0), d * 1e-10);
}

void HydroSolver::compute_slopes(int idc, const int icelln[6], int idim, real_t dq[5], int slope_type) {
    int id_l = icelln[idim * 2], id_r = icelln[idim * 2 + 1];
    if (id_l <= 0) id_l = idc; if (id_r <= 0) id_r = idc;

    real_t ql[5], qc[5], qr[5], u_l[5], u_c[5], u_r[5];
    for(int iv=1; iv<=5; ++iv) { u_l[iv-1]=grid_.uold(id_l, iv); u_c[iv-1]=grid_.uold(idc, iv); u_r[iv-1]=grid_.uold(id_r, iv); }
    ctoprim(u_l, ql, grid_.gamma); ctoprim(u_c, qc, grid_.gamma); ctoprim(u_r, qr, grid_.gamma);

    for (int iv = 0; iv < 5; ++iv) {
        real_t dlft = qc[iv] - ql[iv];
        real_t drgt = qr[iv] - qc[iv];
        if (dlft * drgt <= 0.0) dq[iv] = 0.0;
        else {
            real_t sgn = (dlft >= 0.0) ? 1.0 : -1.0;
            dq[iv] = sgn * std::min(std::abs(dlft), std::abs(drgt));
        }
    }
}

void HydroSolver::trace(const real_t q[], const real_t dq[], real_t dtdx, real_t qm[], real_t qp[], real_t gamma) {
    real_t r = std::max(q[0], 1e-10), u = q[1], p = std::max(q[4], 1e-10);
    real_t drx = dq[0], dux = dq[1], dpx = dq[4];
    real_t sr0 = -u * drx - dux * r;
    real_t sp0 = -u * dpx - dux * gamma * p;
    real_t su0 = -u * dux - dpx / r;

    for (int iv = 0; iv < 5; ++iv) {
        real_t dqi = dq[iv], src = 0.0;
        if (iv == 0) src = sr0; else if (iv == 1) src = su0; else if (iv == 4) src = sp0;
        qp[iv] = q[iv] - 0.5 * dqi + 0.5 * dtdx * src;
        qm[iv] = q[iv] + 0.5 * dqi + 0.5 * dtdx * src;
    }
    if (qp[0] < 1e-10) qp[0] = q[0]; if (qm[0] < 1e-10) qm[0] = q[0];
}

void HydroSolver::godfine1(const std::vector<int>& octs, int ilevel, real_t dt, real_t dx) {
    real_t gamma = grid_.gamma;
    real_t dtdx = dt / dx;
    int slope_type = config_.get_int("hydro_params", "slope_type", 1);

    for (int igrid : octs) {
        real_t qm_oct[8][3][5], qp_oct[8][3][5];
        int igridn[7]; grid_.get_nbor_grids(igrid, igridn);

        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            int icelln[6]; grid_.get_nbor_cells(igridn, ic, icelln, igrid);
            real_t qc[5], u_c[5]; for(int iv=1; iv<=5; ++iv) u_c[iv-1] = grid_.uold(idc, iv);
            ctoprim(u_c, qc, gamma);
            for (int idim = 0; idim < NDIM; ++idim) {
                real_t dq[5], qc_rot[5]; for(int iv=0; iv<5; ++iv) qc_rot[iv]=qc[iv];
                compute_slopes(idc, icelln, idim, dq, slope_type);
                if (idim > 0) { std::swap(qc_rot[1], qc_rot[1+idim]); std::swap(dq[1], dq[1+idim]); }
                trace(qc_rot, dq, dtdx, qm_oct[ic-1][idim], qp_oct[ic-1][idim], gamma);
                if (idim > 0) { std::swap(qm_oct[ic-1][idim][1], qm_oct[ic-1][idim][1+idim]); std::swap(qp_oct[ic-1][idim][1], qp_oct[ic-1][idim][1+idim]); }
            }
        }

        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son[idc] != 0) continue;

            int icelln[6]; grid_.get_nbor_cells(igridn, ic, icelln, igrid);
            real_t flux_sum[5] = {0,0,0,0,0};
            for (int idim = 0; idim < NDIM; ++idim) {
                for (int side = 0; side < 2; ++side) {
                    int id_n = icelln[idim * 2 + side];
                    real_t ql_f[5], qr_f[5], flux[5];
                    int ic_n = -1;
                    if (id_n > grid_.ncoarse) {
                        int ig_n = ((id_n - grid_.ncoarse - 1) % grid_.ngridmax) + 1;
                        if (ig_n == igrid) ic_n = (id_n - grid_.ncoarse - 1) / grid_.ngridmax;
                    }

                    if (ic_n >= 0) {
                        if (side == 0) { for(int iv=0; iv<5; ++iv){ ql_f[iv]=qp_oct[ic_n][idim][iv]; qr_f[iv]=qm_oct[ic-1][idim][iv]; } }
                        else { for(int iv=0; iv<5; ++iv){ ql_f[iv]=qm_oct[ic-1][idim][iv]; qr_f[iv]=qp_oct[ic_n][idim][iv]; } }
                    } else {
                        if (id_n <= 0) id_n = idc;
                        int igridn_n[7];
                        if (id_n <= grid_.ncoarse) { for(int i=0; i<7; ++i) igridn_n[i]=0; }
                        else { int ig_n = ((id_n - grid_.ncoarse - 1) % grid_.ngridmax) + 1; grid_.get_nbor_grids(ig_n, igridn_n); }
                        int ic_n_pos = (id_n <= grid_.ncoarse) ? 1 : ((id_n - grid_.ncoarse - 1) / grid_.ngridmax) + 1;
                        int icelln_n[6]; grid_.get_nbor_cells(igridn_n, ic_n_pos, icelln_n, (id_n <= grid_.ncoarse ? 0 : ((id_n - grid_.ncoarse - 1) % grid_.ngridmax) + 1));
                        
                        real_t dq_n[5], qn[5], u_n[5], qm_n[5], qp_n[5];
                        compute_slopes(id_n, icelln_n, idim, dq_n, slope_type);
                        for(int iv=1; iv<=5; ++iv) u_n[iv-1] = grid_.uold(id_n, iv);
                        ctoprim(u_n, qn, gamma);
                        if (idim > 0) { std::swap(qn[1], qn[1+idim]); std::swap(dq_n[1], dq_n[1+idim]); }
                        trace(qn, dq_n, dtdx, qm_n, qp_n, gamma);
                        if (idim > 0) { std::swap(qm_n[1], qm_n[1+idim]); std::swap(qp_n[1], qp_n[1+idim]); }
                        if (side == 0) { for(int iv=0; iv<5; ++iv){ ql_f[iv]=qp_n[iv]; qr_f[iv]=qm_oct[ic-1][idim][iv]; } }
                        else { for(int iv=0; iv<5; ++iv){ ql_f[iv]=qm_oct[ic-1][idim][iv]; qr_f[iv]=qp_n[iv]; } }
                    }

                    if (idim > 0) { std::swap(ql_f[1], ql_f[1+idim]); std::swap(qr_f[1], qr_f[1+idim]); }
                    RiemannSolver::solve_hllc(ql_f, qr_f, flux, gamma);
                    if (idim > 0) { std::swap(flux[1], flux[1+idim]); }
                    real_t sign = (side == 0) ? 1.0 : -1.0;
                    for (int iv = 0; iv < 5; ++iv) flux_sum[iv] += sign * flux[iv];
                }
            }
            for (int iv = 1; iv <= 5; ++iv) grid_.unew(idc, iv) = grid_.uold(idc, iv) + dtdx * flux_sum[iv-1];
            grid_.unew(idc, 1) = std::max(grid_.unew(idc, 1), 1e-10);
            real_t d_n = grid_.unew(idc, 1), v2_n = 0.0;
            for(int i=1; i<=3; ++i) { real_t v = grid_.unew(idc, 1+i)/d_n; v2_n += v*v; }
            real_t ei_n = grid_.unew(idc, 5) - 0.5*d_n*v2_n;
            if (ei_n < d_n*1e-10/(gamma-1.0)) grid_.unew(idc, 5) = d_n*1e-10/(gamma-1.0) + 0.5*d_n*v2_n;
        }
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
            for (int i = 1; i <= 3; ++i) {
                real_t v = grid_.uold(id, 1 + i) / d;
                v2 += v * v; v_max = std::max(v_max, std::abs(v));
            }
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
            real_t d = grid_.uold(id, 1);
            min_d = std::min(min_d, d);
            real_t v2 = 0.0;
            for (int i = 1; i <= 3; ++i) { real_t v = grid_.uold(id, 1 + i) / std::max(d, 1e-10); v2 += v * v; }
            max_v = std::max(max_v, std::sqrt(v2));
            real_t p = (grid_.uold(id, 5) - 0.5 * std::max(d, 1e-10) * v2) * (grid_.gamma - 1.0);
            real_t T = p / std::max(d, 1e-10);
            min_t = std::min(min_t, T); max_t = std::max(max_t, T);
        }
        igrid = grid_.next[igrid - 1];
    }
}

void HydroSolver::add_gravity_source_terms(int ilevel, real_t dt) {}
void HydroSolver::interpol_hydro(const real_t u1[7][20], real_t u2[8][20]) {}

} // namespace ramses
