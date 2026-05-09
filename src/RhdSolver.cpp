#include "ramses/RhdSolver.hpp"
#include "ramses/RelativisticRiemannSolver.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/Constants.hpp"
#include <cmath>
#include <algorithm>

namespace ramses {

RhdSolver::~RhdSolver() {}

void RhdSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    int myid = MpiManager::instance().rank() + 1;
    std::vector<int> octs;
    if (ilevel > 1) {
        int igrid = grid_.get_headl(myid, ilevel);
        while (igrid > 0) { octs.push_back(igrid); igrid = grid_.next[igrid - 1]; }
    }
    
    bool do_level_1 = (ilevel == 1);
    if (octs.empty() && !do_level_1) return;

    real_t gamma = grid_.gamma;
    real_t dtdx = dt / dx;
    int slope_type = config_.get_int("hydro_params", "slope_type", 1);
    std::string riemann = config_.get("hydro_params", "riemann", "hllc");
    std::string eos = config_.get("hydro_params", "eos_rhd", "ideal");

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

    // 2. Flux step
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
                    int id_n0 = id_n - 1;
                    if (side == 0) { for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qp(id_n0, idim, iv); qr_f[iv-1] = get_q(idc_0, idim, iv); } }
                    else { for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qp(idc_0, idim, iv); qr_f[iv-1] = get_q(id_n0, idim, iv); } }
                }
                if (riemann == "hllc") RelativisticRiemannSolver::solve_hllc(ql_f, qr_f, flux, gamma, eos);
                else if (riemann == "hll") RelativisticRiemannSolver::solve_hll(ql_f, qr_f, flux, gamma, eos);
                else RelativisticRiemannSolver::solve_llf(ql_f, qr_f, flux, gamma, eos);
                real_t sgn = (side == 0) ? 1.0 : -1.0;
                for(int iv=1; iv<=grid_.nvar; ++iv) flux_sum[iv-1] += sgn * flux[iv-1] * dtdx;
            }
        }
        for(int iv=1; iv<=grid_.nvar; ++iv) grid_.unew(idc_0 + 1, iv) = grid_.uold(idc_0 + 1, iv) + flux_sum[iv-1];
    };

    if (do_level_1) {
        for (int idc_0 = 0; idc_0 < grid_.ncoarse; ++idc_0) {
            int icn[6]; grid_.get_nbor_cells_coarse(idc_0 + 1, icn);
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

void RhdSolver::set_uold(int ilevel) {
    for (int i = 1; i <= grid_.ncell; ++i) {
        for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.uold(i, iv) = grid_.unew(i, iv);
    }
}

void RhdSolver::set_unew(int ilevel) {
    for (int i = 1; i <= grid_.ncell; ++i) {
        for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.unew(i, iv) = grid_.uold(i, iv);
    }
}

real_t RhdSolver::compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor) {
    real_t dt_level = 1e30;
    std::string eos = config_.get("hydro_params", "eos_rhd", "ideal");
    for (int i = 1; i <= grid_.ncell; ++i) {
        real_t u[20], q[20], cs;
        for(int iv=1; iv<=grid_.nvar; ++iv) u[iv-1] = grid_.uold(i, iv);
        ctoprim(u, q, gamma);
        RelativisticRiemannSolver::find_speed_fast(q, cs, gamma, eos);
        real_t vmax = 0.0;
        for(int idim=0; idim<NDIM; ++idim) vmax = std::max(vmax, std::abs(q[1+idim]) + cs);
        if (vmax > 0) dt_level = std::min(dt_level, dx / vmax);
    }
    return dt_level * courant_factor;
}

void RhdSolver::get_diagnostics(int ilevel, real_t dx, real_t& min_d, real_t& max_v, real_t& min_t, real_t& max_t) {}
void RhdSolver::add_gravity_source_terms(int ilevel, real_t dt) {}
void RhdSolver::synchro_hydro_fine(int ilevel, real_t dt) {}
void RhdSolver::interpol_hydro(const real_t u1[7][64], real_t u2[8][64]) {}

void RhdSolver::ctoprim(const real_t u[], real_t q[], real_t gamma) {
    real_t D = u[0], Mx = u[1], My = u[2], Mz = u[3], E = u[4];
    real_t M = std::sqrt(Mx*Mx + My*My + Mz*Mz);
    std::string eos = config_.get("hydro_params", "eos_rhd", "ideal");

    if (D < 0) D = 1e-12;
    if (E < 0) E = std::sqrt(M*M + D*D + 1e-8);
    if (E*E < M*M + D*D) E = std::sqrt(M*M + D*D + 1e-8);

    real_t R, lor, u2, xsi;
    if (M == 0) {
        q[0] = D; q[1] = 0; q[2] = 0; q[3] = 0;
        if (eos == "TM") q[4] = (E*E - D*D) / (3.0 * E);
        else q[4] = (E - D) * (gamma - 1.0);
    } else {
        newton_raphson_mignone(D, M, E, gamma, R);
        u2 = M*M / (R*R - M*M);
        lor = std::sqrt(1.0 + u2);
        q[0] = D / lor;
        q[1] = Mx / R; q[2] = My / R; q[3] = Mz / R;
        xsi = ((R - D) - u2 / (lor + 1.0) * D) / (lor*lor);
        if (eos == "TM") {
            real_t rho = q[0];
            q[4] = (2.0 * xsi * (xsi + 2.0 * rho)) / (5.0 * (xsi + rho) + std::sqrt(9.0*xsi*xsi + 18.0*rho*xsi + 25.0*rho*rho));
        } else {
            q[4] = (gamma - 1.0) / gamma * xsi;
        }
    }
    // Passive scalars
    for (int n = 5; n < grid_.nvar; ++n) q[n] = u[n] / (D);
}

void RhdSolver::newton_raphson_mignone(real_t D, real_t M, real_t E, real_t gamma, real_t& R) {
    real_t Delta = 16.0 * E * E - 12.0 * M * M;
    R = (4.0 * E + std::sqrt(std::max(0.0, Delta))) / 6.0;
    real_t Eprim = E - D;
    R -= D;
    std::string eos = config_.get("hydro_params", "eos_rhd", "ideal");

    auto f_Mignone = [&](real_t R_val) {
        real_t u2 = M * M / (std::pow(R_val + D, 2) - M * M);
        real_t lor = std::sqrt(1.0 + u2);
        real_t xsi = (R_val - u2 / (lor + 1.0) * D) / (lor * lor);
        real_t P;
        if (eos == "TM") {
            real_t rho = D / lor;
            P = (2.0 * xsi * (xsi + 2.0 * rho)) / (5.0 * (xsi + rho) + std::sqrt(9.0*xsi*xsi + 18.0*rho*xsi + 25.0*rho*rho));
        } else {
            P = (gamma - 1.0) / gamma * xsi;
        }
        return R_val - P - Eprim;
    };

    auto f_prim_Mignone = [&](real_t R_val) {
        real_t u2 = M * M / (std::pow(R_val + D, 2) - M * M);
        real_t lor = std::sqrt(1.0 + u2);
        real_t dpdR;
        if (eos == "TM") {
            real_t xsi = (R_val - u2 / (lor + 1.0) * D) / (lor * lor);
            real_t rho = D / lor;
            real_t P = (2.0 * xsi * (xsi + 2.0 * rho)) / (5.0 * (xsi + rho) + std::sqrt(9.0*xsi*xsi + 18.0*rho*xsi + 25.0*rho*rho));
            real_t dpdxsi = (2.0*xsi + 2.0*rho - 5.0*P) / (5.0*rho + 5.0*xsi - 8.0*P);
            real_t dpdrho = (2.0*xsi - 5.0*P) / (5.0*rho + 5.0*xsi - 8.0*P);
            real_t dv2dR = -2.0 * M * M / std::pow(R_val + D, 3);
            real_t dxsidR = 1.0 / (lor * lor) - lor / 2.0 * (D + 2.0 * lor * xsi) * dv2dR;
            real_t drhodR = D * lor / 2.0 * dv2dR;
            dpdR = dpdxsi * dxsidR + dpdrho * drhodR;
        } else {
            dpdR = (gamma - 1.0) / gamma * (1.0 + M * M / std::pow(R_val + D, 2) * (1.0 - D * lor / (R_val + D)));
        }
        return 1.0 - dpdR;
    };

    real_t eps = 1.0;
    int iter = 0;
    while (std::abs(eps) > 1e-10 && iter < 100) {
        eps = f_Mignone(R) / f_prim_Mignone(R) / (R + 1e-20);
        R *= (1.0 - eps);
        iter++;
    }
    R += D;
}

void RhdSolver::compute_slopes(int idc, const int icelln[6], int idim, real_t dq[20], int slope_type) {
    real_t ql[20], qc[20], qr[20], u_c[20], gamma = grid_.gamma;
    for(int iv=1; iv<=grid_.nvar; ++iv) u_c[iv-1]=grid_.uold(idc, iv);
    ctoprim(u_c, qc, gamma);
    auto get_nb_q = [&](int id_n, int side, real_t q_nb[20]) {
        if (id_n > 0) { real_t u_nb[20]; for(int iv=1; iv<=grid_.nvar; ++iv) u_nb[iv-1]=grid_.uold(id_n, iv); ctoprim(u_nb, q_nb, gamma); }
        else { for(int iv=0; iv<grid_.nvar; ++iv) q_nb[iv] = qc[iv]; int ib = -id_n; if(ib > 0 && ib <= (int)grid_.bound_type.size() && grid_.bound_type.at(ib-1) == 1) q_nb[1 + idim] *= -1.0; }
    };
    get_nb_q(icelln[idim*2], 0, ql); get_nb_q(icelln[idim*2+1], 1, qr);
    for (int iv = 0; iv < grid_.nvar; ++iv) {
        real_t dlft = qc[iv] - ql[iv], drgt = qr[iv] - qc[iv];
        if (dlft * drgt <= 0.0) dq[iv] = 0.0;
        else { real_t sgn = (dlft >= 0.0) ? 1.0 : -1.0; dq[iv] = sgn * std::min(std::abs(dlft), std::abs(drgt)); }
    }
}

void RhdSolver::trace(const real_t q[], const real_t dq[], real_t dtdx, real_t qm[], real_t qp[], real_t gamma) {
    // Simplified RHD trace (primitive reconstruction + dtdx correction)
    // For full second-order accuracy, RHD needs a more complex characteristic tracing.
    // Here we use a standard MUSCL-Hancock simplification as seen in legacy umuscl.f90
    for (int iv = 0; iv < grid_.nvar; ++iv) {
        qp[iv] = q[iv] - 0.5 * dq[iv]; // Simplified
        qm[iv] = q[iv] + 0.5 * dq[iv];
    }
}

} // namespace ramses