#include "ramses/solvers/hydro/HydroSolver.hpp"
#include "ramses/core/MpiManager.hpp"
#include "ramses/solvers/hydro/RiemannSolver.hpp"
#include "ramses/core/Parameters.hpp"
#include "ramses/core/Constants.hpp"
#include "ramses/core/SolverFactory.hpp"
#include "ramses/solvers/hydro/SlopeLimiter.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

namespace ramses {

HydroSolver::~HydroSolver() {}

void HydroSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    int myid = MpiManager::instance().rank() + 1;
    bool verbose = config_.get_bool("run_params", "verbose", false);
    
    std::vector<int> cell_levels(grid_.ncell + 1, 0);
    for (int ig = 1; ig <= grid_.ngridmax; ++ig) {
        if (grid_.father[ig - 1] > 0) {
            int level = 1;
            int curr_ig = ig;
            while (curr_ig > 0) {
                int father_cell = grid_.father[curr_ig - 1];
                if (father_cell <= grid_.ncoarse) break;
                curr_ig = ((father_cell - grid_.ncoarse - 1) % grid_.ngridmax) + 1;
                level++;
            }
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
                cell_levels[idc] = level;
            }
        }
    }
    if (verbose && ilevel == 5) {
        std::cout << "[DEBUG verify] cell_levels[1057] = " << cell_levels[1057] 
                  << " get_cell_level(1057) = " << grid_.get_cell_level(1057) 
                  << " ncell = " << grid_.ncell << std::endl;
    }

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
    // std::cout << "[DEBUG slope] slope_type = " << slope_type << " config_raw = " << config_.get("hydro_params", "slope_type", "NONE") << std::endl;
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
                bool is_refined = (id_n > 0 && (grid_.son.at(id_n - 1) > 0 || cell_levels[id_n] > ilevel));
                if (is_refined) {
                    // Refined interface: set flux to zero as it will be updated at the finer level.
                    std::fill(flux, flux + grid_.nvar, 0.0);
                } else {
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
                        if (id_n > 0 && !grid_.son[id_n0] && cell_levels[id_n] == ilevel) {
                            // Neighbor is at the same level; use pre-computed traces
                            if (verbose) {
                                static int branch_same = 0;
                                if (branch_same < 10 && ilevel == 4) {
                                    std::cout << "[DEBUG branch] SAME level for idc_0=" << idc_0 << " id_n=" << id_n << " side=" << side << std::endl;
                                    branch_same++;
                                }
                            }
                            if (side == 0) { // Cell i, interface i-1/2. Neighbor is i-1.
                                for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qr(id_n0, idim, iv); qr_f[iv-1] = get_ql(idc_0, idim, iv); }
                            } else { // Cell i, interface i+1/2. Neighbor is i+1.
                                for(int iv=1; iv<=grid_.nvar; ++iv) { ql_f[iv-1] = get_qr(idc_0, idim, iv); qr_f[iv-1] = get_ql(id_n0, idim, iv); }
                            }
                        } else {
                            // Interpolate dynamically from coarse neighbor id_n (at level ilevel - 1)
                            if (verbose) {
                                static int branch_coarse = 0;
                                if (branch_coarse < 10 && ilevel == 4) {
                                    std::cout << "[DEBUG branch] COARSE level for idc_0=" << idc_0 << " id_n=" << id_n << " side=" << side << " cell_levels[id_n]=" << (id_n > 0 ? cell_levels[id_n] : -1) << std::endl;
                                    branch_coarse++;
                                }
                            }
                            real_t q_nb[20] = {0};
                            
                            // 1. Get center primitive variables
                            real_t u_c[20], q_c[20];
                            for(int iv=1; iv<=grid_.nvar; ++iv) u_c[iv-1] = grid_.uold(id_n, iv);
                            ctoprim(u_c, q_c, gamma);
                            
                            // 2. Get neighbor cells of id_n at level ilevel - 1
                            int icn_p[6] = {0};
                            if (id_n <= grid_.ncoarse) {
                                grid_.get_nbor_cells_coarse(id_n, icn_p);
                            } else {
                                int ig_p = ((id_n - grid_.ncoarse - 1) % grid_.ngridmax) + 1;
                                int ic_p = ((id_n - grid_.ncoarse - 1) / grid_.ngridmax) + 1;
                                int ign_p[7]; grid_.get_nbor_grids(ig_p, ign_p);
                                grid_.get_nbor_cells(ign_p, ic_p, icn_p, ig_p);
                            }
                            
                            // 3. Get neighbor variables
                            int id_l = icn_p[idim*2];
                            int id_r = icn_p[idim*2+1];
                            real_t u_l[20] = {0}, u_r[20] = {0};
                            
                            auto get_nb_u_coarse = [&](int id, int side, real_t u_out[20]) {
                                if (id > 0) {
                                    for(int iv=1; iv<=grid_.nvar; ++iv) u_out[iv-1] = grid_.uold(id, iv);
                                } else {
                                    for(int iv=0; iv<grid_.nvar; ++iv) u_out[iv] = u_c[iv];
                                    int ib = -id;
                                    if (ib > 0 && ib <= (int)grid_.bound_type.size() && grid_.bound_type[ib - 1] == 1) {
                                        u_out[1 + idim] *= -1.0;
                                    }
                                }
                            };
                            
                            get_nb_u_coarse(id_l, 0, u_l);
                            get_nb_u_coarse(id_r, 1, u_r);
                            
                            // Swap velocity components for 1D slope computation if idim > 0
                            if (idim > 0) {
                                std::swap(u_l[1], u_l[1+idim]);
                                std::swap(u_c[1], u_c[1+idim]);
                                std::swap(u_r[1], u_r[1+idim]);
                            }
                            
                            int interpol_type = config_.get_int("refine_params", "interpol_type", 1);
                            int interpol_var = config_.get_int("refine_params", "interpol_var", config_.get_int("hydro_params", "interpol_var", 0));
                            real_t dq_coarse[20] = {0};
                            real_t dx_coarse = 2.0 * dx; // Coarse neighbor is 2x larger
                            
                            if (interpol_var == 0) {
                                // Conserved variables interpolation
                                real_t slopes[20] = {0};
                                for (int iv = 0; iv < grid_.nvar; ++iv) {
                                    real_t dlft = u_c[iv] - u_l[iv];
                                    real_t drgt = u_r[iv] - u_c[iv];
                                    if (interpol_type > 0 && dlft * drgt > 0.0) {
                                        real_t sgn = (dlft >= 0.0) ? 1.0 : -1.0;
                                        if (interpol_type == 1) {
                                            slopes[iv] = sgn * 0.5 * std::min(std::abs(dlft), std::abs(drgt));
                                        } else if (interpol_type == 2) {
                                            slopes[iv] = sgn * std::min({std::abs(dlft), std::abs(drgt), 0.25 * std::abs(dlft + drgt)});
                                        } else if (interpol_type == 3) {
                                            slopes[iv] = 0.25 * (dlft + drgt);
                                        }
                                    }
                                    dq_coarse[iv] = slopes[iv] / dx_coarse;
                                }
                                
                                // Swap back conserved vars
                                if (idim > 0) {
                                    std::swap(u_l[1], u_l[1+idim]);
                                    std::swap(u_c[1], u_c[1+idim]);
                                    std::swap(u_r[1], u_r[1+idim]);
                                }
                                
                                // Re-evaluate primitive q_c for trace
                                ctoprim(u_c, q_c, gamma);
                                if (idim > 0) {
                                    std::swap(q_c[1], q_c[1+idim]);
                                }
                                
                                // Compute primitive slopes for the trace step
                                real_t q_l[20], q_r[20];
                                ctoprim(u_l, q_l, gamma);
                                ctoprim(u_r, q_r, gamma);
                                if (idim > 0) {
                                    std::swap(q_l[1], q_l[1+idim]);
                                    std::swap(q_r[1], q_r[1+idim]);
                                }
                                for (int iv = 0; iv < grid_.nvar; ++iv) {
                                    real_t dlft = q_c[iv] - q_l[iv];
                                    real_t drgt = q_r[iv] - q_c[iv];
                                    real_t slope_p = 0.0;
                                    if (interpol_type > 0 && dlft * drgt > 0.0) {
                                        real_t sgn = (dlft >= 0.0) ? 1.0 : -1.0;
                                        if (interpol_type == 1) {
                                            slope_p = sgn * 0.5 * std::min(std::abs(dlft), std::abs(drgt));
                                        } else if (interpol_type == 2) {
                                            slope_p = sgn * std::min({std::abs(dlft), std::abs(drgt), 0.25 * std::abs(dlft + drgt)});
                                        } else if (interpol_type == 3) {
                                            slope_p = 0.25 * (dlft + drgt);
                                        }
                                    }
                                    dq_coarse[iv] = slope_p / dx_coarse;
                                }
                            } else {
                                // Primitive variables interpolation
                                if (idim > 0) {
                                    std::swap(u_l[1], u_l[1+idim]);
                                    std::swap(u_c[1], u_c[1+idim]);
                                    std::swap(u_r[1], u_r[1+idim]);
                                }
                                ctoprim(u_c, q_c, gamma);
                                real_t q_l[20], q_r[20];
                                ctoprim(u_l, q_l, gamma);
                                ctoprim(u_r, q_r, gamma);
                                if (idim > 0) {
                                    std::swap(q_l[1], q_l[1+idim]);
                                    std::swap(q_c[1], q_c[1+idim]);
                                    std::swap(q_r[1], q_r[1+idim]);
                                }
                                for (int iv = 0; iv < grid_.nvar; ++iv) {
                                    real_t dlft = q_c[iv] - q_l[iv];
                                    real_t drgt = q_r[iv] - q_c[iv];
                                    real_t slope = 0.0;
                                    if (interpol_type > 0 && dlft * drgt > 0.0) {
                                        real_t sgn = (dlft >= 0.0) ? 1.0 : -1.0;
                                        if (interpol_type == 1) {
                                            slope = sgn * 0.5 * std::min(std::abs(dlft), std::abs(drgt));
                                        } else if (interpol_type == 2) {
                                            slope = sgn * std::min({std::abs(dlft), std::abs(drgt), 0.25 * std::abs(dlft + drgt)});
                                        } else if (interpol_type == 3) {
                                            slope = 0.25 * (dlft + drgt);
                                        }
                                    }
                                    dq_coarse[iv] = slope / dx_coarse;
                                }
                            }
                            
                            real_t qm_nb[20], qp_nb[20];
                            trace(q_c, dq_coarse, dt, dx_coarse, qm_nb, qp_nb, gamma);
                            
                            // Swap back if idim > 0
                            if (idim > 0) {
                                std::swap(qm_nb[1], qm_nb[1+idim]);
                                std::swap(qp_nb[1], qp_nb[1+idim]);
                            }
                            
                            if (side == 0) { 
                                for(int iv=0; iv<grid_.nvar; ++iv) { ql_f[iv] = qm_nb[iv]; qr_f[iv] = get_ql(idc_0, idim, iv+1); } 
                            } else { 
                                for(int iv=0; iv<grid_.nvar; ++iv) { ql_f[iv] = get_qr(idc_0, idim, iv+1); qr_f[iv] = qp_nb[iv]; } 
                            }
                        }
                    }
                    if (idim > 0) { std::swap(ql_f[1], ql_f[1+idim]); std::swap(qr_f[1], qr_f[1+idim]); }
                    std::string riemann_type = config_.get("hydro_params", "riemann", "llf");
                    if (riemann_type == "hllc") RiemannSolver::solve_hllc(ql_f, qr_f, flux, gamma, nener_, grid_.gamma_rad);
                    else if (riemann_type == "hll") RiemannSolver::solve_hll(ql_f, qr_f, flux, gamma, nener_, grid_.gamma_rad);
                    else if (riemann_type == "exact") RiemannSolver::solve_godunov_nr(ql_f, qr_f, flux, gamma);
                    else if (riemann_type == "acoustic") RiemannSolver::solve_acoustic(ql_f, qr_f, flux, gamma, nener_, grid_.gamma_rad);
                    else RiemannSolver::solve_llf(ql_f, qr_f, flux, gamma, nener_, grid_.gamma_rad);

                    real_t mass_flux = flux[0];
                    real_t u_interface = mass_flux / std::max((mass_flux > 0) ? ql_f[0] : qr_f[0], (real_t)1e-10);
                    for (int iv = NDIM + 2; iv < grid_.nvar; ++iv) {
                        if (iv < NDIM + 2 + nener_) {
                            int ie = iv - (NDIM + 2);
                            real_t p_rad = (mass_flux > 0) ? ql_f[iv] : qr_f[iv];
                            flux[iv] = u_interface * (p_rad / (grid_.gamma_rad[ie] - 1.0));
                        } else {
                            flux[iv] = mass_flux * ((mass_flux > 0) ? ql_f[iv] : qr_f[iv]);
                        }
                    }
                    if (idim > 0) std::swap(flux[1], flux[1+idim]);
                    // Refluxing: update coarser neighbor cell's conservative variables
                    if (id_n > 0 && cell_levels[id_n] < ilevel) {
                        real_t factor = 1.0 / (1 << NDIM);
                        real_t sign = (side == 0) ? 1.0 : -1.0;
                        for (int iv = 1; iv <= grid_.nvar; ++iv) {
                            grid_.unew(id_n, iv) -= dtdx * sign * flux[iv - 1] * factor;
                        }
                    }
                }
                
                real_t sign = (side == 0) ? 1.0 : -1.0;
                for(int iv=0; iv<grid_.nvar; ++iv) flux_sum[iv] += sign * flux[iv];
            }
        }
        for (int iv = 1; iv <= grid_.nvar; ++iv) {
            grid_.unew(idc_0 + 1, iv) = grid_.unew(idc_0 + 1, iv) + dtdx * flux_sum[iv-1];
        }
        if (nener_ > 0) {
            real_t divu = 0.0;
            real_t d_c = std::max(grid_.uold(idc_0 + 1, 1), (real_t)1e-10);
            for (int idim = 0; idim < NDIM; ++idim) {
                int id_l = icn[idim * 2];
                int id_r = icn[idim * 2 + 1];
                
                real_t v_l = 0.0;
                real_t dx_l = dx;
                if (id_l > 0) {
                    real_t d_l = std::max(grid_.uold(id_l, 1), (real_t)1e-10);
                    v_l = grid_.uold(id_l, 2 + idim) / d_l;
                } else {
                    v_l = grid_.uold(idc_0 + 1, 2 + idim) / d_c;
                    int ib = -id_l;
                    if (ib > 0 && ib <= (int)grid_.bound_type.size() && grid_.bound_type[ib - 1] == 1) {
                        v_l *= -1.0;
                    }
                    dx_l = dx * 1.5;
                }
                
                real_t v_r = 0.0;
                real_t dx_r = dx;
                if (id_r > 0) {
                    real_t d_r = std::max(grid_.uold(id_r, 1), (real_t)1e-10);
                    v_r = grid_.uold(id_r, 2 + idim) / d_r;
                } else {
                    v_r = grid_.uold(idc_0 + 1, 2 + idim) / d_c;
                    int ib = -id_r;
                    if (ib > 0 && ib <= (int)grid_.bound_type.size() && grid_.bound_type[ib - 1] == 1) {
                        v_r *= -1.0;
                    }
                    dx_r = dx * 1.5;
                }
                
                divu += (v_r - v_l) / (dx_l + dx_r);
            }
            int iener = NDIM + 2;
            for (int ie = 0; ie < nener_; ++ie) {
                grid_.unew(idc_0 + 1, iener + 1 + ie) -= (grid_.gamma_rad[ie] - 1.0) * grid_.uold(idc_0 + 1, iener + 1 + ie) * divu * dt;
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
        if (verbose) {
            static int print_count = 0;
            if (print_count < 32 && ilevel == 4) {
                std::cout << "[DEBUG flux_sum] idc=" << (idc_0 + 1) << " level=" << ilevel
                          << " flux_sum[0]=" << flux_sum[0] << " flux_sum[1]=" << flux_sum[1] << " flux_sum[2]=" << flux_sum[2] << std::endl;
                print_count++;
            }
        }
    };

    if (do_level_0) {
        for (int idc_0 = 0; idc_0 < grid_.ncoarse; ++idc_0) {
            int icn[6];
            grid_.get_nbor_cells_coarse(idc_0 + 1, icn);
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
    
    for (int iv = iener + 1; iv < grid_.nvar; ++iv) { if (iv < iener + 1 + nener_) q[iv] = u[iv] * (grid_.gamma_rad[iv - (iener + 1)] - 1.0); else q[iv] = u[iv] / d; }
}

void HydroSolver::compute_slopes(int idc, const int icelln[6], int idim, real_t dq[20], int slope_type, real_t dt, real_t dx) {
    if (slope_type == 0) {
        for (int iv = 0; iv < grid_.nvar; ++iv) dq[iv] = 0.0;
        return;
    }
    real_t ql[20], qc[20], qr[20], u_c[20];
    for(int iv=1; iv<=grid_.nvar; ++iv) u_c[iv-1]=grid_.uold(idc, iv);
    ctoprim(u_c, qc, grid_.gamma);
    auto get_nb_q = [&](int id_n, int side, real_t q_nb[20]) {
        if (id_n > 0) { real_t u_nb[20]; for(int iv=1; iv<=grid_.nvar; ++iv) u_nb[iv-1]=grid_.uold(id_n, iv); ctoprim(u_nb, q_nb, grid_.gamma); }
        else { for(int iv=0; iv<grid_.nvar; ++iv) q_nb[iv] = qc[iv]; int ib = -id_n; if(ib > 0 && ib <= (int)grid_.bound_type.size() && grid_.bound_type.at(ib-1) == 1) q_nb[1 + idim] *= -1.0; }
    };
    get_nb_q(icelln[idim*2], 0, ql); get_nb_q(icelln[idim*2+1], 1, qr);

    real_t nu = qc[1 + idim] * dt / dx; // Courant number for Superbee/Ultrabee

    for (int iv = 0; iv < grid_.nvar; ++iv) {
        if (slope_type == 5 && iv != 0) {
            dq[iv] = 0.0; // Ultrabee applies to density only
        } else if (slope_type == 6 && iv != 0) {
            dq[iv] = 0.0; // Centered slope type 6 applies to density only
        } else {
            dq[iv] = SlopeLimiter::compute_slope(ql[iv], qc[iv], qr[iv], slope_type, 1.5, nu) / dx;
        }
    }
}

void HydroSolver::trace(const real_t q[], const real_t dq[], real_t dt, real_t dx, real_t qm[], real_t qp[], real_t gamma) {
    /**
     * MUSCL-Hancock Trace Step for Euler equations in primitive variables W = (rho, u, p)^T.
     * 
     * 1. Governing equations for the time evolution of W:
     *    dW/dt + A(W) * dW/dx = S(W)
     * 
     *    where the primitive Jacobian A(W) is:
     *    A(W) = [   u      rho       0    ]
     *           [   0       u      1/rho  ]
     *           [   0    rho*cs^2    u    ]
     * 
     * 2. Primitive time derivatives:
     *    d(rho)/dt = - u * d(rho)/dx - rho * d(u)/dx
     *                (corresponds to sr0, with dq[0] = d(rho)/dx and dq[1] = d(u)/dx)
     * 
     *    d(u)/dt   = - u * d(u)/dx - (1/rho) * d(p)/dx - (1/rho) * sum( d(p_e)/dx )
     *                (corresponds to su0, with dq[NDIM+1] = d(p)/dx and extra thermal/non-thermal components p_e)
     * 
     *    d(p)/dt   = - u * d(p)/dx - rho * cs^2 * d(u)/dx
     *                (corresponds to sp0, with dq[NDIM+1] = d(p)/dx)
     * 
     * 3. Reconstruct states at cell interfaces for a half time step dt/2 (MUSCL-Hancock):
     *    W_L (qp) = W - 0.5 * dx * dW/dx + 0.5 * dt * dW/dt   (left interface: x - dx/2)
     *    W_R (qm) = W + 0.5 * dx * dW/dx + 0.5 * dt * dW/dt   (right interface: x + dx/2)
     */
    real_t r = std::max(q[0], 1e-10), u = q[1], p = std::max(q[NDIM+1], r * 1e-20);
    real_t sr0 = -u * dq[0] - dq[1] * r;
    real_t su0 = -u * dq[1] - dq[NDIM+1] / r;
    int iener = NDIM + 1;
    for(int ie=0; ie<nener_; ++ie) su0 -= dq[iener+1+ie] / r;
    
    real_t cs2 = RiemannSolver::get_cs2(r, p, gamma, q, nener_, grid_.gamma_rad);
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

        // Ensure states don't explode
        if (iv == 0) {
            qp[iv] = std::max(r * (real_t)1e-2, std::min(r * (real_t)1e2, qp[iv]));
            qm[iv] = std::max(r * (real_t)1e-2, std::min(r * (real_t)1e2, qm[iv]));
        } else if (iv == NDIM + 1) {
            qp[iv] = std::max(p * (real_t)1e-2, std::min(p * (real_t)1e2, qp[iv]));
            qm[iv] = std::max(p * (real_t)1e-2, std::min(p * (real_t)1e2, qm[iv]));
        } else if (iv > 0 && iv <= NDIM) {
            qp[iv] = std::clamp(qp[iv], (real_t)-1e3, (real_t)1e3);
            qm[iv] = std::clamp(qm[iv], (real_t)-1e3, (real_t)1e3);
        }
    }
    // No trace prints here
}

void HydroSolver::interpol_hydro(const real_t u1[7][64], real_t u2[8][64]) {
    int nvar = grid_.nvar; real_t gam = grid_.gamma;
    int interpol_var = config_.get_int("refine_params", "interpol_var", config_.get_int("hydro_params", "interpol_var", 0));
    int interpol_type = config_.get_int("refine_params", "interpol_type", 1);

    real_t q1[7][64] = {0};
    real_t q2[8][64] = {0};

    // Step 1: Convert conserved variables u1 to interpolation variables q1
    if (interpol_var == 0) {
        // Mode 0: conserved variables directly (rho, momentum, E, scalars)
        for (int i = 0; i < 7; ++i) {
            for (int iv = 0; iv < nvar; ++iv) {
                q1[i][iv] = u1[i][iv];
            }
        }
    } else if (interpol_var == 1) {
        // Mode 1: density, momentum, internal energy, scalars
        for (int i = 0; i < 7; ++i) {
            real_t d = std::max(u1[i][0], (real_t)1e-10);
            q1[i][0] = u1[i][0]; // density
            real_t e_kin = 0.0;
            for (int idim = 1; idim <= NDIM; ++idim) {
                q1[i][idim] = u1[i][idim]; // momentum
                e_kin += 0.5 * u1[i][idim] * u1[i][idim] / d;
            }
            real_t e_nonthermal = 0.0;
            int iener = NDIM + 1;
            for (int ie = 0; ie < nener_; ++ie) {
                e_nonthermal += u1[i][iener + 1 + ie];
                q1[i][iener + 1 + ie] = u1[i][iener + 1 + ie];
            }
            q1[i][iener] = u1[i][iener] - e_kin - e_nonthermal; // internal energy
            for (int iv = iener + 1 + nener_; iv < nvar; ++iv) {
                q1[i][iv] = u1[i][iv]; // passive scalars/other vars
            }
        }
    } else {
        // Mode 2: primitive variables (density, velocity, pressure, scalars)
        for (int i = 0; i < 7; ++i) {
            ctoprim(u1[i], q1[i], gam);
        }
    }

    // Step 2: Perform linear interpolation on q1 -> q2
    for (int iv = 0; iv < nvar; ++iv) {
        real_t slopes[3] = {0,0,0};
        for (int idim = 0; idim < NDIM; ++idim) {
            real_t dlft = q1[0][iv] - q1[2*idim+1][iv];
            real_t drgt = q1[2*idim+2][iv] - q1[0][iv];
            if (dlft * drgt > 0.0) {
                real_t sgn = (dlft >= 0.0) ? 1.0 : -1.0;
                int itype = interpol_type;
                if (interpol_type == 4) {
                    // type 3 (central) for velocity, type 2 (moncen) for others
                    if (interpol_var == 2 && iv >= 1 && iv <= NDIM) {
                        itype = 3;
                    } else {
                        itype = 2;
                    }
                }
                if (itype == 1) {
                    slopes[idim] = sgn * 0.5 * std::min(std::abs(dlft), std::abs(drgt));
                } else if (itype == 2) {
                    slopes[idim] = sgn * std::min({std::abs(dlft), std::abs(drgt), 0.25 * std::abs(dlft + drgt)});
                } else if (itype == 3) {
                    slopes[idim] = 0.25 * (dlft + drgt);
                } else {
                    slopes[idim] = 0.0;
                }
            }
        }
        for (int i = 0; i < 8; ++i) {
            int ix = i & 1, iy = (i & 2) >> 1, iz = (i & 4) >> 2;
            q2[i][iv] = q1[0][iv] + (ix-0.5)*slopes[0] + (iy-0.5)*slopes[1] + (iz-0.5)*slopes[2];
        }
    }

    // Step 3: Convert interpolation variables q2 back to conserved variables u2
    if (interpol_var == 0) {
        // Mode 0: conserved variables directly
        for (int i = 0; i < 8; ++i) {
            for (int iv = 0; iv < nvar; ++iv) {
                u2[i][iv] = q2[i][iv];
            }
        }
    } else if (interpol_var == 1) {
        // Mode 1: density, momentum, internal energy -> density, momentum, total energy
        for (int i = 0; i < 8; ++i) {
            real_t d = std::max(q2[i][0], (real_t)1e-10);
            u2[i][0] = d;
            real_t e_kin = 0.0;
            for (int idim = 1; idim <= NDIM; ++idim) {
                u2[i][idim] = q2[i][idim];
                e_kin += 0.5 * q2[i][idim] * q2[i][idim] / d;
            }
            real_t e_nonthermal = 0.0;
            int iener = NDIM + 1;
            for (int ie = 0; ie < nener_; ++ie) {
                u2[i][iener + 1 + ie] = q2[i][iener + 1 + ie];
                e_nonthermal += q2[i][iener + 1 + ie];
            }
            u2[i][iener] = q2[i][iener] + e_kin + e_nonthermal;
            for (int iv = iener + 1 + nener_; iv < nvar; ++iv) {
                u2[i][iv] = q2[i][iv];
            }
        }
    } else {
        // Mode 2: primitive variables -> conserved variables (with momentum conservation)
        int n2d = 1 << NDIM;
        real_t oneover_n2d = 1.0 / n2d;
        for (int d_idx = 1; d_idx <= NDIM; ++d_idx) {
            real_t mom_sum = 0.0;
            for (int i = 0; i < n2d; ++i) {
                mom_sum += q2[i][0] * q2[i][d_idx] * oneover_n2d;
            }
            real_t mom_err = mom_sum - u1[0][d_idx];
            for (int i = 0; i < n2d; ++i) {
                real_t d = std::max(q2[i][0], (real_t)1e-10);
                q2[i][d_idx] = (d * q2[i][d_idx] - mom_err) / d;
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
            for (int ie = 0; ie < nener_; ++ie) { real_t e_rad = q2[i][iener+1+ie] / (grid_.gamma_rad[ie] - 1.0); u2[i][iener+1+ie] = e_rad; e_nonthermal += e_rad; }
            u2[i][iener] = e_thermal + e_kinetic + e_nonthermal;
            for (int iv = iener + 1 + nener_; iv < nvar; ++iv) u2[i][iv] = d * q2[i][iv];
        }
    }
}

real_t HydroSolver::compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor) {
    int myid = MpiManager::instance().rank() + 1;
    real_t dt_max = params::boxlen / 1e-10;
    real_t max_dtdx = 0.1;
    real_t max_dt = max_dtdx * dx;

    auto get_cs = [&](real_t d, real_t p, real_t sum_gamma_p) {
        real_t cs2 = RiemannSolver::get_cs2(d, p, gamma);
        cs2 += sum_gamma_p / d;
        return std::sqrt(std::max(cs2, (real_t)1e-10));
    };

    if (ilevel == 0) {
        for (int id = 1; id <= grid_.ncoarse; ++id) {
            if (grid_.son.at(id - 1) > 0) continue;
            real_t d = std::max(grid_.uold(id, 1), 1e-10), v2 = 0.0;
            for (int i = 1; i <= NDIM; ++i) { real_t v = grid_.uold(id, 1 + i) / d; v2 += v * v; }
            int iener = NDIM + 2;
            
            real_t p = 0;
            real_t sum_gamma_p = 0.0;
            if (params::barotropic_eos) {
                real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
                real_t temp_mu = T2;
                if (params::barotropic_eos_form == "polytrope") temp_mu = T2 * std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
                else if (params::barotropic_eos_form == "double_polytrope") temp_mu = T2 * (1.0 + std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0));
                p = d * temp_mu / (params::units_velocity * params::units_velocity);
            } else {
                real_t e_nonthermal = 0.0;
                for (int ie = 0; ie < nener_; ++ie) {
                    real_t e_rad = grid_.uold(id, iener + 1 + ie);
                    e_nonthermal += e_rad;
                    sum_gamma_p += grid_.gamma_rad[ie] * (grid_.gamma_rad[ie] - 1.0) * e_rad;
                }
                p = std::max((grid_.uold(id, iener) - 0.5 * d * v2 - e_nonthermal) * (gamma - 1.0), d * 1e-10);
                p = std::min(p, d * 1e2);
            }
            
            real_t cs = get_cs(d, p, sum_gamma_p);
            real_t u_wave = NDIM * cs;
            for (int idim = 1; idim <= NDIM; ++idim) {
                u_wave += std::abs(grid_.uold(id, 1 + idim) / d);
            }
            real_t a_grav = 0.0;
            for (int idim = 1; idim <= NDIM; ++idim) {
                a_grav += std::abs(grid_.f(id, idim));
            }
            real_t alpha = a_grav * dx / (u_wave * u_wave);
            alpha = std::max(alpha, (real_t)0.0001);
            real_t dt_cand = (dx / u_wave) * (std::sqrt(1.0 + 2.0 * courant_factor * alpha) - 1.0) / alpha;
            dt_max = std::min(dt_max, dt_cand);
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
                real_t sum_gamma_p = 0.0;
                if (params::barotropic_eos) {
                    real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
                    real_t temp_mu = T2;
                    if (params::barotropic_eos_form == "polytrope") temp_mu = T2 * std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
                    else if (params::barotropic_eos_form == "double_polytrope") temp_mu = T2 * (1.0 + std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0));
                    p = d * temp_mu / (params::units_velocity * params::units_velocity);
                } else {
                    real_t e_nonthermal = 0.0;
                    for (int ie = 0; ie < nener_; ++ie) {
                        real_t e_rad = grid_.uold(id, iener + 1 + ie);
                        e_nonthermal += e_rad;
                        sum_gamma_p += grid_.gamma_rad[ie] * (grid_.gamma_rad[ie] - 1.0) * e_rad;
                    }
                    p = std::max((grid_.uold(id, iener) - 0.5 * d * v2 - e_nonthermal) * (gamma - 1.0), d * 1e-10);
                    p = std::min(p, d * 1e2);
                }
                
                real_t cs = get_cs(d, p, sum_gamma_p);
                real_t u_wave = NDIM * cs;
                for (int idim = 1; idim <= NDIM; ++idim) {
                    u_wave += std::abs(grid_.uold(id, 1 + idim) / d);
                }
                real_t a_grav = 0.0;
                for (int idim = 1; idim <= NDIM; ++idim) {
                    a_grav += std::abs(grid_.f(id, idim));
                }
                real_t alpha = a_grav * dx / (u_wave * u_wave);
                alpha = std::max(alpha, (real_t)0.0001);
                real_t dt_cand = (dx / u_wave) * (std::sqrt(1.0 + 2.0 * courant_factor * alpha) - 1.0) / alpha;
                /*
                if (dt_cand < 1e-8) {
                    std::cout << "[DEBUG] Timestep collapse! idc=" << id << " max_acc=" << max_acc << " dt_cfl=" << dt_cfl << " dt_grav=" << dt_grav << " v_mag=" << v_mag << " cs=" << cs << std::endl;
                }
                */
                if (id == 1 && ilevel == 0) {
                    // we can check flux later
                }
                dt_max = std::min(dt_max, dt_cand);
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
        
        real_t e_nonthermal = 0.0;
        for (int ie = 0; ie < nener_; ++ie) {
            e_nonthermal += grid_.uold(idc, iener + 1 + ie);
        }
        
        real_t e_int = 0;
        if (params::barotropic_eos) {
            real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
            real_t temp_mu = T2;
            if (params::barotropic_eos_form == "polytrope") temp_mu = T2 * std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
            else if (params::barotropic_eos_form == "double_polytrope") temp_mu = T2 * (1.0 + std::pow(d * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0));
            real_t p = d * temp_mu / (params::units_velocity * params::units_velocity);
            e_int = p / (gam - 1.0);
        } else {
            e_int = std::max(grid_.uold(idc, iener) - 0.5*d*v2_old - e_nonthermal, d * 1e-10 / (gam - 1.0));
        }

        for (int idim = 1; idim <= NDIM; ++idim) {
            grid_.uold(idc, 1 + idim) += d * grid_.f(idc, idim) * dt;
            grid_.uold(idc, 1 + idim) = std::clamp(grid_.uold(idc, 1 + idim), -d * 10.0, d * 10.0);
        }
        
        real_t v2_new = 0;
        for(int idim=1; idim<=NDIM; ++idim) { real_t v = grid_.uold(idc, 1+idim)/d; v2_new += v*v; }
        grid_.uold(idc, iener) = e_int + 0.5*d*v2_new + e_nonthermal;
    };

    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= n2d_val; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            apply_synchro(idc);
        }
        igrid = grid_.next[igrid - 1];
    }
}

} // namespace ramses
