#include "ramses/MhdSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Muscl.hpp"
#include "ramses/SlopeLimiter.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

namespace ramses {

void MhdSolver::set_unew(int ilevel) {
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                for (int ivar = 1; ivar <= grid_.nvar; ++ivar) {
                    grid_.unew(ind_cell, ivar) = grid_.uold(ind_cell, ivar);
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void MhdSolver::set_uold(int ilevel) {
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                for (int ivar = 1; ivar <= grid_.nvar; ++ivar) {
                    grid_.uold(ind_cell, ivar) = grid_.unew(ind_cell, ivar);
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

real_t MhdSolver::compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor) {
    real_t dt_max = 1e30;
    real_t smallr = 1e-10;

    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                real_t q[8];
                for(int i=0; i<8; ++i) q[i] = grid_.uold(ind_cell, i+1);
                
                real_t vel_fast;
                find_speed_fast(q, vel_fast, gamma);

                real_t d = std::max(grid_.uold(ind_cell, 1), smallr);
                real_t u = grid_.uold(ind_cell, 2) / d;
                real_t v = grid_.uold(ind_cell, 3) / d;
                real_t w = grid_.uold(ind_cell, 4) / d;
                real_t vel = std::sqrt(u*u + v*v + w*w);
                
                real_t dt_cell = dx / (vel + vel_fast);
                dt_max = std::min(dt_max, dt_cell);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    return dt_max * courant_factor;
}

void MhdSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    std::vector<int> active_octs;
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            active_octs.push_back(igrid);
            igrid = grid_.next[igrid - 1];
        }
    }
    
    if (!active_octs.empty()) {
        godfine1(active_octs, ilevel, dt, dx);
    }
}

void MhdSolver::gather_stencil(int igrid, int ilevel, LocalStencil& stencil) {
    int nbors_father[27];
    grid_.get_3x3x3_father(igrid, nbors_father);
    int nvar = grid_.nvar;
    
    for (int k1 = 0; k1 < 3; ++k1) {
        for (int j1 = 0; j1 < 3; ++j1) {
            for (int i1 = 0; i1 < 3; ++i1) {
                int ifather = nbors_father[i1 + 3*j1 + 9*k1];
                int ison = (ifather > 0 && ifather <= (int)grid_.son.size() - 1) ? grid_.son[ifather] : 0;
                
                for (int k2 = 0; k2 < 2; ++k2) {
                    for (int j2 = 0; j2 < 2; ++j2) {
                        for (int i2 = 0; i2 < 2; ++i2) {
                            int i3 = i1 * 2 + i2; int j3 = j1 * 2 + j2; int k3 = k1 * 2 + k2;
                            
                            for (int iv = 0; iv < 20; ++iv) stencil.uloc[i3][j3][k3][iv] = 0.0;
                            stencil.refined[i3][j3][k3] = false;

                            int icell_pos = 1 + i2 + 2*j2 + 4*k2;
                            if (ison > 0 && icell_pos <= constants::twotondim) {
                                int ind_cell = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + ison;
                                
                                if (ind_cell > 0 && ind_cell <= grid_.ncell) {
                                    for (int iv = 1; iv <= nvar; ++iv) stencil.uloc[i3][j3][k3][iv - 1] = grid_.uold(ind_cell, iv);
                                    stencil.refined[i3][j3][k3] = (grid_.son[ind_cell] > 0);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void MhdSolver::trace(const real_t qloc[6][6][6][20], const real_t bfloc[6][6][6][3][2],
                       const real_t dq[6][6][6][3][20], real_t dt, real_t dx,
                       real_t qm[6][6][6][3][20], real_t qp[6][6][6][3][20]) {
    real_t dtdx = dt / dx;
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    int nvar_pure = grid_.nvar - 3;

    for (int k = 1; k < 5; ++k) {
        for (int j = 1; j < 5; ++j) {
            for (int i = 1; i < 5; ++i) {
                // Cell centered values
                real_t r = qloc[i][j][k][0];
                real_t u = qloc[i][j][k][1];
                real_t v = qloc[i][j][k][2];
                real_t w = qloc[i][j][k][3];
                real_t p = qloc[i][j][k][4];
                real_t A = qloc[i][j][k][5];
                real_t B = qloc[i][j][k][6];
                real_t C = qloc[i][j][k][7];

                for (int idim = 0; idim < NDIM; ++idim) {
                    real_t dr_d = 0.5 * dq[i][j][k][idim][0];
                    real_t du_d = 0.5 * dq[i][j][k][idim][1];
                    real_t dv_d = 0.5 * dq[i][j][k][idim][2];
                    real_t dw_d = 0.5 * dq[i][j][k][idim][3];
                    real_t dp_d = 0.5 * dq[i][j][k][idim][4];
                    real_t dA_d = 0.5 * dq[i][j][k][idim][5];
                    real_t dB_d = 0.5 * dq[i][j][k][idim][6];
                    real_t dC_d = 0.5 * dq[i][j][k][idim][7];

                    // Map variables to normal/transverse
                    real_t vn, vt1, vt2, Bn, Bt1, Bt2;
                    real_t dvn, dvt1, dvt2, dBn, dBt1, dBt2;

                    if (idim == 0) {
                        vn = u; vt1 = v; vt2 = w; Bn = A; Bt1 = B; Bt2 = C;
                        dvn = du_d; dvt1 = dv_d; dvt2 = dw_d; dBn = dA_d; dBt1 = dB_d; dBt2 = dC_d;
                    } else if (idim == 1) {
                        vn = v; vt1 = w; vt2 = u; Bn = B; Bt1 = C; Bt2 = A;
                        dvn = dv_d; dvt1 = dw_d; dvt2 = du_d; dBn = dB_d; dBt1 = dC_d; dBt2 = dA_d;
                    } else {
                        vn = w; vt1 = u; vt2 = v; Bn = C; Bt1 = A; Bt2 = B;
                        dvn = dw_d; dvt1 = du_d; dvt2 = dv_d; dBn = dC_d; dBt1 = dA_d; dBt2 = dB_d;
                    }

                    // MHD Source terms (1D part of the prediction)
                    real_t sr = -vn * dr_d - r * dvn;
                    real_t sp = -vn * dp_d - gamma * p * dvn;
                    real_t svn = -vn * dvn - (dp_d + Bt1 * dBt1 + Bt2 * dBt2) / r;
                    real_t svt1 = -vn * dvt1 + (Bn * dBt1) / r;
                    real_t svt2 = -vn * dvt2 + (Bn * dBt2) / r;
                    real_t sBt1 = -vn * dBt1 + Bn * dvt1 - Bt1 * dvn;
                    real_t sBt2 = -vn * dBt2 + Bn * dvt2 - Bt2 * dvn;
                    // dBn/dt = 0 in 1D source term if using CT

                    // Update states with time prediction
                    real_t r_pred = r + sr * dtdx;
                    real_t p_pred = p + sp * dtdx;
                    real_t vn_pred = vn + svn * dtdx;
                    real_t vt1_pred = vt1 + svt1 * dtdx;
                    real_t vt2_pred = vt2 + svt2 * dtdx;
                    real_t Bt1_pred = Bt1 + sBt1 * dtdx;
                    real_t Bt2_pred = Bt2 + sBt2 * dtdx;

                    // Map back to global qm/qp
                    auto set_qmqp = [&](int iv, real_t val, real_t slope) {
                        qm[i][j][k][idim][iv] = val + slope;
                        qp[i][j][k][idim][iv] = val - slope;
                    };

                    set_qmqp(0, r_pred, dr_d);
                    set_qmqp(4, p_pred, dp_d);
                    
                    if (idim == 0) {
                        set_qmqp(1, vn_pred, du_d); set_qmqp(2, vt1_pred, dv_d); set_qmqp(3, vt2_pred, dw_d);
                        set_qmqp(5, Bn, dA_d); set_qmqp(6, Bt1_pred, dB_d); set_qmqp(7, Bt2_pred, dC_d);
                    } else if (idim == 1) {
                        set_qmqp(2, vn_pred, dv_d); set_qmqp(3, vt1_pred, dw_d); set_qmqp(1, vt2_pred, du_d);
                        set_qmqp(6, Bn, dB_d); set_qmqp(7, Bt1_pred, dC_d); set_qmqp(5, Bt2_pred, dA_d);
                    } else {
                        set_qmqp(3, vn_pred, dw_d); set_qmqp(1, vt1_pred, du_d); set_qmqp(2, vt2_pred, dv_d);
                        set_qmqp(7, Bn, dC_d); set_qmqp(5, Bt1_pred, dA_d); set_qmqp(6, Bt2_pred, dB_d);
                    }
                }
            }
        }
    }
}

void MhdSolver::ctoprim(const real_t u[20], real_t q[20], real_t bf[3][2], real_t gamma) {
    const real_t smallr = 1e-10;
    real_t d = std::max(u[0], smallr);
    q[0] = d;
    q[1] = u[1] / d; // vx
    q[2] = u[2] / d; // vy
    q[3] = u[3] / d; // vz
    
    // Face-centered B
    bf[0][0] = u[5]; // Bx left
    bf[1][0] = u[6]; // By left
    bf[2][0] = u[7]; // Bz left
    
    int nvar_pure = grid_.nvar - 3; // B-right are extra
    bf[0][1] = u[nvar_pure + 0]; // Bx right
    bf[1][1] = u[nvar_pure + 1]; // By right
    bf[2][1] = u[nvar_pure + 2]; // Bz right
    
    // Cell-centered B
    q[5] = 0.5 * (bf[0][0] + bf[0][1]);
    q[6] = 0.5 * (bf[1][0] + bf[1][1]);
    q[7] = 0.5 * (bf[2][0] + bf[2][1]);
    
    real_t e_kin = 0.5 * d * (q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    real_t e_mag = 0.5 * (q[5]*q[5] + q[6]*q[6] + q[7]*q[7]);
    
    // Total energy is u[4]
    real_t e_int = u[4] - e_kin - e_mag;
    q[4] = std::max(e_int * (gamma - 1.0), d * 1e-10); // Pressure
    
    // Passive scalars
    for (int iv = 8; iv < nvar_pure; ++iv) {
        q[iv] = u[iv] / d;
    }
}

void MhdSolver::cmpflxm(const real_t qm[6][6][6][3][20], const real_t qp[6][6][6][3][20],
                        int idim, real_t gamma, real_t flux[6][6][6][20]) {
    int nvar_pure = grid_.nvar - 3;
    for (int k = 1; k < 5; ++k) {
        for (int j = 1; j < 5; ++j) {
            for (int i = 1; i < 5; ++i) {
                int ni = i + (idim == 0 ? 1 : 0);
                int nj = j + (idim == 1 ? 1 : 0);
                int nk = k + (idim == 2 ? 1 : 0);
                if (ni >= 6 || nj >= 6 || nk >= 6) continue;

                real_t qL_mhd[8], qR_mhd[8], f_mhd[9];
                // q: rho, p, vn, Bn, vt1, Bt1, vt2, Bt2
                qL_mhd[0] = qm[i][j][k][idim][0];
                qL_mhd[1] = qm[i][j][k][idim][4];
                qL_mhd[2] = qm[i][j][k][idim][1 + idim];
                qL_mhd[3] = qm[i][j][k][idim][5 + idim]; // Bn
                qL_mhd[4] = qm[i][j][k][idim][1 + (idim + 1) % 3];
                qL_mhd[5] = qm[i][j][k][idim][5 + (idim + 1) % 3];
                qL_mhd[6] = qm[i][j][k][idim][1 + (idim + 2) % 3];
                qL_mhd[7] = qm[i][j][k][idim][5 + (idim + 2) % 3];

                qR_mhd[0] = qp[ni][nj][nk][idim][0];
                qR_mhd[1] = qp[ni][nj][nk][idim][4];
                qR_mhd[2] = qp[ni][nj][nk][idim][1 + idim];
                qR_mhd[3] = qp[ni][nj][nk][idim][5 + idim];
                qR_mhd[4] = qp[ni][nj][nk][idim][1 + (idim + 1) % 3];
                qR_mhd[5] = qp[ni][nj][nk][idim][5 + (idim + 1) % 3];
                qR_mhd[6] = qp[ni][nj][nk][idim][1 + (idim + 2) % 3];
                qR_mhd[7] = qp[ni][nj][nk][idim][5 + (idim + 2) % 3];

                hlld(qL_mhd, qR_mhd, f_mhd, gamma);

                flux[i][j][k][0] = f_mhd[0];
                flux[i][j][k][4] = f_mhd[1];
                flux[i][j][k][1 + idim] = f_mhd[2];
                flux[i][j][k][1 + (idim + 1) % 3] = f_mhd[4];
                flux[i][j][k][5 + (idim + 1) % 3] = f_mhd[5];
                flux[i][j][k][1 + (idim + 2) % 3] = f_mhd[6];
                flux[i][j][k][5 + (idim + 2) % 3] = f_mhd[7];
                
                // Passive scalars
                for(int iv=8; iv<nvar_pure; ++iv) {
                    flux[i][j][k][iv] = (f_mhd[0] > 0) ? f_mhd[0] * qm[i][j][k][idim][iv] : f_mhd[0] * qp[ni][nj][nk][idim][iv];
                }
            }
        }
    }
}

void MhdSolver::godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx) {
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    real_t dt_dx = dt / dx;
    int nvar = grid_.nvar;
    int nvar_pure = nvar - 3;

    for (int igrid : ind_grid) {
        gather_stencil(igrid, ilevel, *stencil_ptr_);
        auto& stencil = *stencil_ptr_;

        real_t qloc[6][6][6][20];
        real_t bfloc[6][6][6][3][2];
        for(int k=0; k<6; ++k) for(int j=0; j<6; ++j) for(int i=0; i<6; ++i) 
            ctoprim(stencil.uloc[i][j][k], qloc[i][j][k], bfloc[i][j][k], gamma);

        // 1. Compute slopes
        real_t dq[6][6][6][3][20] = {0.0};
        for(int k=1; k<5; ++k) for(int j=1; j<5; ++j) for(int i=1; i<5; ++i) {
            for(int iv=0; iv<nvar_pure; ++iv) {
                dq[i][j][k][0][iv] = SlopeLimiter::compute_slope(qloc[i-1][j][k][iv], qloc[i][j][k][iv], qloc[i+1][j][k][iv], 1);
                if (NDIM > 1) dq[i][j][k][1][iv] = SlopeLimiter::compute_slope(qloc[i][j-1][k][iv], qloc[i][j][k][iv], qloc[i][j+1][k][iv], 1);
                if (NDIM > 2) dq[i][j][k][2][iv] = SlopeLimiter::compute_slope(qloc[i][j][k-1][iv], qloc[i][j][k][iv], qloc[i][j][k+1][iv], 1);
            }
        }

        // 2. Trace states
        real_t qm[6][6][6][3][20], qp[6][6][6][3][20];
        trace(qloc, bfloc, dq, dt, dx, qm, qp);

        // 3. Compute fluxes
        real_t fluxes[3][6][6][6][20] = {0.0};
        for (int idim = 0; idim < NDIM; ++idim) {
            cmpflxm(qm, qp, idim, gamma, fluxes[idim]);
        }

        // 4. Constrained Transport (EMF)
        real_t emfx[6][6][6] = {0.0}, emfy[6][6][6] = {0.0}, emfz[6][6][6] = {0.0};
        for(int k=1; k<5; ++k) for(int j=1; j<5; ++j) for(int i=1; i<5; ++i) {
            // Ez = u*By - v*Bx
            // This is a very simplified arithmetic average for EMF
            // In a real code, we use Upwind or a more complex average
            auto get_ez = [&](int ii, int jj, int kk) {
                real_t u_avg = 0.25*(qloc[ii-1][jj-1][kk][1] + qloc[ii-1][jj][kk][1] + qloc[ii][jj-1][kk][1] + qloc[ii][jj][kk][1]);
                real_t v_avg = 0.25*(qloc[ii-1][jj-1][kk][2] + qloc[ii-1][jj][kk][2] + qloc[ii][jj-1][kk][2] + qloc[ii][jj][kk][2]);
                real_t Ax_avg = 0.5*(bfloc[ii][jj-1][kk][0][0] + bfloc[ii][jj][kk][0][0]);
                real_t Ay_avg = 0.5*(bfloc[ii-1][jj][kk][1][0] + bfloc[ii][jj][kk][1][0]);
                return u_avg * Ay_avg - v_avg * Ax_avg;
            };
            emfz[i][j][k] = get_ez(i, j, k);

            if (NDIM > 2) {
                // Ey = w*Bx - u*Bz
                auto get_ey = [&](int ii, int jj, int kk) {
                    real_t u_avg = 0.25*(qloc[ii-1][jj][kk-1][1] + qloc[ii-1][jj][kk][1] + qloc[ii][jj][kk-1][1] + qloc[ii][jj][kk][1]);
                    real_t w_avg = 0.25*(qloc[ii-1][jj][kk-1][3] + qloc[ii-1][jj][kk][3] + qloc[ii][jj][kk-1][3] + qloc[ii][jj][kk][3]);
                    real_t Ax_avg = 0.5*(bfloc[ii][jj][kk-1][0][0] + bfloc[ii][jj][kk][0][0]);
                    real_t Az_avg = 0.5*(bfloc[ii-1][jj][kk][2][0] + bfloc[ii][jj][kk][2][0]);
                    return w_avg * Ax_avg - u_avg * Az_avg;
                };
                emfy[i][j][k] = get_ey(i, j, k);

                // Ex = v*Bz - w*By
                auto get_ex = [&](int ii, int jj, int kk) {
                    real_t v_avg = 0.25*(qloc[ii][jj-1][kk-1][2] + qloc[ii][jj-1][kk][2] + qloc[ii][jj][kk-1][2] + qloc[ii][jj][kk][2]);
                    real_t w_avg = 0.25*(qloc[ii][jj-1][kk-1][3] + qloc[ii][jj-1][kk][3] + qloc[ii][jj][kk-1][3] + qloc[ii][jj][kk][3]);
                    real_t Ay_avg = 0.5*(bfloc[ii][jj][kk-1][1][0] + bfloc[ii][jj][kk][0][0]);
                    real_t Az_avg = 0.5*(bfloc[ii][jj-1][kk][2][0] + bfloc[ii][jj][kk][2][0]);
                    return v_avg * Az_avg - w_avg * Ay_avg;
                };
                emfx[i][j][k] = get_ex(i, j, k);
            }
        }

        // 4. Update unew
        for (int k2 = 0; k2 < (NDIM > 2 ? 2 : 1); ++k2) {
        for (int j2 = 0; j2 < (NDIM > 1 ? 2 : 1); ++j2) {
        for (int i2 = 0; i2 < 2; ++i2) {
            int i = 2 + i2; int j = 2 + j2; int k = 2 + k2;
            int icell_pos = 1 + i2 + 2*j2 + 4*k2;
            if (icell_pos <= constants::twotondim) {
                int ind_cell = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + igrid;
                
                for (int iv = 1; iv <= nvar_pure; ++iv) {
                    for (int idim = 0; idim < NDIM; ++idim) {
                        int ni = i - (idim==0?1:0); int nj = j - (idim==1?1:0); int nk = k - (idim==2?1:0);
                        grid_.unew(ind_cell, iv) += (fluxes[idim][ni][nj][nk][iv-1] - fluxes[idim][i][j][k][iv-1]) * dt_dx;
                    }
                }

                // CT update for B-left faces (unew 6,7,8)
                // dBx = dt/dz (Ey(k) - Ey(k+1)) - dt/dy (Ez(j) - Ez(j+1))
                real_t dBx_l = (emfy[i][j][k] - emfy[i][j][k+1]) * dt_dx - (emfz[i][j][k] - emfz[i][j+1][k]) * dt_dx;
                grid_.unew(ind_cell, 6) += dBx_l;
                if (NDIM > 1) {
                    real_t dBy_l = (emfz[i][j][k] - emfz[i+1][j][k]) * dt_dx - (emfx[i][j][k] - emfx[i][j][k+1]) * dt_dx;
                    grid_.unew(ind_cell, 7) += dBy_l;
                }
                if (NDIM > 2) {
                    real_t dBz_l = (emfx[i][j][k] - emfx[i][j+1][k]) * dt_dx - (emfy[i][j][k] - emfy[i+1][j][k]) * dt_dx;
                    grid_.unew(ind_cell, 8) += dBz_l;
                }

                // CT update for B-right faces (unew nvar-2, nvar-1, nvar)
                real_t dBx_r = (emfy[i+1][j][k] - emfy[i+1][j][k+1]) * dt_dx - (emfz[i+1][j][k] - emfz[i+1][j+1][k]) * dt_dx;
                grid_.unew(ind_cell, nvar_pure + 1) += dBx_r;
                if (NDIM > 1) {
                    real_t dBy_r = (emfz[i][j+1][k] - emfz[i+1][j+1][k]) * dt_dx - (emfx[i][j+1][k] - emfx[i][j+1][k+1]) * dt_dx;
                    grid_.unew(ind_cell, nvar_pure + 2) += dBy_r;
                }
                if (NDIM > 2) {
                    real_t dBz_r = (emfx[i][j][k+1] - emfx[i][j+1][k+1]) * dt_dx - (emfy[i][j][k+1] - emfy[i+1][j][k+1]) * dt_dx;
                    grid_.unew(ind_cell, nvar_pure + 3) += dBz_r;
                }
            }
        }
        }
        }
    }
}
void MhdSolver::hlld(const real_t* qleft, const real_t* qright, real_t* fgdnv, real_t gamma) {
    // Ported HLLD implementation...
    real_t entho = 1.0 / (gamma - 1.0);
    real_t zero = 0.0;
    real_t half = 0.5;
    real_t one = 1.0;

    real_t A = half * (qleft[3] + qright[3]);
    real_t sgnm = (A >= 0) ? 1.0 : -1.0;

    real_t ql[8], qr[8];
    for(int i=0; i<8; ++i) { ql[i]=qleft[i]; qr[i]=qright[i]; }
    ql[3]=A; qr[3]=A;

    real_t rl = ql[0], Pl = ql[1], ul = ql[2];
    real_t vl = ql[4], Bl = ql[5], wl = ql[6], Cl = ql[7];
    real_t ecinl = half * (ul * ul + vl * vl + wl * wl) * rl;
    real_t emagl = half * (A * A + Bl * Bl + Cl * Cl);
    real_t etotl = Pl * entho + ecinl + emagl;
    real_t Ptotl = Pl + emagl;
    real_t vdotBl = ul * A + vl * Bl + wl * Cl;
    real_t eintl = Pl * entho;

    real_t rr = qr[0], Pr = qr[1], ur = qr[2];
    real_t vr = qr[4], Br = qr[5], wr = qr[6], Cr = qr[7];
    real_t ecinr = half * (ur * ur + vr * vr + wr * wr) * rr;
    real_t emagr = half * (A * A + Br * Br + Cr * Cr);
    real_t etotr = Pr * entho + ecinr + emagr;
    real_t Ptotr = Pr + emagr;
    real_t vdotBr = ur * A + vr * Br + wr * Cr;
    real_t eintr = Pr * entho;

    real_t cfastl, cfastr;
    find_speed_fast(ql, cfastl, gamma);
    find_speed_fast(qr, cfastr, gamma);

    real_t SL = std::min(ul, ur) - std::max(cfastl, cfastr);
    real_t SR = std::max(ul, ur) + std::max(cfastl, cfastr);

    real_t rcl = rl * (ul - SL);
    real_t rcr = rr * (SR - ur);
    real_t ustar = (rcr * ur + rcl * ul + (Ptotl - Ptotr)) / (rcr + rcl);
    real_t Ptotstar = (rcr * Ptotl + rcl * Ptotr + rcl * rcr * (ul - ur)) / (rcr + rcl);

    real_t rstarl = rl * (SL - ul) / (SL - ustar);
    real_t estar_l = rl * (SL - ul) * (SL - ustar) - A * A;
    real_t el = rl * (SL - ul) * (SL - ul) - A * A;
    real_t vstarl, Bstarl, wstarl, Cstarl;
    if (std::abs(estar_l) < 1e-6 * A * A) {
        vstarl = vl; Bstarl = Bl; wstarl = wl; Cstarl = Cl;
    } else {
        vstarl = vl - A * Bl * (ustar - ul) / estar_l;
        Bstarl = Bl * el / estar_l;
        wstarl = wl - A * Cl * (ustar - ul) / estar_l;
        Cstarl = Cl * el / estar_l;
    }
    real_t vdotBstarl = ustar * A + vstarl * Bstarl + wstarl * Cstarl;
    real_t etotstarl = ((SL - ul) * etotl - Ptotl * ul + Ptotstar * ustar + A * (vdotBl - vdotBstarl)) / (SL - ustar);
    real_t sqrrstarl = std::sqrt(rstarl);
    real_t calfvenl = std::abs(A) / sqrrstarl;
    real_t SAL = ustar - calfvenl;

    real_t rstarr = rr * (SR - ur) / (SR - ustar);
    real_t estar_r = rr * (SR - ur) * (SR - ustar) - A * A;
    real_t er = rr * (SR - ur) * (SR - ur) - A * A;
    real_t vstarr, Bstarr, wstarr, Cstarr;
    if (std::abs(estar_r) < 1e-6 * A * A) {
        vstarr = vr; Bstarr = Br; wstarr = wr; Cstarr = Cr;
    } else {
        vstarr = vr - A * Br * (ustar - ur) / estar_r;
        Bstarr = Br * er / estar_r;
        wstarr = wr - A * Cr * (ustar - ur) / estar_r;
        Cstarr = Cr * er / estar_r;
    }
    real_t vdotBstarr = ustar * A + vstarr * Bstarr + wstarr * Cstarr;
    real_t etotstarr = ((SR - ur) * etotr - Ptotr * ur + Ptotstar * ustar + A * (vdotBr - vdotBstarr)) / (SR - ustar);
    real_t sqrrstarr = std::sqrt(rstarr);
    real_t calfvenr = std::abs(A) / sqrrstarr;
    real_t SAR = ustar + calfvenr;

    real_t vstarstar = (sqrrstarl * vstarl + sqrrstarr * vstarr + sgnm * (Bstarr - Bstarl)) / (sqrrstarl + sqrrstarr);
    real_t wstarstar = (sqrrstarl * wstarl + sqrrstarr * wstarr + sgnm * (Cstarr - Cstarl)) / (sqrrstarl + sqrrstarr);
    real_t Bstarstar = (sqrrstarl * Bstarr + sqrrstarr * Bstarl + sgnm * sqrrstarl * sqrrstarr * (vstarr - vstarl)) / (sqrrstarl + sqrrstarr);
    real_t Cstarstar = (sqrrstarl * Cstarr + sqrrstarr * Cstarl + sgnm * sqrrstarl * sqrrstarr * (wstarr - wstarl)) / (sqrrstarl + sqrrstarr);
    real_t vdotBstarstar = ustar * A + vstarstar * Bstarstar + wstarstar * Cstarstar;
    real_t etotstarstarl = etotstarl - sgnm * sqrrstarl * (vdotBstarl - vdotBstarstar);
    real_t etotstarstarr = etotstarr + sgnm * sqrrstarr * (vdotBstarr - vdotBstarstar);

    real_t ro, uo, vo, wo, Bo, Co, Ptoto, etoto, vdotBo, einto;
    if (SL > 0.0) {
        ro = rl; uo = ul; vo = vl; wo = wl; Bo = Bl; Co = Cl; Ptoto = Ptotl; etoto = etotl; vdotBo = vdotBl; einto = eintl;
    } else if (SAL > 0.0) {
        ro = rstarl; uo = ustar; vo = vstarl; wo = wstarl; Bo = Bstarl; Co = Cstarl; Ptoto = Ptotstar; etoto = etotstarstarl; vdotBo = vdotBstarl; einto = eintl*(SL-ul)/(SL-ustar);
    } else if (ustar > 0.0) {
        ro = rstarl; uo = ustar; vo = vstarstar; wo = wstarstar; Bo = Bstarstar; Co = Cstarstar; Ptoto = Ptotstar; etoto = etotstarstarl; vdotBo = vdotBstarstar; einto = eintl*(SL-ul)/(SL-ustar);
    } else if (SAR > 0.0) {
        ro = rstarr; uo = ustar; vo = vstarstar; wo = wstarstar; Bo = Bstarstar; Co = Cstarstar; Ptoto = Ptotstar; etoto = etotstarstarr; vdotBo = vdotBstarstar; einto = eintr*(SR-ur)/(SR-ustar);
    } else if (SR > 0.0) {
        ro = rstarr; uo = ustar; vo = vstarr; wo = wstarr; Bo = Bstarr; Co = Cstarr; Ptoto = Ptotstar; etoto = etotstarr; vdotBo = vdotBstarr; einto = eintr*(SR-ur)/(SR-ustar);
    } else {
        ro = rr; uo = ur; vo = vr; wo = wr; Bo = Br; Co = Cr; Ptoto = Ptotr; etoto = etotr; vdotBo = vdotBr; einto = eintr;
    }

    fgdnv[0] = ro * uo;
    fgdnv[1] = (etoto + Ptoto) * uo - A * vdotBo;
    fgdnv[2] = ro * uo * uo + Ptoto - A * A;
    fgdnv[3] = 0.0;
    fgdnv[4] = ro * uo * vo - A * Bo;
    fgdnv[5] = Bo * uo - A * vo;
    fgdnv[6] = ro * uo * wo - A * Co;
    fgdnv[7] = Co * uo - A * wo;
    fgdnv[8] = uo * einto;
}

void MhdSolver::find_mhd_flux(const real_t* qvar, real_t* cvar, real_t* ff, real_t gamma) {
    real_t entho = 1.0 / (gamma - 1.0);
    real_t d = qvar[0], P = qvar[1], u = qvar[2], A = qvar[3];
    real_t v = qvar[4], B = qvar[5], w = qvar[6], C = qvar[7];
    real_t ecin = 0.5 * (u * u + v * v + w * w) * d;
    real_t emag = 0.5 * (A * A + B * B + C * C);
    real_t etot = P * entho + ecin + emag;
    real_t Ptot = P + emag;

    ff[0] = d * u;
    ff[1] = (etot + Ptot) * u - A * (A * u + B * v + C * w);
    ff[2] = d * u * u + Ptot - A * A;
    ff[3] = 0.0;
    ff[4] = d * u * v - A * B;
    ff[5] = B * u - A * v;
    ff[6] = d * u * w - A * C;
    ff[7] = C * u - A * w;
    ff[8] = P * entho * u;
}

void MhdSolver::find_speed_fast(const real_t* qvar, real_t& vel_info, real_t gamma) {
    real_t d = qvar[0], P = qvar[1], u = qvar[2], A = qvar[3];
    real_t v = qvar[4], B = qvar[5], w = qvar[6], C = qvar[7];
    real_t B2 = A * A + B * B + C * C;
    real_t c2 = gamma * P / d;
    real_t d2 = 0.5 * (B2 / d + c2);
    vel_info = std::sqrt(d2 + std::sqrt(std::max(0.0, d2 * d2 - c2 * A * A / d)));
}

} // namespace ramses
