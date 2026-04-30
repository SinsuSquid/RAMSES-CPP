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

real_t MhdSolver::compute_max_div_b(int ilevel, real_t dx) {
    real_t max_div = 0.0;
    int nvar_pure = grid_.nvar - 3;
    
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                
                // Div B = (Bx_r - Bx_l) / dx + (By_r - By_l) / dy + (Bz_r - Bz_l) / dz
                real_t div = (grid_.uold(ind_cell, nvar_pure + 1) - grid_.uold(ind_cell, 6));
                if (NDIM > 1) div += (grid_.uold(ind_cell, nvar_pure + 2) - grid_.uold(ind_cell, 7));
                if (NDIM > 2) div += (grid_.uold(ind_cell, nvar_pure + 3) - grid_.uold(ind_cell, 8));
                
                max_div = std::max(max_div, std::abs(div) / dx);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    return max_div;
}

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
                // Map conservative to primitives for Courant check
                real_t d = std::max(grid_.uold(ind_cell, 1), smallr);
                q[0] = d;
                q[1] = grid_.uold(ind_cell, 5); // Total energy? No, find_speed_fast expects P at index 1
                // Wait, find_speed_fast(q, vel_info, gamma) expects:
                // q[0]=d, q[1]=P, q[2]=u, q[3]=A, q[4]=v, q[5]=B, q[6]=w, q[7]=C
                
                real_t u = grid_.uold(ind_cell, 2) / d;
                real_t v = grid_.uold(ind_cell, 3) / d;
                real_t w = grid_.uold(ind_cell, 4) / d;
                real_t A = grid_.uold(ind_cell, 6);
                real_t B = grid_.uold(ind_cell, 7);
                real_t C = grid_.uold(ind_cell, 8);
                
                // Total energy is index 5 (1-based)
                real_t etot = grid_.uold(ind_cell, 5);
                real_t ekin = 0.5 * d * (u*u + v*v + w*w);
                real_t emag = 0.5 * (A*A + B*B + C*C);
                real_t p = std::max((etot - ekin - emag) * (gamma - 1.0), d * 1e-10);
                
                real_t q_mhd[8] = {d, p, u, A, v, B, w, C};
                real_t vel_fast;
                find_speed_fast(q_mhd, vel_fast, gamma);

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
                            } else {
                                // Fill dummy directions with active cell data for NDIM < 3
                                int active_i2 = i2;
                                int active_j2 = (NDIM > 1) ? j2 : 0;
                                int active_k2 = (NDIM > 2) ? k2 : 0;
                                int active_icell_pos = 1 + active_i2 + 2*active_j2 + 4*active_k2;
                                if (ison > 0 && active_icell_pos <= constants::twotondim) {
                                    int ind_cell = grid_.ncoarse + (active_icell_pos - 1) * grid_.ngridmax + ison;
                                    if (ind_cell > 0 && ind_cell <= grid_.ncell) {
                                        for (int iv = 1; iv <= nvar; ++iv) stencil.uloc[i3][j3][k3][iv - 1] = grid_.uold(ind_cell, iv);
                                    }
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
                       const real_t dq[6][6][6][3][20], const real_t dbf[6][6][6][3][2],
                       real_t dt, real_t dx,
                       real_t qm[6][6][6][3][20], real_t qp[6][6][6][3][20]) {
    real_t dtdx = dt / dx;
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    int nvar_pure = grid_.nvar - 3;
    real_t smallp = 1e-10; // TODO: use parameter

    // 1. Compute edge-centered electric fields (Ez = uBy - vBx)
    real_t ez[6][6][6] = {0.0}, ey[6][6][6] = {0.0}, ex[6][6][6] = {0.0};
    for(int k=1; k<6; ++k) for(int j=1; j<6; ++j) for(int i=1; i<6; ++i) {
        real_t u = 0.25 * (qloc[i-1][j-1][k][1] + qloc[i-1][j][k][1] + qloc[i][j-1][k][1] + qloc[i][j][k][1]);
        real_t v = 0.25 * (qloc[i-1][j-1][k][2] + qloc[i-1][j][k][2] + qloc[i][j-1][k][2] + qloc[i][j][k][2]);
        real_t Bx = 0.5 * (bfloc[i][j-1][k][0][0] + bfloc[i][j][k][0][0]);
        real_t By = 0.5 * (bfloc[i-1][j][k][1][0] + bfloc[i][j][k][1][0]);
        ez[i][j][k] = u * By - v * Bx;

        if (NDIM > 2) {
            real_t w = 0.25 * (qloc[i-1][j][k-1][3] + qloc[i-1][j][k][3] + qloc[i][j][k-1][3] + qloc[i][j][k][3]);
            real_t u3 = 0.25 * (qloc[i-1][j][k-1][1] + qloc[i-1][j][k][1] + qloc[i][j][k-1][1] + qloc[i][j][k][1]);
            real_t Bz = 0.5 * (bfloc[i-1][j][k][2][0] + bfloc[i][j][k][2][0]);
            real_t Bx3 = 0.5 * (bfloc[i][j][k-1][0][0] + bfloc[i][j][k][0][0]);
            ey[i][j][k] = w * Bx3 - u3 * Bz;

            real_t v4 = 0.25 * (qloc[i][j-1][k-1][2] + qloc[i][j-1][k][2] + qloc[i][j][k-1][2] + qloc[i][j][k][2]);
            real_t w4 = 0.25 * (qloc[i][j-1][k-1][3] + qloc[i][j-1][k][3] + qloc[i][j][k-1][3] + qloc[i][j][k][3]);
            real_t By4 = 0.5 * (bfloc[i][j][k-1][1][0] + bfloc[i][j][k][0][0]); // TYPO in By4?
            // Corrected By4: By is index 1.
            By4 = 0.5 * (bfloc[i][j][k-1][1][0] + bfloc[i][j][k][1][0]);
            real_t Bz4 = 0.5 * (bfloc[i][j-1][k][2][0] + bfloc[i][j][k][2][0]);
            ex[i][j][k] = v4 * Bz4 - w4 * By4;
        }
    }

    for (int k = 1; k < 5; ++k) {
        for (int j = 1; j < 5; ++j) {
            for (int i = 1; i < 5; ++i) {
                // Cell centered values
                real_t r = qloc[i][j][k][0];
                real_t u = qloc[i][j][k][1];
                real_t v = qloc[i][j][k][2];
                real_t w = qloc[i][j][k][3];
                real_t p = qloc[i][j][k][4];
                real_t A = qloc[i][j][k][5]; // Bx
                real_t B = qloc[i][j][k][6]; // By
                real_t C = qloc[i][j][k][7]; // Bz

                // TVD slopes
                real_t drx = 0.5 * dq[i][j][k][0][0]; real_t dry = 0.5 * dq[i][j][k][1][0]; real_t drz = 0.5 * dq[i][j][k][2][0];
                real_t dux = 0.5 * dq[i][j][k][0][1]; real_t duy = 0.5 * dq[i][j][k][1][1]; real_t duz = 0.5 * dq[i][j][k][2][1];
                real_t dvx = 0.5 * dq[i][j][k][0][2]; real_t dvy = 0.5 * dq[i][j][k][1][2]; real_t dvz = 0.5 * dq[i][j][k][2][2];
                real_t dwx = 0.5 * dq[i][j][k][0][3]; real_t dwy = 0.5 * dq[i][j][k][1][3]; real_t dwz = 0.5 * dq[i][j][k][2][3];
                real_t dpx = 0.5 * dq[i][j][k][0][4]; real_t dpy = 0.5 * dq[i][j][k][1][4]; real_t dpz = 0.5 * dq[i][j][k][2][4];
                real_t dAx = 0.5 * dq[i][j][k][0][5]; real_t dAy = 0.5 * dq[i][j][k][1][5]; real_t dAz = 0.5 * dq[i][j][k][2][5];
                real_t dBx = 0.5 * dq[i][j][k][0][6]; real_t dBy = 0.5 * dq[i][j][k][1][6]; real_t dBz = 0.5 * dq[i][j][k][2][6];
                real_t dCx = 0.5 * dq[i][j][k][0][7]; real_t dCy = 0.5 * dq[i][j][k][1][7]; real_t dCz = 0.5 * dq[i][j][k][2][7];

                // Face-centered B prediction using EMFs
                real_t BLx = bfloc[i][j][k][0][0], BRx = bfloc[i+1][j][k][0][0];
                real_t BLy = bfloc[i][j][k][1][0], BRy = bfloc[i][j+1][k][1][0];
                real_t BLz = bfloc[i][j][k][2][0], BRz = bfloc[i][j][k+1][2][0];

                if (NDIM > 1) {
                    real_t sBLx = (ez[i][j+1][k] - ez[i][j][k]) * dtdx * 0.5;
                    real_t sBRx = (ez[i+1][j+1][k] - ez[i+1][j][k]) * dtdx * 0.5;
                    real_t sBLy = -(ez[i+1][j][k] - ez[i][j][k]) * dtdx * 0.5;
                    real_t sBRy = -(ez[i+1][j+1][k] - ez[i][j+1][k]) * dtdx * 0.5;
                    if (NDIM > 2) {
                        sBLx += -(ey[i][j][k+1] - ey[i][j][k]) * dtdx * 0.5;
                        sBRx += -(ey[i+1][j][k+1] - ey[i+1][j][k]) * dtdx * 0.5;
                        sBLy += (ex[i][j][k+1] - ex[i][j][k]) * dtdx * 0.5;
                        sBRy += (ex[i][j+1][k+1] - ex[i][j+1][k]) * dtdx * 0.5;
                    }
                    BLx += sBLx; BRx += sBRx; BLy += sBLy; BRy += sBRy;
                }
                if (NDIM > 2) {
                    real_t sBLz = (ey[i+1][j][k] - ey[i][j][k]) * dtdx * 0.5 - (ex[i][j+1][k] - ex[i][j][k]) * dtdx * 0.5;
                    real_t sBRz = (ey[i+1][j][k+1] - ey[i][j][k+1]) * dtdx * 0.5 - (ex[i][j+1][k+1] - ex[i][j][k+1]) * dtdx * 0.5;
                    BLz += sBLz; BRz += sBRz;
                }

                // Cell-centered source terms (Full unsplit prediction)
                real_t sr0 = (-u*drx - r*dux) * dtdx;
                real_t su0 = (-u*dux - (dpx + B*dBx + C*dCx)/r) * dtdx;
                real_t sv0 = (-u*dvx + A*dBx/r) * dtdx;
                real_t sw0 = (-u*dwx + A*dCx/r) * dtdx;
                real_t sp0 = (-u*dpx - gamma*p*dux) * dtdx;
                real_t sAx0 = 0.0, sAy0 = (-u*dBx + A*dvx - B*dux) * dtdx, sAz0 = (-u*dCx + A*dwx - C*dux) * dtdx;

                if (NDIM > 1) {
                    sr0 += (-v*dry - r*dvy) * dtdx;
                    su0 += (-v*duy + B*dAy/r) * dtdx;
                    sv0 += (-v*dvy - (dpy + A*dAy + C*dCy)/r) * dtdx;
                    sw0 += (-v*dwy + B*dCy/r) * dtdx;
                    sp0 += (-v*dpy - gamma*p*dvy) * dtdx;
                    sAx0 += (-v*dAy + B*duy - A*dvy) * dtdx;
                    sAz0 += (-v*dCy + B*dwy - C*dvy) * dtdx;
                }
                if (NDIM > 2) {
                    sr0 += (-w*drz - r*dwz) * dtdx;
                    su0 += (-w*duz + C*dAz/r) * dtdx;
                    sv0 += (-w*dvz + C*dBz/r) * dtdx;
                    sw0 += (-w*dwz - (dpz + A*dAz + B*dBz)/r) * dtdx;
                    sp0 += (-w*dpz - gamma*p*dwz) * dtdx;
                    sAx0 += (-w*dAz + C*duz - A*dwz) * dtdx;
                    sAy0 += (-w*dBz + C*dvz - B*dwz) * dtdx;
                }

                // Apply source terms
                r += sr0; u += su0; v += sv0; w += sw0; p += sp0; A = 0.5*(BLx+BRx); B = 0.5*(BLy+BRy); C = 0.5*(BLz+BRz);
                
                // Construct interface states
                auto set_states = [&](int idim, real_t dr, real_t du, real_t dv, real_t dw, real_t dp, real_t dA, real_t dB, real_t dC, real_t BL_val, real_t BR_val) {
                    qp[i][j][k][idim][0] = std::max(1e-10, r - dr);
                    qp[i][j][k][idim][1] = u - du;
                    qp[i][j][k][idim][2] = v - dv;
                    qp[i][j][k][idim][3] = w - dw;
                    qp[i][j][k][idim][4] = std::max(smallp, p - dp);

                    qm[i][j][k][idim][0] = std::max(1e-10, r + dr);
                    qm[i][j][k][idim][1] = u + du;
                    qm[i][j][k][idim][2] = v + dv;
                    qm[i][j][k][idim][3] = w + dw;
                    qm[i][j][k][idim][4] = std::max(smallp, p + dp);

                    // Map B-fields correctly
                    if (idim == 0) {
                        qp[i][j][k][idim][5] = BL_val; qm[i][j][k][idim][5] = BR_val;
                        qp[i][j][k][idim][6] = B - dB; qm[i][j][k][idim][6] = B + dB;
                        qp[i][j][k][idim][7] = C - dC; qm[i][j][k][idim][7] = C + dC;
                    } else if (idim == 1) {
                        qp[i][j][k][idim][5] = A - dA; qm[i][j][k][idim][5] = A + dA;
                        qp[i][j][k][idim][6] = BL_val; qm[i][j][k][idim][6] = BR_val;
                        qp[i][j][k][idim][7] = C - dC; qm[i][j][k][idim][7] = C + dC;
                    } else {
                        qp[i][j][k][idim][5] = A - dA; qm[i][j][k][idim][5] = A + dA;
                        qp[i][j][k][idim][6] = B - dB; qm[i][j][k][idim][6] = B + dB;
                        qp[i][j][k][idim][7] = BL_val; qm[i][j][k][idim][7] = BR_val;
                    }
                };

                set_states(0, drx, dux, dvx, dwx, dpx, dAx, dBx, dCx, BLx, BRx);
                if (NDIM > 1) set_states(1, dry, duy, dvy, dwy, dpy, dAy, dBy, dCy, BLy, BRy);
                if (NDIM > 2) set_states(2, drz, duz, dvz, dwz, dpz, dAz, dBz, dCz, BLz, BRz);
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
                        const real_t bfloc[6][6][6][3][2],
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
                qL_mhd[3] = bfloc[i][j][k][idim][1]; // Bn (staggered)
                qL_mhd[4] = qm[i][j][k][idim][1 + (idim + 1) % 3];
                qL_mhd[5] = qm[i][j][k][idim][5 + (idim + 1) % 3];
                qL_mhd[6] = qm[i][j][k][idim][1 + (idim + 2) % 3];
                qL_mhd[7] = qm[i][j][k][idim][5 + (idim + 2) % 3];

                qR_mhd[0] = qp[ni][nj][nk][idim][0];
                qR_mhd[1] = qp[ni][nj][nk][idim][4];
                qR_mhd[2] = qp[ni][nj][nk][idim][1 + idim];
                qR_mhd[3] = bfloc[ni][nj][nk][idim][0]; // Bn (staggered)
                qR_mhd[4] = qp[ni][nj][nk][idim][1 + (idim + 1) % 3];
                qR_mhd[5] = qp[ni][nj][nk][idim][5 + (idim + 1) % 3];
                qR_mhd[6] = qp[ni][nj][nk][idim][1 + (idim + 2) % 3];
                qR_mhd[7] = qp[ni][nj][nk][idim][5 + (idim + 2) % 3];

                std::string riemann = config_.get("hydro_params", "riemann", "hlld");
                if (riemann == "llf") {
                    llf(qL_mhd, qR_mhd, f_mhd, gamma);
                } else {
                    hlld(qL_mhd, qR_mhd, f_mhd, gamma);
                }

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

        // 1. Compute slopes for cell-centered variables
        real_t dq[6][6][6][3][20] = {0.0};
        int slope_type = config_.get_int("hydro_params", "slope_type", 1);
        for(int k=1; k<5; ++k) for(int j=1; j<5; ++j) for(int i=1; i<5; ++i) {
            for(int iv=0; iv<nvar_pure; ++iv) {
                dq[i][j][k][0][iv] = SlopeLimiter::compute_slope(qloc[i-1][j][k][iv], qloc[i][j][k][iv], qloc[i+1][j][k][iv], slope_type);
                if (NDIM > 1) dq[i][j][k][1][iv] = SlopeLimiter::compute_slope(qloc[i][j-1][k][iv], qloc[i][j][k][iv], qloc[i][j+1][k][iv], slope_type);
                if (NDIM > 2) dq[i][j][k][2][iv] = SlopeLimiter::compute_slope(qloc[i][j][k-1][iv], qloc[i][j][k][iv], qloc[i][j][k+1][iv], slope_type);
            }
        }

        // 2. Compute slopes for face-centered B-fields
        real_t dbf[6][6][6][3][2] = {0.0};
        for(int k=1; k<5; ++k) for(int j=1; j<5; ++j) for(int i=1; i<5; ++i) {
            // Bx faces are at (i, j, k) and (i+1, j, k)
            for(int i_face=i; i_face<=i+1; ++i_face) {
                if (NDIM > 1) dbf[i_face][j][k][0][0] = SlopeLimiter::compute_slope(bfloc[i_face][j-1][k][0][0], bfloc[i_face][j][k][0][0], bfloc[i_face][j+1][k][0][0], slope_type);
                if (NDIM > 2) dbf[i_face][j][k][0][1] = SlopeLimiter::compute_slope(bfloc[i_face][j][k-1][0][0], bfloc[i_face][j][k][0][0], bfloc[i_face][j][k+1][0][0], slope_type);
            }
            // By faces
            for(int j_face=j; j_face<=j+1; ++j_face) {
                if (NDIM > 1) dbf[i][j_face][k][1][0] = SlopeLimiter::compute_slope(bfloc[i-1][j_face][k][1][0], bfloc[i][j_face][k][1][0], bfloc[i+1][j_face][k][1][0], slope_type);
                if (NDIM > 2) dbf[i][j_face][k][1][1] = SlopeLimiter::compute_slope(bfloc[i][j_face][k-1][1][0], bfloc[i][j_face][k][1][0], bfloc[i][j_face][k+1][1][0], slope_type);
            }
            // Bz faces
            for(int k_face=k; k_face<=k+1; ++k_face) {
                if (NDIM > 1) dbf[i][j][k_face][2][0] = SlopeLimiter::compute_slope(bfloc[i-1][j][k_face][2][0], bfloc[i][j][k_face][2][0], bfloc[i+1][j][k_face][2][0], slope_type);
                if (NDIM > 2) dbf[i][j][k_face][2][1] = SlopeLimiter::compute_slope(bfloc[i][j-1][k_face][2][0], bfloc[i][j][k_face][2][0], bfloc[i][j+1][k_face][2][0], slope_type);
            }
        }

        // 3. Trace states
        real_t qm[6][6][6][3][20], qp[6][6][6][3][20];
        trace(qloc, bfloc, dq, dbf, dt, dx, qm, qp);

        // 4. Compute fluxes
        real_t fluxes[3][6][6][6][20] = {0.0};
        for (int idim = 0; idim < NDIM; ++idim) {
            cmpflxm(qm, qp, bfloc, idim, gamma, fluxes[idim]);
        }

        // 4. Constrained Transport (EMF)
        real_t emfx[6][6][6] = {0.0}, emfy[6][6][6] = {0.0}, emfz[6][6][6] = {0.0};
        for(int k=1; k<5; ++k) for(int j=1; j<5; ++j) for(int i=1; i<5; ++i) {
            // emfz at (i-1/2, j-1/2, k)
            emfz[i][j][k] = 0.25 * (fluxes[0][i-1][j][k][6] + fluxes[0][i-1][j-1][k][6]
                                   - fluxes[1][i][j-1][k][5] - fluxes[1][i-1][j-1][k][5]);

            if (NDIM > 2) {
                // emfy at (i-1/2, j, k-1/2)
                emfy[i][j][k] = 0.25 * (fluxes[2][i][j][k-1][5] + fluxes[2][i-1][j][k-1][5]
                                       - fluxes[0][i-1][j][k][7] - fluxes[0][i-1][j][k-1][7]);
                // emfx at (i, j-1/2, k-1/2)
                emfx[i][j][k] = 0.25 * (fluxes[1][i][j-1][k][7] + fluxes[1][i][j-1][k-1][7]
                                       - fluxes[2][i][j][k-1][6] - fluxes[2][i][j-1][k-1][6]);
            }
        }

        // 5. Update unew
        for (int k2 = 0; k2 < (NDIM > 2 ? 2 : 1); ++k2) {
        for (int j2 = 0; j2 < (NDIM > 1 ? 2 : 1); ++j2) {
        for (int i2 = 0; i2 < 2; ++i2) {
            int i = 2 + i2; int j = 2 + j2; int k = 2 + k2;
            int icell_pos = 1 + i2 + 2*j2 + 4*k2;
            if (icell_pos <= constants::twotondim) {
                int ind_cell = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + igrid;
                
                // Update hydro variables (1-5)
                for (int iv = 1; iv <= 5; ++iv) {
                    for (int idim = 0; idim < NDIM; ++idim) {
                        int ni = i - (idim==0?1:0); int nj = j - (idim==1?1:0); int nk = k - (idim==2?1:0);
                        grid_.unew(ind_cell, iv) += (fluxes[idim][ni][nj][nk][iv-1] - fluxes[idim][i][j][k][iv-1]) * dt_dx;
                    }
                }
                // Update passive scalars (9 to nvar_pure)
                for (int iv = 9; iv <= nvar_pure; ++iv) {
                    for (int idim = 0; idim < NDIM; ++idim) {
                        int ni = i - (idim==0?1:0); int nj = j - (idim==1?1:0); int nk = k - (idim==2?1:0);
                        grid_.unew(ind_cell, iv) += (fluxes[idim][ni][nj][nk][iv-1] - fluxes[idim][i][j][k][iv-1]) * dt_dx;
                    }
                }

                // CT update for B-left faces (unew 6,7,8)
                if (NDIM > 1) {
                    real_t dBx_l = 0.0;
                    if (NDIM > 2) dBx_l += (emfy[i][j][k] - emfy[i][j][k+1]) * dt_dx;
                    dBx_l -= (emfz[i][j][k] - emfz[i][j+1][k]) * dt_dx;
                    grid_.unew(ind_cell, 6) += dBx_l;
                }
                if (NDIM > 1) {
                    real_t dBy_l = (emfz[i][j][k] - emfz[i+1][j][k]) * dt_dx;
                    if (NDIM > 2) dBy_l -= (emfx[i][j][k] - emfx[i][j][k+1]) * dt_dx;
                    grid_.unew(ind_cell, 7) += dBy_l;
                }
                if (NDIM > 2) {
                    real_t dBz_l = (emfx[i][j][k] - emfx[i][j+1][k]) * dt_dx - (emfy[i][j][k] - emfy[i+1][j][k]) * dt_dx;
                    grid_.unew(ind_cell, 8) += dBz_l;
                }

                // CT update for B-right faces (unew nvar_pure + 1, 2, 3)
                if (NDIM > 1) {
                    real_t dBx_r = 0.0;
                    if (NDIM > 2) dBx_r += (emfy[i+1][j][k] - emfy[i+1][j][k+1]) * dt_dx;
                    dBx_r -= (emfz[i+1][j][k] - emfz[i+1][j+1][k]) * dt_dx;
                    grid_.unew(ind_cell, nvar_pure + 1) += dBx_r;
                }
                if (NDIM > 1) {
                    real_t dBy_r = (emfz[i+1][j][k] - emfz[i+1][j+1][k]) * dt_dx;
                    if (NDIM > 2) dBy_r -= (emfx[i+1][j][k] - emfx[i+1][j][k+1]) * dt_dx;
                    grid_.unew(ind_cell, nvar_pure + 2) += dBy_r;
                }
                if (NDIM > 2) {
                    real_t dBz_r = (emfx[i][j][k+1] - emfx[i+1][j][k+1]) * dt_dx - (emfy[i][j+1][k+1] - emfy[i+1][j+1][k+1]) * dt_dx;
                    grid_.unew(ind_cell, nvar_pure + 3) += dBz_r;
                }
            }
        }
        }
        }
    }
}
void MhdSolver::llf(const real_t* ql, const real_t* qr, real_t* fgdnv, real_t gamma) {
    real_t fl[9], fr[9], cl, cr;
    find_mhd_flux(ql, nullptr, fl, gamma);
    find_mhd_flux(qr, nullptr, fr, gamma);
    find_speed_fast(ql, cl, gamma);
    find_speed_fast(qr, cr, gamma);
    
    real_t s_max = std::max(std::abs(ql[2]) + cl, std::abs(qr[2]) + cr);
    
    real_t entho = 1.0 / (gamma - 1.0);
    auto get_cons = [&](const real_t* q, real_t* c) {
        real_t d = q[0], p = q[1], u = q[2], A = q[3], v = q[4], B = q[5], w = q[6], C = q[7];
        c[0] = d; c[1] = d * u; c[2] = d * v; c[3] = d * w;
        real_t ekin = 0.5 * d * (u*u + v*v + w*w);
        real_t emag = 0.5 * (A*A + B*B + C*C);
        c[4] = p * entho + ekin + emag;
        c[5] = A; c[6] = B; c[7] = C;
        c[8] = p * entho;
    };
    
    real_t cll[9], crr[9];
    get_cons(ql, cll); get_cons(qr, crr);
    
    for (int i = 0; i < 9; ++i) {
        fgdnv[i] = 0.5 * (fl[i] + fr[i] - s_max * (crr[i] - cll[i]));
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
