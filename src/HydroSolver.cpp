#include "ramses/HydroSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/RiemannSolver.hpp"
#include "ramses/SlopeLimiter.hpp"
#include "ramses/Muscl.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

namespace ramses {

void HydroSolver::set_unew(int ilevel) {
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ind = 1; ind <= constants::twotondim; ++ind) {
                int ind_cell = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;
                for (int ivar = 1; ivar <= grid_.nvar; ++ivar) {
                    grid_.unew(ind_cell, ivar) = grid_.uold(ind_cell, ivar);
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void HydroSolver::set_uold(int ilevel) {
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ind = 1; ind <= constants::twotondim; ++ind) {
                int ind_cell = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;
                for (int ivar = 1; ivar <= grid_.nvar; ++ivar) {
                    grid_.uold(ind_cell, ivar) = grid_.unew(ind_cell, ivar);
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void HydroSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
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

void HydroSolver::ctoprim(const real_t u[], real_t q[], real_t gamma) {
    const real_t smallr = 1e-10;
    real_t d = std::max(u[0], smallr);
    q[0] = d;
    real_t vel2 = 0.0;
    for (int idim = 1; idim <= NDIM; ++idim) {
        q[idim] = u[idim] / d;
        vel2 += q[idim] * q[idim];
    }
    real_t e_kin = 0.5 * d * vel2;
    real_t e_int = u[NDIM + 1] - e_kin;
    q[NDIM + 1] = std::max(e_int * (gamma - 1.0), d * 1e-10);
    
    // Extra variables
    for (int ivar = NDIM + 2; ivar < grid_.nvar; ++ivar) {
        q[ivar] = u[ivar];
    }
}

void HydroSolver::interpol_hydro(const real_t u1[7][20], real_t u2[8][20]) {
    int nvar = grid_.nvar;
    real_t w[3][20] = {0.0};
    for (int idim = 0; idim < NDIM; ++idim) {
        for (int iv = 0; iv < nvar; ++iv) {
            real_t dlft = 0.5 * (u1[0][iv] - u1[2 * idim + 1][iv]);
            real_t drgt = 0.5 * (u1[2 * idim + 2][iv] - u1[0][iv]);
            if (dlft * drgt <= 0.0) w[idim][iv] = 0.0;
            else w[idim][iv] = (std::abs(dlft) < std::abs(drgt)) ? dlft : drgt;
        }
    }
    for (int ind = 0; ind < constants::twotondim; ++ind) {
        int ix = (ind & 1); int iy = (ind & 2) >> 1; int iz = (ind & 4) >> 2;
        real_t xc[3] = {static_cast<real_t>(ix) - 0.5f, static_cast<real_t>(iy) - 0.5f, static_cast<real_t>(iz) - 0.5f};
        for (int iv = 0; iv < nvar; ++iv) {
            u2[ind][iv] = u1[0][iv];
            for (int idim = 0; idim < NDIM; ++idim) u2[ind][iv] += 2.0 * w[idim][iv] * xc[idim];
        }
    }
}

void HydroSolver::gather_stencil(int igrid, int ilevel, LocalStencil& stencil) {
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
                            
                            for (int iv = 0; iv < nvar; ++iv) stencil.uloc[i3][j3][k3][iv] = 0.0;
                            stencil.refined[i3][j3][k3] = false;

                            int icell_pos = 1 + i2 + 2*j2 + 4*k2;
                            if (ison > 0 && icell_pos <= constants::twotondim) {
                                int ind_cell = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + ison;
                                
                                if (ind_cell > 0 && ind_cell <= grid_.ncell) {
                                    for (int iv = 1; iv <= nvar; ++iv) stencil.uloc[i3][j3][k3][iv - 1] = grid_.uold(ind_cell, iv);
                                    stencil.refined[i3][j3][k3] = (grid_.son[ind_cell] > 0);
                                }
                            } else if (ifather > 0 && ifather <= grid_.ncell) {
                                // Neighbor is not refined, use father's data
                                for (int iv = 1; iv <= nvar; ++iv) stencil.uloc[i3][j3][k3][iv - 1] = grid_.uold(ifather, iv);
                                stencil.refined[i3][j3][k3] = false;
                            }
                        }
                    }
                }
            }
        }
    }
}

real_t HydroSolver::compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor) {
    real_t dt_max = 1e30; const real_t smallr = 1e-10;
    
    // We loop over all cells. For a rigorous Courant condition on AMR, 
    // the max stable dt for the root level must consider the finest level's constraints.
    for (int i = 1; i <= grid_.ncell; ++i) {
        if (grid_.son[i] != 0) continue; // Only compute dt on leaf cells
        
        // Find which level this cell belongs to
        int cell_level = 1;
        if (i > grid_.ncoarse) {
            int igrid = (i - grid_.ncoarse - 1) % grid_.ngridmax + 1;
            int ifather = grid_.father[igrid - 1];
            cell_level = 2;
            while (ifather > grid_.ncoarse) {
                cell_level++;
                int igrid_father = (ifather - grid_.ncoarse - 1) % grid_.ngridmax + 1;
                ifather = grid_.father[igrid_father - 1];
            }
        }
        
        real_t cell_dx = params::boxlen / static_cast<real_t>(params::nx * (1 << (cell_level - 1)));

        real_t d = std::max(grid_.uold(i, 1), smallr);
        real_t vel2 = 0.0;
        real_t vel_max = 0.0;
        real_t g_max = 0.0;
        for (int idim = 1; idim <= NDIM; ++idim) {
            real_t v = grid_.uold(i, 1 + idim) / d;
            vel2 += v * v;
            vel_max = std::max(vel_max, std::abs(v));
            g_max = std::max(g_max, std::abs(grid_.f(i, idim)));
        }
        real_t e_kin = 0.5 * d * vel2;
        real_t e_int = grid_.uold(i, NDIM + 2) - e_kin;
        real_t p = std::max(e_int * (gamma - 1.0), d * 1e-10);
        real_t cs = std::sqrt(gamma * p / d);
        
        // dt based on hydro
        real_t dt_hydro = courant_factor * cell_dx / (vel_max + cs + 1e-20);
        // dt based on gravity
        real_t dt_grav = 1e30;
        if (g_max > 0) dt_grav = std::sqrt(cell_dx / g_max);
        
        real_t dt_cell = std::min(dt_hydro, dt_grav);
        
        // Scale up the cell's dt to the root level dt
        real_t dt_root = dt_cell * (1 << (cell_level - 1));
        dt_max = std::min(dt_max, dt_root);
    }
    return dt_max;
}

void HydroSolver::godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx) {
    real_t gamma = 1.4; 
    real_t dt_dx = dt / dx;
    int nvar = grid_.nvar;

    for (int igrid : ind_grid) {
        gather_stencil(igrid, ilevel, *stencil_ptr_);
        auto& stencil = *stencil_ptr_;

        real_t qloc[6][6][6][20];
        for(int k=0; k<6; ++k) for(int j=0; j<6; ++j) for(int i=0; i<6; ++i) 
            ctoprim(stencil.uloc[i][j][k], qloc[i][j][k], gamma);

        real_t dq[6][6][6][3][20];
        for(int k=1; k<5; ++k) for(int j=1; j<5; ++j) for(int i=1; i<5; ++i) {
            for(int iv=0; iv<nvar; ++iv) {
                dq[i][j][k][0][iv] = SlopeLimiter::compute_slope(qloc[i-1][j][k][iv], qloc[i][j][k][iv], qloc[i+1][j][k][iv], 1);
                if (NDIM > 1) dq[i][j][k][1][iv] = SlopeLimiter::compute_slope(qloc[i][j-1][k][iv], qloc[i][j][k][iv], qloc[i][j+1][k][iv], 1);
                if (NDIM > 2) dq[i][j][k][2][iv] = SlopeLimiter::compute_slope(qloc[i][j][k-1][iv], qloc[i][j][k][iv], qloc[i][j][k+1][iv], 1);
            }
        }

        real_t flux[6][6][6][3][20] = {0.0};
        for (int idim = 0; idim < NDIM; ++idim) {
            for(int k=1; k<5; ++k) for(int j=1; j<5; ++j) for(int i=1; i<5; ++i) {
                real_t ql_interface[20], qr_interface[20], s0[20] = {0.0};
                real_t qm_tmp[20], qp_tmp[20];
                
                Muscl::predict(qloc[i][j][k], dq[i][j][k][idim], s0, dt_dx, qm_tmp, qp_tmp, nvar);
                ql_interface[0] = qm_tmp[0];
                ql_interface[1] = qm_tmp[1+idim];
                if (NDIM > 1) ql_interface[2] = qm_tmp[1+(idim+1)%NDIM];
                if (NDIM > 2) ql_interface[3] = qm_tmp[1+(idim+2)%NDIM];
                ql_interface[NDIM + 1] = qm_tmp[NDIM + 1];

                int ni = i + (idim==0?1:0); int nj = j + (idim==1?1:0); int nk = k + (idim==2?1:0);
                if (ni < 6 && nj < 6 && nk < 6) {
                    Muscl::predict(qloc[ni][nj][nk], dq[ni][nj][nk][idim], s0, dt_dx, qm_tmp, qp_tmp, nvar);
                    qr_interface[0] = qp_tmp[0];
                    qr_interface[1] = qp_tmp[1+idim];
                    if (NDIM > 1) qr_interface[2] = qp_tmp[1+(idim+1)%NDIM];
                    if (NDIM > 2) qr_interface[3] = qp_tmp[1+(idim+2)%NDIM];
                    qr_interface[NDIM + 1] = qp_tmp[NDIM + 1];

                    real_t f_tmp[20];
                    if (params::iriemann == 2) RiemannSolver::solve_hllc(ql_interface, qr_interface, f_tmp, gamma);
                    else RiemannSolver::solve_llf(ql_interface, qr_interface, f_tmp, gamma);
                    
                    flux[i][j][k][idim][0] = f_tmp[0];
                    flux[i][j][k][idim][1+idim] = f_tmp[1];
                    if (NDIM > 1) flux[i][j][k][idim][1+(idim+1)%NDIM] = f_tmp[2];
                    if (NDIM > 2) flux[i][j][k][idim][1+(idim+2)%NDIM] = f_tmp[3];
                    flux[i][j][k][idim][NDIM+1] = f_tmp[NDIM+1];
                    
                    for (int iv = NDIM + 2; iv < nvar; ++iv) {
                        flux[i][j][k][idim][iv] = (f_tmp[0] > 0) ? f_tmp[0] * ql_interface[iv] : f_tmp[0] * qr_interface[iv];
                    }
                }
            }
        }

        for (int k2 = 0; k2 < (NDIM > 2 ? 2 : 1); ++k2) {
        for (int j2 = 0; j2 < (NDIM > 1 ? 2 : 1); ++j2) {
        for (int i2 = 0; i2 < 2; ++i2) {
            int i = 2 + i2; int j = 2 + j2; int k = 2 + k2;
            int icell_pos = 1 + i2 + 2*j2 + 4*k2;
            if (icell_pos <= constants::twotondim) {
                int ind_cell = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + igrid;
                for (int iv = 1; iv <= nvar; ++iv) {
                    for (int idim = 0; idim < NDIM; ++idim) {
                        int ni = i - (idim==0?1:0); int nj = j - (idim==1?1:0); int nk = k - (idim==2?1:0);
                        grid_.unew(ind_cell, iv) += (flux[ni][nj][nk][idim][iv-1] - flux[i][j][k][idim][iv-1]) * dt_dx;
                    }
                }
            }
        }
        }
        }
    }
}

void HydroSolver::add_gravity_source_terms(int ilevel, real_t dt) {
    if (grid_.count_grids_at_level(ilevel) == 0) return;

    real_t fact = 0.5 * dt;
    const real_t smallr = 1e-10;
    
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ind = 1; ind <= constants::twotondim; ++ind) {
                int ind_cell = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;
                
                real_t d = std::max(grid_.unew(ind_cell, 1), smallr);
                real_t u = 0.0, v = 0.0, w = 0.0;
                
                if (NDIM > 0) u = grid_.unew(ind_cell, 2) / d;
                if (NDIM > 1) v = grid_.unew(ind_cell, 3) / d;
                if (NDIM > 2) w = grid_.unew(ind_cell, 4) / d;
                
                real_t e_kin = 0.5 * d * (u*u + v*v + w*w);
                real_t e_prim = grid_.unew(ind_cell, NDIM + 2) - e_kin; // internal energy + thermal + non-thermal
                
                if (NDIM > 0) {
                    u = u + grid_.f(ind_cell, 1) * fact;
                    grid_.unew(ind_cell, 2) = d * u;
                }
                if (NDIM > 1) {
                    v = v + grid_.f(ind_cell, 2) * fact;
                    grid_.unew(ind_cell, 3) = d * v;
                }
                if (NDIM > 2) {
                    w = w + grid_.f(ind_cell, 3) * fact;
                    grid_.unew(ind_cell, 4) = d * w;
                }
                
                e_kin = 0.5 * d * (u*u + v*v + w*w);
                grid_.unew(ind_cell, NDIM + 2) = e_prim + e_kin;
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void HydroSolver::get_diagnostics(int ilevel, real_t dx, real_t& mind, real_t& maxv) {
    mind = 1e30; maxv = 0.0;
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                if (grid_.son[ind_cell] != 0) continue; // Only leaf cells

                real_t d = grid_.unew(ind_cell, 1);
                mind = std::min(mind, d);
                real_t u = grid_.unew(ind_cell, 2) / std::max(1e-10, d);
                real_t v = grid_.unew(ind_cell, 3) / std::max(1e-10, d);
                real_t w = grid_.unew(ind_cell, 4) / std::max(1e-10, d);
                maxv = std::max(maxv, std::sqrt(u*u + v*v + w*w));
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

} // namespace ramses
