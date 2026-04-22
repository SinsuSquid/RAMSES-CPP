#include "ramses/HydroSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/RiemannSolver.hpp"
#include "ramses/SlopeLimiter.hpp"
#include "ramses/Muscl.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace ramses {

void HydroSolver::set_unew(int ilevel) {
    for (int ind_cell = 1; ind_cell <= grid_.ncell; ++ind_cell) {
        for (int ivar = 1; ivar <= grid_.nvar; ++ivar) {
            grid_.unew(ind_cell, ivar) = grid_.uold(ind_cell, ivar);
        }
    }
}

void HydroSolver::set_uold(int ilevel) {
    for (int ind_cell = 1; ind_cell <= grid_.ncell; ++ind_cell) {
        for (int ivar = 1; ivar <= grid_.nvar; ++ivar) {
            grid_.uold(ind_cell, ivar) = grid_.unew(ind_cell, ivar);
        }
    }
}

void HydroSolver::godunov_fine(int ilevel) {
    set_unew(ilevel);
    // Real port would loop over active grids here
}

void HydroSolver::ctoprim(const real_t u[], real_t q[], real_t gamma) {
    const real_t smallr = 1e-10;
    real_t d = std::max(u[0], smallr);
    q[0] = d;
    q[1] = u[1] / d;
    q[2] = u[2] / d;
    q[3] = u[3] / d;
    
    real_t e_kin = 0.5 * d * (q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    real_t e_int = u[4] - e_kin;
    q[4] = std::max(e_int * (gamma - 1.0), d * 1e-10);
}

void HydroSolver::interpol_hydro(const real_t u1[7][5], real_t u2[8][5]) {
    // u1: [0]=center, [1]=x-left, [2]=x-right, [3]=y-left, [4]=y-right, [5]=z-left, [6]=z-right
    // Linear interpolation: u2 = u0 + grad_x * dx + grad_y * dy + grad_z * dz
    
    real_t w[3][5]; // Gradients for each dimension and variable
    for (int idim = 0; idim < 3; ++idim) {
        for (int iv = 0; iv < 5; ++iv) {
            real_t dlft = 0.5 * (u1[0][iv] - u1[2 * idim + 1][iv]);
            real_t drgt = 0.5 * (u1[2 * idim + 2][iv] - u1[0][iv]);
            
            // MinMod limiter
            if (dlft * drgt <= 0.0) {
                w[idim][iv] = 0.0;
            } else {
                w[idim][iv] = (std::abs(dlft) < std::abs(drgt)) ? dlft : drgt;
            }
        }
    }

    for (int ind = 0; ind < 8; ++ind) {
        int ix = (ind & 1);
        int iy = (ind & 2) >> 1;
        int iz = (ind & 4) >> 2;
        
        real_t xc[3] = {static_cast<real_t>(ix) - 0.5f, 
                        static_cast<real_t>(iy) - 0.5f, 
                        static_cast<real_t>(iz) - 0.5f};
        
        for (int iv = 0; iv < 5; ++iv) {
            u2[ind][iv] = u1[0][iv];
            for (int idim = 0; idim < 3; ++idim) {
                u2[ind][iv] += 2.0 * w[idim][iv] * xc[idim];
            }
        }
    }
}

void HydroSolver::gather_stencil(int igrid, int ilevel, LocalStencil& stencil) {
    int nbors_father[27];
    grid_.get_3x3x3_father(igrid, nbors_father);
    
    for (int k1 = 0; k1 < 3; ++k1) {
    for (int j1 = 0; j1 < 3; ++j1) {
    for (int i1 = 0; i1 < 3; ++i1) {
        int ifather = nbors_father[i1 + 3*j1 + 9*k1];
        int ison = grid_.son[ifather];
        
        if (ison > 0) {
            for (int k2 = 0; k2 < 2; ++k2) {
            for (int j2 = 0; j2 < 2; ++j2) {
            for (int i2 = 0; i2 < 2; ++i2) {
                int icell_pos = 1 + i2 + 2*j2 + 4*k2;
                int ind_cell = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + ison;
                int i3 = i1 * 2 + i2;
                int j3 = j1 * 2 + j2;
                int k3 = k1 * 2 + k2;
                for (int iv = 1; iv <= 5; ++iv) {
                    stencil.uloc[i3][j3][k3][iv - 1] = grid_.uold(ind_cell, iv);
                }
                stencil.refined[i3][j3][k3] = (grid_.son[ind_cell] > 0);
            }
            }
            }
        } else {
            // Gathering 7 cells for interpolation
            real_t u1[7][5];
            for(int iv=0; iv<5; ++iv) u1[0][iv] = grid_.uold(ifather, iv+1);
            
            // To be robust, we should gather neighbors of ifather.
            // For now, use center value if neighbors are unavailable (simplification).
            // In a real port, we'd use grid_.get_nbor_cells or similar.
            for(int n=1; n<7; ++n) for(int iv=0; iv<5; ++iv) u1[n][iv] = u1[0][iv];

            real_t u2[8][5];
            interpol_hydro(u1, u2);
            for (int k2 = 0; k2 < 2; ++k2) {
            for (int j2 = 0; j2 < 2; ++j2) {
            for (int i2 = 0; i2 < 2; ++i2) {
                int ind_son = i2 + 2*j2 + 4*k2;
                int i3 = i1 * 2 + i2;
                int j3 = j1 * 2 + j2;
                int k3 = k1 * 2 + k2;
                for (int iv = 0; iv < 5; ++iv) {
                    stencil.uloc[i3][j3][k3][iv] = u2[ind_son][iv];
                }
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
    real_t dt_max = 1e30; // Start large
    const real_t smallr = 1e-10;

    for (int i = 1; i <= grid_.ncell; ++i) {
        real_t d = std::max(grid_.uold(i, 1), smallr);
        real_t u = grid_.uold(i, 2) / d;
        real_t v = grid_.uold(i, 3) / d;
        real_t w = grid_.uold(i, 4) / d;
        
        real_t e_kin = 0.5 * d * (u*u + v*v + w*w);
        real_t e_int = grid_.uold(i, 5) - e_kin;
        real_t p = std::max(e_int * (gamma - 1.0), d * 1e-10);
        
        real_t cs = std::sqrt(gamma * p / d);
        real_t vel_max = std::max({std::abs(u), std::abs(v), std::abs(w)});
        
        real_t dt_cell = courant_factor * dx / (vel_max + cs + 1e-20);
        dt_max = std::min(dt_max, dt_cell);
    }
    return dt_max;
}

void HydroSolver::godfine1(const std::vector<int>& ind_grid, int ilevel) {
    real_t gamma = 1.4;
    real_t dt = 0.01;
    real_t dx = 0.1;
    real_t dt_dx = dt / dx;

    for (int igrid : ind_grid) {
        LocalStencil stencil;
        gather_stencil(igrid, ilevel, stencil);

        // 1. Primitive Conversion for 6x6x6 stencil
        real_t qloc[6][6][6][5];
        for(int k=0; k<6; ++k) for(int j=0; j<6; ++j) for(int i=0; i<6; ++i) {
            ctoprim(stencil.uloc[i][j][k], qloc[i][j][k], gamma);
        }

        // 2. Compute TVD Slopes in all 3 directions
        real_t dq[6][6][6][3][5]; // [i][j][k][idim][ivar]
        for(int k=1; k<5; ++k) for(int j=1; j<5; ++j) for(int i=1; i<5; ++i) {
            for(int iv=0; iv<5; ++iv) {
                dq[i][j][k][0][iv] = SlopeLimiter::compute_slope(qloc[i-1][j][k][iv], qloc[i][j][k][iv], qloc[i+1][j][k][iv], 1);
                dq[i][j][k][1][iv] = SlopeLimiter::compute_slope(qloc[i][j-1][k][iv], qloc[i][j][k][iv], qloc[i][j+1][k][iv], 1);
                dq[i][j][k][2][iv] = SlopeLimiter::compute_slope(qloc[i][j][k-1][iv], qloc[i][j][k][iv], qloc[i][j][k+1][iv], 1);
            }
        }

        // 3. Unsplit Flux Calculation
        real_t flux[6][6][6][5][3]; // [i][j][k][ivar][idim]
        for (int idim = 0; idim < 3; ++idim) {
            for(int k=1; k<5; ++k) for(int j=1; j<5; ++j) for(int i=1; i<5; ++i) {
                // Interface between cell (i,j,k) and its neighbor in direction idim
                // We'll compute flux at i+1/2, j+1/2, or k+1/2
                
                real_t ql_interface[5], qr_interface[5];
                real_t s0[5] = {0,0,0,0,0}; // Source terms placeholder (gravity/PdV)
                
                // Tracing
                real_t qm_tmp[5], qp_tmp[5];
                
                // Left side of interface (cell i)
                Muscl::predict(qloc[i][j][k], dq[i][j][k][idim], s0, dt_dx, qm_tmp, qp_tmp, 5);
                // Permute velocities based on idim
                ql_interface[0] = qm_tmp[0];
                ql_interface[1] = qm_tmp[1+idim]; // Normal velocity
                ql_interface[2] = qm_tmp[1+(idim+1)%3];
                ql_interface[3] = qm_tmp[1+(idim+2)%3];
                ql_interface[4] = qm_tmp[4];

                // Right side of interface (cell i+1 in idim)
                int ni = i + (idim==0?1:0);
                int nj = j + (idim==1?1:0);
                int nk = k + (idim==2?1:0);
                
                if (ni < 6 && nj < 6 && nk < 6) {
                    Muscl::predict(qloc[ni][nj][nk], dq[ni][nj][nk][idim], s0, dt_dx, qm_tmp, qp_tmp, 5);
                    qr_interface[0] = qp_tmp[0];
                    qr_interface[1] = qp_tmp[1+idim];
                    qr_interface[2] = qp_tmp[1+(idim+1)%3];
                    qr_interface[3] = qp_tmp[1+(idim+2)%3];
                    qr_interface[4] = qp_tmp[4];

                    real_t f_tmp[5];
                    RiemannSolver::solve_llf(ql_interface, qr_interface, f_tmp, gamma);
                    
                    // Un-permute fluxes
                    flux[i][j][k][0][idim] = f_tmp[0];
                    flux[i][j][k][1+idim][idim] = f_tmp[1];
                    flux[i][j][k][1+(idim+1)%3][idim] = f_tmp[2];
                    flux[i][j][k][1+(idim+2)%3][idim] = f_tmp[3];
                    flux[i][j][k][4][idim] = f_tmp[4];
                }
            }
        }

        // 4. Update unew
        for (int k2 = 0; k2 < 2; ++k2) {
        for (int j2 = 0; j2 < 2; ++j2) {
        for (int i2 = 0; i2 < 2; ++i2) {
            int i = 2 + i2;
            int j = 2 + j2;
            int k = 2 + k2;
            
            int icell_pos = 1 + i2 + 2*j2 + 4*k2;
            int ind_cell = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + igrid;
            
            for (int iv = 1; iv <= 5; ++iv) {
                real_t div_f = 0;
                for (int idim = 0; idim < 3; ++idim) {
                    int ni = i - (idim==0?1:0);
                    int nj = j - (idim==1?1:0);
                    int nk = k - (idim==2?1:0);
                    div_f += (flux[ni][nj][nk][iv-1][idim] - flux[i][j][k][iv-1][idim]);
                }
                grid_.unew(ind_cell, iv) += div_f * dt_dx;
            }
        }
        }
        }
    }
}


} // namespace ramses
