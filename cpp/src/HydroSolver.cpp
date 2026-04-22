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
    // Straight injection for now
    for (int i = 0; i < 8; ++i) {
        for (int iv = 0; iv < 5; ++iv) {
            u2[i][iv] = u1[0][iv];
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
            // Placeholder: Interpolation from father cell
            real_t u1[7][5];
            for(int iv=0; iv<5; ++iv) u1[0][iv] = grid_.uold(ifather, iv+1);
            // Neighboring fathers not gathered here for simplicity
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

        // 2. Slopes in X direction (Simplified for 1D POC)
        real_t dq[6][6][6][5];
        for(int k=2; k<4; ++k) for(int j=2; j<4; ++j) for(int i=1; i<5; ++i) {
            for(int iv=0; iv<5; ++iv) {
                dq[i][j][k][iv] = SlopeLimiter::compute_slope(qloc[i-1][j][k][iv], qloc[i][j][k][iv], qloc[i+1][j][k][iv], 1);
            }
        }

        // 3. Trace & Riemann in X
        // Interfaces in X: i=2.5, 3.5, 4.5
        real_t flux_x[6][6][6][5]; // Flux at interface i+1/2
        for(int k=2; k<4; ++k) for(int j=2; j<4; ++j) for(int i=1; i<4; ++i) {
            real_t qm[5], qp[5], ql_interface[5], qr_interface[5];
            real_t s0[5] = {0,0,0,0,0}; // Source terms placeholder
            
            // Predict state at i+1/2 from cell i
            Muscl::predict(qloc[i][j][k], dq[i][j][k], s0, dt_dx, qm, qp, 5);
            for(int iv=0; iv<5; ++iv) ql_interface[iv] = qm[iv];

            // Predict state at i+1/2 from cell i+1
            Muscl::predict(qloc[i+1][j][k], dq[i+1][j][k], s0, dt_dx, qm, qp, 5);
            for(int iv=0; iv<5; ++iv) qr_interface[iv] = qp[iv];

            RiemannSolver::solve_llf(ql_interface, qr_interface, flux_x[i][j][k], gamma);
        }

        // 4. Update unew (Oct-internal cells are i,j,k = 2,3)
        for (int k2 = 0; k2 < 2; ++k2) {
        for (int j2 = 0; j2 < 2; ++j2) {
        for (int i2 = 0; i2 < 2; ++i2) {
            int i = 2 + i2;
            int j = 2 + j2;
            int k = 2 + k2;
            
            int icell_pos = 1 + i2 + 2*j2 + 4*k2;
            int ind_cell = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + igrid;
            
            for (int iv = 0; iv < 5; ++iv) {
                grid_.unew(ind_cell, iv + 1) += (flux_x[i-1][j][k][iv] - flux_x[i][j][k][iv]) * dt_dx;
            }
        }
        }
        }
    }
}

} // namespace ramses
