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
    // Top level entry
    set_unew(ilevel);
    
    // In a real port, we would loop over active grids.
    // For this 1D POC, we'll just demonstrate the logical flow.
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
    q[4] = std::max(e_int * (gamma - 1.0), d * 1e-10); // pressure floor
}

void HydroSolver::godfine1(const std::vector<int>& ind_grid, int ilevel) {
    // ind_grid contains 1-based oct indices
    real_t gamma = 1.4;
    
    for (int igrid : ind_grid) {
        // For each oct, we need a 6x6x6 stencil
        // Inner cells of oct: (ix,iy,iz) = [1,2] x [1,2] x [1,2] (in a 4x4x4 block)
        // RAMSES uses IU1=-1, IU2=4 for buffer
        
        // This is a complex assembly step. For the port, we'll implement 
        // a simplified 1D sweep to prove the flow.
        
        for (int icell_pos = 1; icell_pos <= constants::twotondim; ++icell_pos) {
            int ind_cell = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + igrid;
            
            // 1. Get neighbors
            int igridn[7];
            grid_.get_nbor_grids(igrid, igridn);
            int icelln[6];
            grid_.get_nbor_cells(igridn, icell_pos, icelln);
            
            // 2. Compute fluxes in each direction
            for (int idim = 0; idim < NDIM; ++idim) {
                int left_cell = icelln[idim * 2];
                int right_cell = icelln[idim * 2 + 1];
                
                if (left_cell > 0 && right_cell > 0) {
                    real_t ql[5], qc[5], qr[5];
                    ctoprim(grid_.uold.data() + (left_cell - 1) * 5, ql, gamma);
                    ctoprim(grid_.uold.data() + (ind_cell - 1) * 5, qc, gamma);
                    ctoprim(grid_.uold.data() + (right_cell - 1) * 5, qr, gamma);
                    
                    // MUSCL Reconstruction
                    real_t dql[5], dqr[5];
                    for (int iv = 0; iv < 5; ++iv) {
                        // Left interface state
                        real_t slope_l = SlopeLimiter::compute_slope(0, ql[iv], qc[iv], 1); // Simplification
                        real_t slope_r = SlopeLimiter::compute_slope(ql[iv], qc[iv], qr[iv], 1);
                        
                        // Prediction step (dt=0 for now)
                        // ...
                    }
                }
            }
        }
    }
}

} // namespace ramses
