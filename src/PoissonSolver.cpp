#include "ramses/PoissonSolver.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace ramses {

void PoissonSolver::solve(int ilevel) {
    if (grid_.ncoarse == 0) return;
    if (res.size() != grid_.phi.size()) res.assign(grid_.phi.size(), 0.0);
    
    // Perform Multigrid V-Cycles
    for (int iter = 0; iter < 10; ++iter) {
        vcycle(ilevel);
    }
}

void PoissonSolver::vcycle(int ilevel) {
    smooth(ilevel);
    if (ilevel > 1) {
        restrict(ilevel);
        vcycle(ilevel - 1);
        prolong(ilevel);
    }
    smooth(ilevel);
}

void PoissonSolver::smooth(int ilevel) {
    real_t dx = 0.5 / static_cast<real_t>(1 << (ilevel - 1));
    real_t dx2 = dx * dx;
    real_t inv_6 = 1.0 / 6.0;

    // Loop over active grids at this level
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            int igridn[7];
            grid_.get_nbor_grids(igrid, igridn);

            for (int ic = 1; ic <= 8; ++ic) {
                int icell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                
                // Get neighbors for 7-point stencil
                int icelln[6];
                grid_.get_nbor_cells(igridn, ic, icelln);

                real_t nb_sum = 0;
                for (int n = 0; n < 6; ++n) {
                    if (icelln[n] > 0) nb_sum += grid_.phi[icelln[n]];
                    else nb_sum += grid_.phi[icell]; // Zero-gradient placeholder
                }

                // Gauss-Seidel update: phi = (sum(nb) - dx^2 * rho) / 6
                grid_.phi[icell] = (nb_sum - dx2 * grid_.rho[icell]) * inv_6;
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void PoissonSolver::restrict(int ilevel) {
    // Average fine residuals to coarse RHS
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            int ifather_cell = grid_.father[igrid - 1];
            if (ifather_cell > 0) {
                real_t sum_rho = 0;
                for (int ic = 1; ic <= 8; ++ic) {
                    int icell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    sum_rho += grid_.rho[icell];
                }
                grid_.rho[ifather_cell] = sum_rho * 0.125;
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void PoissonSolver::prolong(int ilevel) {
    // Interpolate coarse potential correction to fine
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            int ifather_cell = grid_.father[igrid - 1];
            if (ifather_cell > 0) {
                real_t phi_c = grid_.phi[ifather_cell];
                for (int ic = 1; ic <= 8; ++ic) {
                    int icell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    grid_.phi[icell] += phi_c; // Basic correction injection
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

} // namespace ramses
