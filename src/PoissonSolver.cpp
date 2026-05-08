#include "ramses/PoissonSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>

namespace ramses {

PoissonSolver::~PoissonSolver() {}

void PoissonSolver::solve(int ilevel, real_t aexp, real_t omega_m, real_t rho_tot) {
    if (grid_.ncoarse == 0) return;
    if (res.size() != grid_.phi.size()) res.assign(grid_.phi.size(), 0.0);
    
    real_t fourpi = 8.0 * std::atan(1.0); 
    if (config_.get_bool("run_params", "cosmo", false)) {
        real_t nx_loc = 1.0; 
        real_t scale = params::boxlen / nx_loc;
        fourpi = 1.5 * omega_m * aexp * scale;
    }

    // Perform Multigrid V-Cycles
    for (int iter = 0; iter < 10; ++iter) {
        vcycle(ilevel, fourpi, rho_tot);
    }
}

void PoissonSolver::vcycle(int ilevel, real_t fourpi, real_t rho_tot) {
    smooth(ilevel, fourpi, rho_tot);
    if (ilevel > 1) {
        restrict(ilevel);
        vcycle(ilevel - 1, fourpi, rho_tot);
        prolong(ilevel);
    }
    smooth(ilevel, fourpi, rho_tot);
}

void PoissonSolver::smooth(int ilevel, real_t fourpi, real_t rho_tot) {
    real_t dx = params::boxlen / (real_t)(params::nx * (1 << (ilevel - 1)));
    real_t dx2 = dx * dx;
    real_t inv_twondim = 1.0 / static_cast<real_t>(2 * NDIM);

    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            int igridn[7];
            grid_.get_nbor_grids(igrid, igridn);

            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int icell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                int icelln[6];
                grid_.get_nbor_cells(igridn, ic, icelln, igrid);


                real_t nb_sum = 0;
                for (int n = 0; n < 2 * NDIM; ++n) {
                    if (icelln[n] > 0) nb_sum += grid_.phi[icelln[n] - 1];
                    else nb_sum += grid_.phi[icell - 1]; 
                }
                grid_.phi[icell - 1] = (nb_sum - dx2 * fourpi * (grid_.rho[icell - 1] - rho_tot)) * inv_twondim;
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void PoissonSolver::restrict(int ilevel) {
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            int ifather_cell = grid_.father[igrid - 1];
            if (ifather_cell > 0) {
                real_t sum_rho = 0;
                for (int ic = 1; ic <= constants::twotondim; ++ic) {
                    int icell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    sum_rho += grid_.rho[icell - 1];
                }
                grid_.rho[ifather_cell - 1] = sum_rho / static_cast<real_t>(constants::twotondim);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void PoissonSolver::prolong(int ilevel) {
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            int ifather_cell = grid_.father[igrid - 1];
            if (ifather_cell > 0) {
                real_t phi_c = grid_.phi[ifather_cell - 1];
                for (int ic = 1; ic <= constants::twotondim; ++ic) {
                    int icell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    grid_.phi[icell - 1] += phi_c;
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void PoissonSolver::compute_force(int ilevel) {
    real_t dx = params::boxlen / (real_t)(params::nx * (1 << (ilevel - 1)));
    real_t inv_2dx = 0.5 / dx;

    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            int igridn[7];
            grid_.get_nbor_grids(igrid, igridn);

            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int icell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                int icelln[6];
                grid_.get_nbor_cells(igridn, ic, icelln, igrid);


                for (int idim = 1; idim <= NDIM; ++idim) {
                    real_t phi_left = (icelln[2 * (idim - 1)] > 0) ? grid_.phi[icelln[2 * (idim - 1)] - 1] : grid_.phi[icell - 1];
                    real_t phi_right = (icelln[2 * (idim - 1) + 1] > 0) ? grid_.phi[icelln[2 * (idim - 1) + 1] - 1] : grid_.phi[icell - 1];
                    grid_.f(icell, idim) = -(phi_right - phi_left) * inv_2dx;
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

} // namespace ramses
