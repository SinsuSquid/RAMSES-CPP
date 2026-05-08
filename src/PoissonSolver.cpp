#include "ramses/PoissonSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>

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

            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int icell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                
                // Get neighbors for 7-point stencil
                int icelln[6];
                grid_.get_nbor_cells(igridn, ic, icelln, igrid);

                real_t nb_sum = 0;
                for (int n = 0; n < constants::twondim; ++n) {
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
                for (int ic = 1; ic <= constants::twotondim; ++ic) {
                    int icell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    sum_rho += grid_.rho[icell];
                }
                grid_.rho[ifather_cell] = sum_rho / static_cast<real_t>(constants::twotondim);
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
                for (int ic = 1; ic <= constants::twotondim; ++ic) {
                    int icell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    grid_.phi[icell] += phi_c; // Basic correction injection
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void PoissonSolver::compute_force(int ilevel) {
    int gravity_type = config_.get_int("poisson_params", "gravity_type", 0);
    if (gravity_type == 3) {
        std::string param_str = config_.get("poisson_params", "gravity_params", "");
        real_t a1 = 1.42e-3;
        real_t a2 = 5.49e-4;
        real_t z0 = 0.18e3;
        
        // Very basic parsing for 3 floats
        if (!param_str.empty()) {
            std::replace(param_str.begin(), param_str.end(), 'd', 'e');
            std::replace(param_str.begin(), param_str.end(), 'D', 'e');
            std::replace(param_str.begin(), param_str.end(), ',', ' ');
            std::stringstream ss(param_str);
            ss >> a1 >> a2 >> z0;
        }
        
        real_t scale_l = config_.get_double("units_params", "units_length", 1.0);
        real_t scale_t = config_.get_double("units_params", "units_time", 1.0);
        
        real_t kpc2cm = 3.085677581282e21;
        real_t pc2cm = 3.085677581282e18;
        real_t Myr2sec = 3.15576e13;

        a1 = a1 * kpc2cm / (Myr2sec * Myr2sec) / scale_l * (scale_t * scale_t);
        a2 = a2 / (Myr2sec * Myr2sec) * (scale_t * scale_t);
        z0 = z0 * pc2cm / scale_l;

        if (ilevel == 1) {
            for (int idc = 1; idc <= grid_.ncoarse; ++idc) {
                // TODO: calculate x for coarse cells
                for (int idim = 1; idim <= NDIM; ++idim) grid_.f(idc, idim) = 0.0;
            }
        }

        for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
            int igrid = (ilevel > 1) ? grid_.headl(icpu, ilevel - 1) : 0;
            while (igrid > 0) {
                for (int ic = 1; ic <= 8; ++ic) {
                    if (ic > constants::twotondim) break;
                    int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    
                    real_t x = grid_.xg[(NDIM - 1) * grid_.ngridmax + igrid - 1];
                    int ix = (ic - 1) & 1;
                    int iy = ((ic - 1) & 2) >> 1;
                    int iz = ((ic - 1) & 4) >> 2;
                    real_t dx = 0.5 / static_cast<real_t>(1 << (ilevel - 1));
                    if (NDIM == 1) x += (static_cast<real_t>(ix) - 0.5) * dx;
                    else if (NDIM == 2) x += (static_cast<real_t>(iy) - 0.5) * dx;
                    else if (NDIM == 3) x += (static_cast<real_t>(iz) - 0.5) * dx;

                    real_t rz = x * params::boxlen - 0.5 * params::boxlen;
                    grid_.f(ind_cell, NDIM) = -a1 * rz / std::sqrt(rz * rz + z0 * z0) - a2 * rz;
                    for (int idim = 1; idim < NDIM; ++idim) {
                        grid_.f(ind_cell, idim) = 0.0;
                    }
                }
                igrid = grid_.next[igrid - 1];
            }
        }
    } else {
        // -grad(phi) would be implemented here for self-gravity.
    }
}

} // namespace ramses
