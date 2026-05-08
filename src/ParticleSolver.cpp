#include "ramses/ParticleSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace ramses {

void ParticleSolver::assign_mass_fine(int ilevel) {
    if (grid_.npart == 0) return;

    real_t dx = params::boxlen / (real_t)(params::nx * (1 << (ilevel - 1)));
    real_t vol = std::pow(dx, NDIM);
    
    for (int ip = 0; ip < grid_.npart; ++ip) {
        if (grid_.levelp[ip] != ilevel) continue;

        real_t xp = grid_.xp[0 * grid_.npartmax + ip];
        real_t yp = (NDIM > 1) ? grid_.xp[1 * grid_.npartmax + ip] : 0.5 * params::boxlen;
        real_t zp = (NDIM > 2) ? grid_.xp[2 * grid_.npartmax + ip] : 0.5 * params::boxlen;

        // CIC mass assignment (Simplified: just NGP with safety for now)
        int icell = find_cell_by_coords(xp, yp, zp);
        if (icell > 0) {
            grid_.rho[icell - 1] += grid_.mp[ip] / vol;
        }
    }
}

void ParticleSolver::move_fine(int ilevel, real_t dt) {
    if (grid_.npart == 0) return;

    for (int ip = 0; ip < grid_.npart; ++ip) {
        if (grid_.levelp[ip] != ilevel) continue;

        // 1. Interpolate force to particle position (NGP for now)
        real_t xp = grid_.xp[0 * grid_.npartmax + ip];
        real_t yp = (NDIM > 1) ? grid_.xp[1 * grid_.npartmax + ip] : 0.5 * params::boxlen;
        real_t zp = (NDIM > 2) ? grid_.xp[2 * grid_.npartmax + ip] : 0.5 * params::boxlen;

        int icell = find_cell_by_coords(xp, yp, zp);
        if (icell > 0) {
            for (int idim = 1; idim <= NDIM; ++idim) {
                real_t f = grid_.f(icell, idim);
                // 2. Update velocity (v = v + f*dt)
                grid_.vp[(idim - 1) * grid_.npartmax + ip] += f * dt;
                // 3. Update position (x = x + v*dt)
                grid_.xp[(idim - 1) * grid_.npartmax + ip] += grid_.vp[(idim - 1) * grid_.npartmax + ip] * dt;
            }
        }
        
        // Periodic boundary wrapping
        for (int idim = 0; idim < NDIM; ++idim) {
            if (grid_.xp[idim * grid_.npartmax + ip] < 0) grid_.xp[idim * grid_.npartmax + ip] += params::boxlen;
            if (grid_.xp[idim * grid_.npartmax + ip] >= params::boxlen) grid_.xp[idim * grid_.npartmax + ip] -= params::boxlen;
        }
    }
}

int ParticleSolver::find_cell_by_coords(real_t x, real_t y, real_t z) {
    // 1. Find coarse cell
    real_t dx_coarse = params::boxlen / std::max({params::nx, params::ny, params::nz});
    int ix = std::clamp(static_cast<int>(x / dx_coarse), 0, params::nx - 1);
    int iy = std::clamp(static_cast<int>(y / dx_coarse), 0, params::ny - 1);
    int iz = std::clamp(static_cast<int>(z / dx_coarse), 0, params::nz - 1);
    int curr_cell = iz * params::nx * params::ny + iy * params::nx + ix + 1;
    
    // 2. Traverse tree
    int ilevel = 1;
    while (grid_.son[curr_cell - 1] > 0) {
        int ig = grid_.son[curr_cell - 1];
        ilevel++;
        real_t dx_level = params::boxlen / (real_t)(params::nx * (1 << (ilevel - 1)));
        
        real_t xg = grid_.xg[0 * grid_.ngridmax + ig - 1];
        real_t yg = (NDIM > 1) ? grid_.xg[1 * grid_.ngridmax + ig - 1] : 0.5 * params::boxlen;
        real_t zg = (NDIM > 2) ? grid_.xg[2 * grid_.ngridmax + ig - 1] : 0.5 * params::boxlen;
        
        int ixc = (x > xg) ? 1 : 0;
        int iyc = (y > yg) ? 1 : 0;
        int izc = (z > zg) ? 1 : 0;
        int ic = izc * 4 + iyc * 2 + ixc + 1;
        curr_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
    }
    return curr_cell;
}

} // namespace ramses
