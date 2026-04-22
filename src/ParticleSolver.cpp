#include "ramses/ParticleSolver.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace ramses {

void ParticleSolver::move_fine(int ilevel, real_t dt) {
    std::vector<int> ind_part;
    for (int ipart = 1; ipart <= ps_.npartmax; ++ipart) {
        if (ps_.idp[ipart - 1] > 0 && ps_.levelp[ipart - 1] == ilevel) {
            ind_part.push_back(ipart);
            if (ind_part.size() >= 100) {
                move_particles(ind_part, dt);
                ind_part.clear();
            }
        }
    }
    if (!ind_part.empty()) move_particles(ind_part, dt);
}

void ParticleSolver::move_particles(const std::vector<int>& ind_part, real_t dt) {
    const real_t boxlen = 1.0;
    for (int ipart : ind_part) {
        for (int idim = 1; idim <= NDIM; ++idim) {
            real_t& x = ps_.get_xp(ipart, idim);
            real_t v = ps_.get_vp(ipart, idim);
            x += v * dt;
            if (x < 0.0) x += boxlen;
            if (x >= boxlen) x -= boxlen;
        }
    }
}

void ParticleSolver::assign_mass(int ilevel) {
    // Zero out density field
    std::fill(grid_.rho.begin(), grid_.rho.end(), 0.0);

    // CIC Mass Assignment (simplified)
    for (int ipart = 1; ipart <= ps_.npartmax; ++ipart) {
        if (ps_.idp[ipart - 1] > 0 && ps_.levelp[ipart - 1] == ilevel) {
            real_t x = ps_.get_xp(ipart, 1);
            real_t y = ps_.get_xp(ipart, 2);
            real_t z = ps_.get_xp(ipart, 3);
            
            // Find cell index (simplified 1-based index)
            int ix = static_cast<int>(x * params::nx);
            int iy = static_cast<int>(y * params::ny);
            int iz = static_cast<int>(z * params::nz);
            int icell = 1 + ix + iy * params::nx + iz * params::nx * params::ny;
            
            if (icell <= grid_.ncell) {
                grid_.rho[icell] += ps_.mp[ipart - 1];
            }
        }
    }
}

} // namespace ramses
