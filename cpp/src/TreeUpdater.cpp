#include "ramses/TreeUpdater.hpp"
#include "ramses/Parameters.hpp"
#include <iostream>

namespace ramses {

void TreeUpdater::refine_coarse() {
    // Simplified loop over coarse grid
    for (int k = 0; k < params::nz; ++k) {
        for (int j = 0; j < params::ny; ++j) {
            for (int i = 0; i < params::nx; ++i) {
                int ind = 1 + i + j * params::nx + k * params::nx * params::ny;
                if (grid_.flag2[ind] == 1 && grid_.flag1[ind] == 1 && grid_.son[ind] == 0) {
                    make_grid_coarse(ind, 0, false);
                }
            }
        }
    }
}

void TreeUpdater::make_grid_coarse(int ind_cell, int ibound, bool boundary_region) {
    if (grid_.numbf <= 0) {
        std::cerr << "No more free memory for grid refinement!" << std::endl;
        return; 
    }

    // 1. Pop from free list
    int igrid = grid_.headf;
    grid_.headf = grid_.next[igrid - 1];
    if (grid_.headf > 0) {
        grid_.prev[grid_.headf - 1] = 0;
    }
    grid_.numbf--;

    // 2. Compute grid center
    int nxny = params::nx * params::ny;
    int iz = (ind_cell - 1) / nxny;
    int iy = (ind_cell - 1 - iz * nxny) / params::nx;
    int ix = (ind_cell - 1 - iz * nxny - iy * params::nx);
    
    if (NDIM > 0) grid_.get_xg(igrid, 1) = static_cast<real_t>(ix) + 0.5;
    if (NDIM > 1) grid_.get_xg(igrid, 2) = static_cast<real_t>(iy) + 0.5;
    if (NDIM > 2) grid_.get_xg(igrid, 3) = static_cast<real_t>(iz) + 0.5;

    // 3. Connect to father cell
    grid_.son[ind_cell] = igrid;
    grid_.father[igrid - 1] = ind_cell;

    // 4. Connect to neighboring father cells (with periodic boundaries)
    int ixn[6], iyn[6], izn[6];
    for (int n = 0; n < 6; ++n) {
        ixn[n] = ix; iyn[n] = iy; izn[n] = iz;
    }

    if (params::nx > 0) {
        ixn[0] = (ix > 0) ? ix - 1 : params::nx - 1; // Left
        ixn[1] = (ix < params::nx - 1) ? ix + 1 : 0; // Right
    }
    if (params::ny > 1) {
        iyn[2] = (iy > 0) ? iy - 1 : params::ny - 1;
        iyn[3] = (iy < params::ny - 1) ? iy + 1 : 0;
    }
    if (params::nz > 1) {
        izn[4] = (iz > 0) ? iz - 1 : params::nz - 1;
        izn[5] = (iz < params::nz - 1) ? iz + 1 : 0;
    }

    for (int n = 1; n <= constants::twondim; ++n) {
        int idim = (n - 1) / 2;
        int side = (n - 1) % 2;
        int n_ix = ix, n_iy = iy, n_iz = iz;
        if (idim == 0) n_ix = (side == 0) ? ixn[0] : ixn[1];
        if (idim == 1) n_iy = (side == 0) ? iyn[2] : iyn[3];
        if (idim == 2) n_iz = (side == 0) ? izn[4] : izn[5];
        
        grid_.nbor[(n - 1) * grid_.ngridmax + (igrid - 1)] = 1 + n_ix + n_iy * params::nx + n_iz * nxny;
    }

    // 5. Update cpu_map (simplified for now: match father)
    int icpu = grid_.cpu_map[ind_cell];
    for (int j = 1; j <= constants::twotondim; ++j) {
        int iskip = grid_.ncoarse + (j - 1) * grid_.ngridmax;
        grid_.cpu_map[iskip + igrid] = icpu;
    }

    // 6. Connect to level 1 linked list
    if (!boundary_region) {
        if (grid_.numbl(icpu, 1) > 0) {
            int tail = grid_.taill(icpu, 1);
            grid_.next[igrid - 1] = 0;
            grid_.prev[igrid - 1] = tail;
            grid_.next[tail - 1] = igrid;
            grid_.taill(icpu, 1) = igrid;
            grid_.numbl(icpu, 1)++;
        } else {
            grid_.next[igrid - 1] = 0;
            grid_.prev[igrid - 1] = 0;
            grid_.headl(icpu, 1) = igrid;
            grid_.taill(icpu, 1) = igrid;
            grid_.numbl(icpu, 1) = 1;
        }
    }
}

void TreeUpdater::mark_cells(int ilevel) {
    const real_t err_grad_p = 0.05; // Threshold for pressure gradient
    const real_t floor_p = 1e-10;
    const real_t gamma = 1.4;

    for (int i = 1; i <= grid_.ncell; ++i) {
        // Simplified gradient check: if neighbors were available, we'd check (P_right - P_left)/P_center
        // For now, let's mark cells with high pressure (the blast center)
        
        real_t d = std::max(grid_.uold(i, 1), 1e-10);
        real_t e_kin = 0.5 * (grid_.uold(i, 2)*grid_.uold(i, 2) + 
                             grid_.uold(i, 3)*grid_.uold(i, 3) + 
                             grid_.uold(i, 4)*grid_.uold(i, 4)) / d;
        real_t p = (grid_.uold(i, 5) - e_kin) * (gamma - 1.0);

        if (p > 1e-3) { // Blast region
            grid_.flag1[i] = 1;
            grid_.flag2[i] = 1;
        } else {
            grid_.flag1[i] = 0;
            grid_.flag2[i] = 0;
        }
    }
}


void TreeUpdater::make_grid_fine(int ind_grid_father, int icell_pos, int ilevel, int ibound, bool boundary_region) {
    // Placeholder for Phase 4 extension
}

} // namespace ramses
