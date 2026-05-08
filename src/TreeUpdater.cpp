#include "ramses/TreeUpdater.hpp"
#include "ramses/AmrGrid.hpp"
#include "ramses/Constants.hpp"
#include "ramses/MpiManager.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

namespace ramses {

TreeUpdater::TreeUpdater(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

int get_nbor_of_coarse(const AmrGrid& grid, int ic, int idim, int side) {
    int nx = grid.nx, ny = grid.ny, nz = grid.nz;
    int idx = ic - 1;
    int iz = idx / (nx * ny); int rem = idx % (nx * ny);
    int iy = rem / nx; int ix = rem % nx;
    int ixyz[3] = {ix, iy, iz};
    int n[3] = {nx, ny, nz};
    if (side == 0) ixyz[idim] = (ixyz[idim] - 1 + n[idim]) % n[idim];
    else ixyz[idim] = (ixyz[idim] + 1) % n[idim];
    return ixyz[2] * nx * ny + ixyz[1] * nx + ixyz[0] + 1;
}

void TreeUpdater::make_grid_fine(int ilevel) {
    int myid = MpiManager::instance().rank() + 1, n2d = (1 << NDIM);
    std::vector<int> cells_to_refine;
    
    if (ilevel == 1) {
        for (int ic = 1; ic <= grid_.ncoarse; ++ic) {
            if (grid_.flag1[ic - 1] == 1 && grid_.son[ic - 1] == 0) cells_to_refine.push_back(ic);
        }
    } else {
        int ig = grid_.get_headl(myid, ilevel - 1);
        while (ig > 0) {
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                if (grid_.flag1[idc] == 1 && grid_.son[idc] == 0) cells_to_refine.push_back(idc + 1);
            }
            ig = grid_.next[ig - 1];
        }
    }

    for (int ic_coarse : cells_to_refine) {
        int idc = ic_coarse - 1;
        int new_ig = grid_.get_free_grid(); if (new_ig == 0) return;
        grid_.son[idc] = new_ig; grid_.father[new_ig - 1] = idc + 1;
        
        if (ilevel == 1) {
            for (int idim = 0; idim < 3; ++idim) {
                for (int side = 0; side < 2; ++side) {
                    grid_.nbor[(idim * 2 + side) * grid_.ngridmax + new_ig - 1] = get_nbor_of_coarse(grid_, ic_coarse, idim, side);
                }
            }
            int ixyz[3], idx = ic_coarse - 1; ixyz[2] = idx / (grid_.nx * grid_.ny); idx %= (grid_.nx * grid_.ny); ixyz[1] = idx / grid_.nx; ixyz[0] = idx % grid_.nx;
            real_t dx_coarse = grid_.boxlen / std::max({grid_.nx, grid_.ny, grid_.nz});
            for (int d = 1; d <= NDIM; ++d) grid_.xg[(d - 1) * grid_.ngridmax + new_ig - 1] = (ixyz[d - 1] + 0.5) * dx_coarse;
        } else {
            int ig = ((idc - grid_.ncoarse) % grid_.ngridmax) + 1;
            int ic = ((idc - grid_.ncoarse) / grid_.ngridmax) + 1;
            int ign[7]; grid_.get_nbor_grids(ig, ign);
            int icn_nb[6]; grid_.get_nbor_cells(ign, ic, icn_nb, ig);
            for (int i = 0; i < 6; ++i) grid_.nbor[i * grid_.ngridmax + new_ig - 1] = icn_nb[i];
            for (int d = 1; d <= NDIM; ++d) {
                real_t dx_level = grid_.boxlen / (real_t)(grid_.nx * (1 << (ilevel - 1)));
                int ix = (ic - 1) & 1, iy = ((ic - 1) & 2) >> 1, iz = ((ic - 1) & 4) >> 2; int ixyz[3] = {ix, iy, iz};
                grid_.xg[(d - 1) * grid_.ngridmax + (new_ig - 1)] = grid_.xg[(d - 1) * grid_.ngridmax + (ig - 1)] + (ixyz[d - 1] - 0.5) * dx_level * 0.5;
            }
        }
        
        real_t u1[7][64] = {0}, u2[8][64] = {0};
        for(int iv=1; iv<=grid_.nvar; ++iv) u1[0][iv-1] = grid_.uold(ic_coarse, iv);
        int icn_ref[6] = {0};
        if (ilevel > 1) {
            int ig_ref = ((ic_coarse - 1 - grid_.ncoarse) % grid_.ngridmax) + 1;
            int ic_ref = ((ic_coarse - 1 - grid_.ncoarse) / grid_.ngridmax) + 1;
            int ign_ref[7] = {0}; grid_.get_nbor_grids(ig_ref, ign_ref);
            grid_.get_nbor_cells(ign_ref, ic_ref, icn_ref, ig_ref);
        } else {
            for(int idim=0; idim<NDIM; ++idim) {
                icn_ref[idim*2] = get_nbor_of_coarse(grid_, ic_coarse, idim, 0);
                icn_ref[idim*2+1] = get_nbor_of_coarse(grid_, ic_coarse, idim, 1);
            }
        }
        for(int idim=0; idim<NDIM; ++idim) {
            int id_l = icn_ref[idim*2], id_r = icn_ref[idim*2+1];
            if (id_l > 0) for(int iv=1; iv<=grid_.nvar; ++iv) u1[2*idim+1][iv-1] = grid_.uold(id_l, iv); else for(int iv=1; iv<=grid_.nvar; ++iv) u1[2*idim+1][iv-1] = u1[0][iv-1];
            if (id_r > 0) for(int iv=1; iv<=grid_.nvar; ++iv) u1[2*idim+2][iv-1] = grid_.uold(id_r, iv); else for(int iv=1; iv<=grid_.nvar; ++iv) u1[2*idim+2][iv-1] = u1[0][iv-1];
        }
        if (interpol_hook_) interpol_hook_(u1, u2); else { for(int i=0; i<8; ++i) for(int iv=0; iv<grid_.nvar; ++iv) u2[i][iv] = u1[0][iv]; }
        
        for (int isc = 1; isc <= n2d; ++isc) {
            int id_child = grid_.ncoarse + (isc - 1) * grid_.ngridmax + new_ig - 1;
            for (int iv = 1; iv <= grid_.nvar; ++iv) {
                grid_.uold_vec[(iv-1)*grid_.ncell + id_child] = u2[isc-1][iv-1];
                grid_.unew_vec[(iv-1)*grid_.ncell + id_child] = u2[isc-1][iv-1];
            }
        }
        grid_.add_to_level_list(new_ig, ilevel);
    }
}

void TreeUpdater::remove_grid_fine(int ilevel) {
    int myid = 1, n2d = (1 << NDIM);
    if (ilevel >= grid_.nlevelmax) return;
    int ig = grid_.get_headl(myid, ilevel);
    while (ig > 0) {
        int next_ig = grid_.next[ig - 1], id_p = grid_.father[ig - 1];
        if (id_p > 0) {
            bool has_c = false;
            for (int ic = 1; ic <= n2d; ++ic) if (grid_.son[grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1] != 0) { has_c = true; break; }
            if (!has_c && grid_.flag1[id_p - 1] == 0 && ilevel > config_.get_int("amr_params", "levelmin", 1)) {
                grid_.son[id_p - 1] = 0;
                int p = grid_.prev[ig - 1], n = grid_.next[ig - 1];
                if (p > 0) grid_.next[p - 1] = n; else grid_.headl_vec[(ilevel - 1) * grid_.ncpu + (myid - 1)] = n;
                if (n > 0) grid_.prev[n - 1] = p; else grid_.taill_vec[(ilevel - 1) * grid_.ncpu + (myid - 1)] = p;
                grid_.numbl_vec[(ilevel - 1) * grid_.ncpu + (myid - 1)]--;
                grid_.free_grid(ig);
            }
        }
        ig = next_ig;
    }
}

void TreeUpdater::restrict_fine(int ilevel) {
    int myid = 1, n2d = (1 << NDIM);
    int ig = grid_.get_headl(myid, ilevel);
    real_t twotondim = std::pow(2.0, NDIM);
    while (ig > 0) {
        int id_p = grid_.father[ig - 1];
        if (id_p > 0) {
            for (int iv = 1; iv <= grid_.nvar; ++iv) {
                real_t sum = 0;
                for (int ic = 1; ic <= n2d; ++ic) sum += grid_.uold(grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1, iv);
                grid_.uold(id_p, iv) = sum / twotondim;
            }
        }
        ig = grid_.next[ig - 1];
    }
}

void TreeUpdater::flag_fine(int ilevel, real_t ed, real_t ep, real_t ev, real_t eb2, const std::vector<real_t>& evar, int nexp) {
    int myid = 1, n2d = (1 << NDIM), lmin = config_.get_int("amr_params", "levelmin", 1);
    if (ilevel == 1) {
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            if (ilevel < lmin) { grid_.flag1[i-1] = 1; continue; }
            grid_.flag1[i-1] = 0;
            real_t d = grid_.uold(i, 1), p = (grid_.uold(i, 5) - 0.5*d*(0.0))*(grid_.gamma-1.0); // Simplified
            // basic refinement criteria
            if (ed > 0 && d > ed) grid_.flag1[i-1] = 1;
        }
    } else {
        int ig = grid_.get_headl(myid, ilevel - 1);
        while (ig > 0) {
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                if (ilevel < lmin) { grid_.flag1[idc] = 1; continue; }
                grid_.flag1[idc] = 0;
                real_t d = grid_.uold(idc + 1, 1);
                if (ed > 0 && d > ed) grid_.flag1[idc] = 1;
            }
            ig = grid_.next[ig - 1];
        }
    }
}

} // namespace ramses
