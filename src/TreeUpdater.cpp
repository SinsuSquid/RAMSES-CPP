#include "ramses/TreeUpdater.hpp"
#include "ramses/AmrGrid.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

namespace ramses {

TreeUpdater::TreeUpdater(AmrGrid& grid) : grid_(grid) {}

// Helper for coarse neighbors
static int get_nbor_of_coarse(const AmrGrid& grid, int icell, int idim, int side) {
    int idx = icell - 1;
    int nx = grid.nx, ny = grid.ny, nz = grid.nz;
    int iz = idx / (nx * ny);
    int iy = (idx % (nx * ny)) / nx;
    int ix = idx % nx;
    
    if (idim == 0) {
        if (side == 0) ix = (ix == 0) ? nx - 1 : ix - 1;
        else ix = (ix == nx - 1) ? 0 : ix + 1;
    } else if (idim == 1) {
        if (side == 0) iy = (iy == 0) ? ny - 1 : iy - 1;
        else iy = (iy == ny - 1) ? 0 : iy + 1;
    } else if (idim == 2) {
        if (side == 0) iz = (iz == 0) ? nz - 1 : iz - 1;
        else iz = (iz == nz - 1) ? 0 : iz + 1;
    }
    return 1 + ix + iy * nx + iz * nx * ny;
}

void TreeUpdater::make_grid_fine(int ilevel) {
    if (ilevel >= grid_.nlevelmax) return;

    std::vector<int> grids_to_refine;
    std::vector<int> coarse_cells_to_refine;

    if (ilevel == 1) {
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            if (grid_.flag1[i-1] > 0 && grid_.son[i-1] == 0) {
                coarse_cells_to_refine.push_back(i);
            }
        }
    } else {
        int myid = grid_.ncpu > 0 ? grid_.ncpu : 1;
        for (int icpu = 1; icpu <= myid; ++icpu) {
            int igrid = grid_.get_headl(icpu, ilevel);
            while (igrid > 0) {
                bool refine = false;
                for (int ic = 1; ic <= constants::twotondim; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid - 1;
                    if (grid_.flag1[idc] > 0 && grid_.son[idc] == 0) { refine = true; break; }
                }
                if (refine) grids_to_refine.push_back(igrid);
                igrid = grid_.next[igrid - 1];
            }
        }
    }

    // Refine coarse cells
    for (int ic_coarse : coarse_cells_to_refine) {
        int new_ig = grid_.get_free_grid();
        if (new_ig == 0) return;
        grid_.son[ic_coarse - 1] = new_ig;
        grid_.father[new_ig - 1] = ic_coarse;
        // Set neighbors for new level 2 grid
        for (int idim = 0; idim < 3; ++idim) {
            for (int side = 0; side < 2; ++side) {
                grid_.nbor[(idim * 2 + side) * grid_.ngridmax + new_ig - 1] = get_nbor_of_coarse(grid_, ic_coarse, idim, side);
            }
        }
        // Coordinates for child grid
        for (int d = 1; d <= NDIM; ++d) {
            int ixyz[3];
            int idx = ic_coarse - 1;
            ixyz[2] = idx / (grid_.nx * grid_.ny); idx %= (grid_.nx * grid_.ny);
            ixyz[1] = idx / grid_.nx;
            ixyz[0] = idx % grid_.nx;
            real_t dx_coarse = grid_.boxlen / std::max({grid_.nx, grid_.ny, grid_.nz});
            grid_.xg[(d-1)*grid_.ngridmax + new_ig - 1] = (ixyz[d-1] + 0.5) * dx_coarse;
        }
        // Data interpolation (simple injection for now)
        for (int isc = 1; isc <= constants::twotondim; ++isc) {
            int id_child = grid_.ncoarse + (isc - 1) * grid_.ngridmax + new_ig - 1;
            for (int iv = 1; iv <= grid_.nvar; ++iv) {
                grid_.uold_vec[(iv-1)*grid_.ncell + id_child] = grid_.uold_vec[(iv-1)*grid_.ncell + ic_coarse - 1];
                grid_.unew_vec[(iv-1)*grid_.ncell + id_child] = grid_.uold_vec[(iv-1)*grid_.ncell + ic_coarse - 1];
            }
        }
        grid_.add_to_level_list(new_ig, ilevel + 1);
    }

    // Refine grids
    for (int ig : grids_to_refine) {
        int ign[7]; grid_.get_nbor_grids(ig, ign);
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
            if (grid_.flag1[idc] > 0 && grid_.son[idc] == 0) {
                int new_ig = grid_.get_free_grid();
                if (new_ig == 0) return;
                grid_.son[idc] = new_ig;
                grid_.father[new_ig - 1] = idc + 1; // 1-based
                // Set neighbors for new child grid
                int icn_nb[6]; grid_.get_nbor_cells(ign, ic, icn_nb, ig);
                for (int i = 0; i < 6; ++i) {
                    grid_.nbor[i * grid_.ngridmax + new_ig - 1] = icn_nb[i];
                }
                for (int d = 1; d <= NDIM; ++d) {
                    real_t dx_level = grid_.boxlen / (real_t)(grid_.nx * (1 << (ilevel - 1)));
                    int ix = (ic - 1) & 1, iy = ((ic - 1) & 2) >> 1, iz = ((ic - 1) & 4) >> 2;
                    int ixyz[3] = {ix, iy, iz};
                    grid_.xg[(d - 1) * grid_.ngridmax + (new_ig - 1)] = grid_.xg[(d - 1) * grid_.ngridmax + (ig - 1)] + (ixyz[d-1] - 0.5) * dx_level;
                }
                // Data interpolation
                int icn[6]; grid_.get_nbor_cells(ign, ic, icn, ig);
                real_t u1[7][20] = {0};
                for(int iv=1; iv<=grid_.nvar; ++iv) u1[0][iv-1] = grid_.uold(idc+1, iv);
                for(int idim=0; idim<NDIM; ++idim) {
                    int id_l = icn[idim*2], id_r = icn[idim*2+1];
                    if (id_l > 0) for(int iv=1; iv<=grid_.nvar; ++iv) u1[2*idim+1][iv-1] = grid_.uold(id_l, iv);
                    else for(int iv=1; iv<=grid_.nvar; ++iv) u1[2*idim+1][iv-1] = u1[0][iv-1];
                    if (id_r > 0) for(int iv=1; iv<=grid_.nvar; ++iv) u1[2*idim+2][iv-1] = grid_.uold(id_r, iv);
                    else for(int iv=1; iv<=grid_.nvar; ++iv) u1[2*idim+2][iv-1] = u1[0][iv-1];
                }
                real_t u2[8][20] = {0};
                if (interpol_hook_) interpol_hook_(u1, u2);
                else { for(int i=0; i<8; ++i) for(int iv=0; iv<grid_.nvar; ++iv) u2[i][iv] = u1[0][iv]; }
                for (int isc = 1; isc <= constants::twotondim; ++isc) {
                    int id_child = grid_.ncoarse + (isc - 1) * grid_.ngridmax + new_ig - 1;
                    for (int iv = 1; iv <= grid_.nvar; ++iv) {
                        grid_.uold_vec[(iv-1)*grid_.ncell + id_child] = u2[isc-1][iv-1];
                        grid_.unew_vec[(iv-1)*grid_.ncell + id_child] = u2[isc-1][iv-1];
                    }
                }
                grid_.add_to_level_list(new_ig, ilevel + 1);
            }
        }
    }
}

void TreeUpdater::restrict_fine(int ilevel) {
    if (ilevel >= grid_.nlevelmax) return;
    int myid = grid_.ncpu > 0 ? grid_.ncpu : 1;
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid - 1;
            int ig_son = grid_.son[idc];
            if (ig_son > 0) {
                for (int iv = 1; iv <= grid_.nvar; ++iv) {
                    real_t sum = 0;
                    for (int is = 1; is <= constants::twotondim; ++is) {
                        sum += grid_.uold(grid_.ncoarse + (is - 1) * grid_.ngridmax + ig_son, iv);
                    }
                    grid_.uold(idc + 1, iv) = sum / constants::twotondim;
                }
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

void TreeUpdater::flag_fine(int ilevel, real_t err_grad_d, real_t err_grad_p, real_t err_grad_v) {
    if (ilevel == 1) {
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            grid_.flag1[i-1] = 0;
            if (grid_.son[i-1] != 0) continue;
            // Simplified flagging for coarse cells (always refine to levelmin?)
        }
    } else {
        int myid = grid_.ncpu > 0 ? grid_.ncpu : 1;
        int igrid = grid_.get_headl(myid, ilevel);
        while (igrid > 0) {
            int ign[7]; grid_.get_nbor_grids(igrid, ign);
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid - 1;
                grid_.flag1[idc] = 0;
                if (grid_.son[idc] != 0) continue;
                int icn[6]; grid_.get_nbor_cells(ign, ic, icn, igrid);
                real_t d_c = grid_.uold(idc + 1, 1);
                for (int idim = 0; idim < NDIM; ++idim) {
                    int idl = icn[idim * 2], idr = icn[idim * 2 + 1];
                    if (idl > 0 && idr > 0) {
                        real_t d_l = grid_.uold(idl, 1), d_r = grid_.uold(idr, 1);
                        real_t grad = std::abs(d_r - d_l) / std::max({d_l, d_r, 1e-10});
                        if (grad > err_grad_d) { 
                            grid_.flag1[idc] = 1; 
                            break; 
                        }
                    }
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

} // namespace ramses
