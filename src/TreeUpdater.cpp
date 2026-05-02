#include "ramses/TreeUpdater.hpp"
#include "ramses/AmrGrid.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace ramses {

TreeUpdater::TreeUpdater(AmrGrid& grid) : grid_(grid) {}

void TreeUpdater::make_grid_fine(int ilevel) {
    if (ilevel >= grid_.nlevelmax) return;

    std::vector<int> grids_to_refine;
    int myid = grid_.ncpu > 0 ? grid_.ncpu : 1; // Simplification for serial
    for (int icpu = 1; icpu <= myid; ++icpu) {
        int igrid = grid_.get_headl(icpu, ilevel);
        while (igrid > 0) {
            bool refine = false;
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                if (grid_.flag1[idc] > 0 && grid_.son[idc] == 0) {
                    refine = true;
                    break;
                }
            }
            if (refine) grids_to_refine.push_back(igrid);
            igrid = grid_.next[igrid - 1];
        }
    }

    for (int ig : grids_to_refine) {
        int ign[7]; grid_.get_nbor_grids(ig, ign);
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
            if (grid_.flag1[idc] > 0 && grid_.son[idc] == 0) {
                int new_ig = grid_.get_free_grid();
                if (new_ig == 0) {
                    std::cerr << "[TreeUpdater] Error: No free grids available!" << std::endl;
                    return;
                }
                grid_.son[idc] = new_ig;
                grid_.father[new_ig - 1] = idc;
                
                // Initialize child coordinates
                for (int d = 1; d <= NDIM; ++d) {
                    real_t dx_level = grid_.boxlen / std::pow(2.0, ilevel);
                    int ix = (ic - 1) & 1, iy = ((ic - 1) & 2) >> 1, iz = ((ic - 1) & 4) >> 2;
                    real_t offset = (d == 1) ? (ix - 0.5) : (d == 2) ? (iy - 0.5) : (iz - 0.5);
                    grid_.xg[(d - 1) * grid_.ngridmax + (new_ig - 1)] = grid_.xg[(d - 1) * grid_.ngridmax + (ig - 1)] + offset * dx_level;
                }

                // Gather neighbors for interpolation
                int icn[6]; grid_.get_nbor_cells(ign, ic, icn, ig);
                real_t u1[7][20] = {0};
                for(int iv=1; iv<=grid_.nvar; ++iv) u1[0][iv-1] = grid_.uold(idc, iv);
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

                for (int isc = 1; ic <= constants::twotondim; ++ic) {
                    int id_child = grid_.ncoarse + (isc - 1) * grid_.ngridmax + new_ig;
                    for (int iv = 1; iv <= grid_.nvar; ++iv) {
                        grid_.uold(id_child, iv) = u2[isc-1][iv-1];
                        grid_.unew(id_child, iv) = u2[isc-1][iv-1];
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
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            int ig_son = grid_.son[idc];
            if (ig_son > 0) {
                for (int iv = 1; iv <= grid_.nvar; ++iv) {
                    real_t sum = 0;
                    for (int is = 1; is <= constants::twotondim; ++is) {
                        sum += grid_.uold(grid_.ncoarse + (is - 1) * grid_.ngridmax + ig_son, iv);
                    }
                    grid_.uold(idc, iv) = sum / constants::twotondim;
                }
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

void TreeUpdater::flag_fine(int ilevel, real_t err_grad_d, real_t err_grad_p, real_t err_grad_v) {
    int myid = grid_.ncpu > 0 ? grid_.ncpu : 1;
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        int ign[7]; grid_.get_nbor_grids(igrid, ign);
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            grid_.flag1[idc] = 0;
            if (grid_.son[idc] != 0) continue;

            int icn[6]; grid_.get_nbor_cells(ign, ic, icn, igrid);
            real_t d_c = grid_.uold(idc, 1);
            for (int idim = 0; idim < NDIM; ++idim) {
                int idl = icn[idim * 2], idr = icn[idim * 2 + 1];
                if (idl > 0 && idr > 0) {
                    real_t d_l = grid_.uold(idl, 1), d_r = grid_.uold(idr, 1);
                    if (std::abs(d_r - d_l) / std::max({d_l, d_r, 1e-10}) > err_grad_d) {
                        grid_.flag1[idc] = 1;
                        break;
                    }
                }
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

} // namespace ramses
