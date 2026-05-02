#include "ramses/TreeUpdater.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/MpiManager.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace ramses {

void TreeUpdater::refine_coarse() {
    real_t dx = 1.0 / static_cast<real_t>(params::nx);
    for (int ic = 1; ic <= grid_.ncoarse; ++ic) {
        if (grid_.flag1[ic] == 1 && grid_.son[ic] == 0) {
            if (grid_.numbf <= 0) return;
            int igrid = grid_.headf;
            grid_.headf = grid_.next[igrid - 1];
            if (grid_.headf > 0) grid_.prev[grid_.headf - 1] = 0;
            else grid_.tailf = 0;
            grid_.numbf--;

            grid_.son[ic] = igrid;
            grid_.father[igrid - 1] = ic;

            int iz = (ic - 1) / (params::nx * params::ny);
            int iy = (ic - 1 - iz * params::nx * params::ny) / params::nx;
            int ix = (ic - 1 - iy * params::nx - iz * params::nx * params::ny);

            grid_.get_xg(igrid, 1) = (static_cast<real_t>(ix) + 0.5f) * dx;
            if (NDIM > 1) grid_.get_xg(igrid, 2) = (static_cast<real_t>(iy) + 0.5f) * dx;
            if (NDIM > 2) grid_.get_xg(igrid, 3) = (static_cast<real_t>(iz) + 0.5f) * dx;

            // Populate nbor grids
            for (int n = 1; n <= constants::twondim; ++n) {
                int dim = (n - 1) / 2;
                int side = (n - 1) % 2;
                real_t xn[3] = { grid_.get_xg(igrid, 1), grid_.get_xg(igrid, 2), grid_.get_xg(igrid, 3) };
                xn[dim] += (side == 0 ? -1.0 : 1.0) * dx;
                int inb = grid_.find_cell_by_coords(xn, 1, grid_.nboundary == 0);
                grid_.get_nbor(igrid, n) = inb;
            }

            int icpu = grid_.cpu_map[ic];
            for (int jc = 1; jc <= constants::twotondim; ++jc) {
                int cell_idx = grid_.ncoarse + (jc - 1) * grid_.ngridmax + igrid;
                grid_.cpu_map[cell_idx] = icpu;
                for (int iv = 1; iv <= grid_.nvar; ++iv) {
                    grid_.uold(cell_idx, iv) = grid_.uold(ic, iv);
                    grid_.unew(cell_idx, iv) = grid_.uold(ic, iv);
                }
                grid_.son[cell_idx] = 0;
            }

            int n_ilevel = 2;
            if (grid_.numbl(icpu, n_ilevel) > 0) {
                int tail = grid_.taill(icpu, n_ilevel);
                grid_.next[tail - 1] = igrid;
                grid_.prev[igrid - 1] = tail;
                grid_.next[igrid - 1] = 0;
                grid_.taill(icpu, n_ilevel) = igrid;
            } else {
                grid_.headl(icpu, n_ilevel) = igrid;
                grid_.taill(icpu, n_ilevel) = igrid;
                grid_.prev[igrid - 1] = 0;
                grid_.next[igrid - 1] = 0;
            }
            grid_.numbl(icpu, n_ilevel)++;
        }
    }
}

void TreeUpdater::refine_fine(int ilevel) {
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ind = 1; ind <= constants::twotondim; ++ind) {
            int idc = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;
            if (grid_.flag1[idc] == 1 && grid_.son[idc] == 0) {
                make_grid_fine(igrid, ind, ilevel, 0, false);
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

void TreeUpdater::make_grid_fine(int ind_grid_father, int icell_pos, int ilevel, int ibound, bool boundary_region) {
    if (grid_.numbf <= 0) return;
    int igrid = grid_.headf;
    grid_.headf = grid_.next[igrid - 1];
    if (grid_.headf > 0) grid_.prev[grid_.headf - 1] = 0;
    else grid_.tailf = 0;
    grid_.numbf--;

    int ind_cell_father = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + ind_grid_father;
    grid_.son[ind_cell_father] = igrid;
    grid_.father[igrid - 1] = ind_cell_father;

    real_t offset = 1.0 / static_cast<real_t>(params::nx * (1 << ilevel));
    int ix = (icell_pos - 1) & 1;
    int iy = ((icell_pos - 1) & 2) >> 1;
    int iz = ((icell_pos - 1) & 4) >> 2;

    grid_.get_xg(igrid, 1) = grid_.get_xg(ind_grid_father, 1) + (static_cast<real_t>(ix) - 0.5f) * 2.0f * offset;
    if (NDIM > 1) grid_.get_xg(igrid, 2) = grid_.get_xg(ind_grid_father, 2) + (static_cast<real_t>(iy) - 0.5f) * 2.0f * offset;
    if (NDIM > 2) grid_.get_xg(igrid, 3) = grid_.get_xg(ind_grid_father, 3) + (static_cast<real_t>(iz) - 0.5f) * 2.0f * offset;

    // Populate nbor grids
    real_t dx_father = 2.0f * offset;
    for (int n = 1; n <= constants::twondim; ++n) {
        int dim = (n - 1) / 2;
        int side = (n - 1) % 2;
        real_t xn[3] = { grid_.get_xg(igrid, 1), grid_.get_xg(igrid, 2), grid_.get_xg(igrid, 3) };
        xn[dim] += (side == 0 ? -1.0 : 1.0) * dx_father;
        int inb = grid_.find_cell_by_coords(xn, ilevel, grid_.nboundary == 0);
        grid_.get_nbor(igrid, n) = inb;
    }

    int icpu = grid_.cpu_map[ind_cell_father];
    for (int ic = 1; ic <= constants::twotondim; ++ic) {
        int cell_idx = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
        grid_.cpu_map[cell_idx] = icpu;
        for (int iv = 1; iv <= grid_.nvar; ++iv) {
            grid_.uold(cell_idx, iv) = grid_.uold(ind_cell_father, iv);
            grid_.unew(cell_idx, iv) = grid_.uold(ind_cell_father, iv);
        }
        grid_.son[cell_idx] = 0;
    }

    int n_ilevel = ilevel + 1;
    if (grid_.numbl(icpu, n_ilevel) > 0) {
        int tail = grid_.taill(icpu, n_ilevel);
        grid_.next[tail - 1] = igrid;
        grid_.prev[igrid - 1] = tail;
        grid_.next[igrid - 1] = 0;
        grid_.taill(icpu, n_ilevel) = igrid;
    } else {
        grid_.headl(icpu, n_ilevel) = igrid;
        grid_.taill(icpu, n_ilevel) = igrid;
        grid_.prev[igrid - 1] = 0;
        grid_.next[igrid - 1] = 0;
    }
    grid_.numbl(icpu, n_ilevel)++;
}

void TreeUpdater::mark_cells(int ilevel) {
    real_t err_grad_d = config_.get_double("refine_params", "err_grad_d", -1.0);
    real_t err_grad_p = config_.get_double("refine_params", "err_grad_p", -1.0);
    real_t err_grad_u = config_.get_double("refine_params", "err_grad_u", -1.0);
    real_t floor_d = 1e-10;

    int myid = MpiManager::instance().rank() + 1;
    
    if (ilevel == 1) {
        for (int idc = 1; idc <= grid_.ncoarse; ++idc) {
            bool refine = false;
            if (params::levelmin > 1) {
                refine = true;
            } else {
                real_t x[3]; grid_.get_cell_center(idc, x);
                for (int idim = 0; idim < NDIM; ++idim) {
                    real_t xp[3] = {x[0], x[1], x[2]};
                    real_t dx_c = 1.0 / static_cast<real_t>(params::nx);
                    xp[idim] -= dx_c; int id_l = grid_.find_cell_by_coords(xp, 1, grid_.nboundary == 0);
                    xp[idim] += 2.0 * dx_c; int id_r = grid_.find_cell_by_coords(xp, 1, grid_.nboundary == 0);
                    
                    if (id_l <= 0) id_l = idc; if (id_r <= 0) id_r = idc;

                    if (err_grad_d >= 0.0) {
                        real_t dl = grid_.uold(id_l, 1), dc = grid_.uold(idc, 1), dr = grid_.uold(id_r, 1);
                        real_t err = 2.0 * std::max(std::abs(dr - dc) / (dr + dc + floor_d), std::abs(dc - dl) / (dc + dl + floor_d));
                        if (err > err_grad_d) refine = true;
                    }
                    if (err_grad_u >= 0.0 && !refine) {
                        real_t ul = 0, uc = 0, ur = 0;
                        real_t dl = std::max(grid_.uold(id_l, 1), floor_d), dc = std::max(grid_.uold(idc, 1), floor_d), dr = std::max(grid_.uold(id_r, 1), floor_d);
                        for(int i=1; i<=NDIM; ++i) {
                            ul += (grid_.uold(id_l, 1+i)/dl)*(grid_.uold(id_l, 1+i)/dl);
                            uc += (grid_.uold(idc, 1+i)/dc)*(grid_.uold(idc, 1+i)/dc);
                            ur += (grid_.uold(id_r, 1+i)/dr)*(grid_.uold(id_r, 1+i)/dr);
                        }
                        ul = std::sqrt(ul); uc = std::sqrt(uc); ur = std::sqrt(ur);
                        real_t err = 2.0 * std::max(std::abs(ur - uc) / (ur + uc + floor_d), std::abs(uc - ul) / (uc + ul + floor_d));
                        if (err > err_grad_u) refine = true;
                    }
                    if (refine) break;
                }
            }
            grid_.flag1[idc] = refine ? 1 : 0;
        }
        return;
    }

    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        int igridn[7]; grid_.get_nbor_grids(igrid, igridn);
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            bool refine = false;
            
            if (ilevel < params::levelmin) {
                refine = true;
            } else {
                int icelln[6]; grid_.get_nbor_cells(igridn, ic, icelln, igrid);
                for (int idim = 0; idim < NDIM; ++idim) {
                    int id_l = icelln[idim * 2], id_r = icelln[idim * 2 + 1];
                    if (id_l <= 0) id_l = idc;
                    if (id_r <= 0) id_r = idc;

                    if (err_grad_d >= 0.0) {
                        real_t dl = grid_.uold(id_l, 1), dc = grid_.uold(idc, 1), dr = grid_.uold(id_r, 1);
                        real_t err = 2.0 * std::max(std::abs(dr - dc) / (dr + dc + floor_d), std::abs(dc - dl) / (dc + dl + floor_d));
                        if (err > err_grad_d) refine = true;
                    }
                    if (err_grad_p >= 0.0 && !refine) {
                        real_t pl = grid_.uold(id_l, 5), pc = grid_.uold(idc, 5), pr = grid_.uold(id_r, 5);
                        real_t err = 2.0 * std::max(std::abs(pr - pc) / (pr + pc + floor_d), std::abs(pc - pl) / (pc + pl + floor_d));
                        if (err > err_grad_p) refine = true;
                    }
                    if (err_grad_u >= 0.0 && !refine) {
                        real_t ul = 0, uc = 0, ur = 0;
                        real_t dl = std::max(grid_.uold(id_l, 1), floor_d), dc = std::max(grid_.uold(idc, 1), floor_d), dr = std::max(grid_.uold(id_r, 1), floor_d);
                        for(int i=1; i<=NDIM; ++i) {
                            ul += (grid_.uold(id_l, 1+i)/dl)*(grid_.uold(id_l, 1+i)/dl);
                            uc += (grid_.uold(idc, 1+i)/dc)*(grid_.uold(idc, 1+i)/dc);
                            ur += (grid_.uold(id_r, 1+i)/dr)*(grid_.uold(id_r, 1+i)/dr);
                        }
                        ul = std::sqrt(ul); uc = std::sqrt(uc); ur = std::sqrt(ur);
                        real_t err = 2.0 * std::max(std::abs(ur - uc) / (ur + uc + floor_d), std::abs(uc - ul) / (uc + ul + floor_d));
                        if (err > err_grad_u) refine = true;
                    }
                    if (refine) break;
                }
            }
            grid_.flag1[idc] = refine ? 1 : 0;
        }
        igrid = grid_.next[igrid - 1];
    }
}

void TreeUpdater::mark_all(int ilevel) {
    if (ilevel == 0) {
        for (int ic = 1; ic <= grid_.ncoarse; ++ic) grid_.flag1[ic] = 1;
        return;
    }
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            grid_.flag1[idc] = 1;
        }
        igrid = grid_.next[igrid - 1];
    }
}

} // namespace ramses
