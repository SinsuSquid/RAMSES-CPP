#include "ramses/TreeUpdater.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/MpiManager.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

namespace ramses {

void TreeUpdater::refine_coarse() {
    for (int i = 1; i <= grid_.ncoarse; ++i) {
        if (grid_.flag1[i] == 1 && grid_.flag2[i] == 1 && grid_.son[i] == 0) {
            make_grid_coarse(i, 0, false);
        }
    }
}

void TreeUpdater::make_grid_coarse(int ind_cell, int ibound, bool boundary_region) {
    if (grid_.numbf <= 0) return;
    int igrid = grid_.headf;
    grid_.headf = grid_.next[igrid - 1];
    if (grid_.headf > 0) grid_.prev[grid_.headf - 1] = 0;
    grid_.numbf--;

    grid_.son[ind_cell] = igrid;
    grid_.father[igrid - 1] = ind_cell;

    int nxny = params::nx * params::ny;
    int iz = (ind_cell - 1) / nxny;
    int iy = (ind_cell - 1 - iz * nxny) / params::nx;
    int ix = (ind_cell - 1 - iy * params::nx - iz * nxny);

    if (NDIM > 0) grid_.get_xg(igrid, 1) = (static_cast<real_t>(ix) + 0.5f) / static_cast<real_t>(params::nx);
    if (NDIM > 1) grid_.get_xg(igrid, 2) = (static_cast<real_t>(iy) + 0.5f) / static_cast<real_t>(params::ny);
    if (NDIM > 2) grid_.get_xg(igrid, 3) = (static_cast<real_t>(iz) + 0.5f) / static_cast<real_t>(params::nz);

    for (int n = 1; n <= constants::twondim; ++n) grid_.get_nbor(igrid, n) = 0;

    int icpu_father = grid_.cpu_map[ind_cell];
    for (int j = 1; j <= constants::twotondim; ++j) {
        int cell_idx = grid_.ncoarse + (j - 1) * grid_.ngridmax + igrid;
        grid_.cpu_map[cell_idx] = icpu_father;
        for (int iv = 1; iv <= grid_.nvar; ++iv) {
            grid_.uold(cell_idx, iv) = grid_.uold(ind_cell, iv);
            grid_.unew(cell_idx, iv) = grid_.uold(ind_cell, iv);
        }
    }

    if (!boundary_region) {
        int myid = MpiManager::instance().rank() + 1;
        // In shared coarse refinement, all ranks must keep the grid in their list?
        // No, only the rank that OWNS the father cell should keep it.
        // But for Level 0-1, we want all ranks to have it.
        // Simplified: use icpu_father.
        if (grid_.numbl(icpu_father, 1) > 0) {
            int tail = grid_.taill(icpu_father, 1);
            grid_.next[tail - 1] = igrid;
            grid_.prev[igrid - 1] = tail;
            grid_.next[igrid - 1] = 0;
            grid_.taill(icpu_father, 1) = igrid;
        } else {
            grid_.headl(icpu_father, 1) = igrid;
            grid_.taill(icpu_father, 1) = igrid;
            grid_.prev[igrid - 1] = 0;
            grid_.next[igrid - 1] = 0;
        }
        grid_.numbl(icpu_father, 1)++;
    }
}

void TreeUpdater::refine_fine(int ilevel) {
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        int next_grid = grid_.next[igrid - 1];
        for (int ind = 1; ind <= constants::twotondim; ++ind) {
            int ind_cell = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;
            if (grid_.flag1[ind_cell] == 1 && grid_.flag2[ind_cell] == 1 && grid_.son[ind_cell] == 0) {
                make_grid_fine(igrid, ind, ilevel, 0, false);
            }
        }
        igrid = next_grid;
    }
}

void TreeUpdater::make_grid_fine(int ind_grid_father, int icell_pos, int ilevel, int ibound, bool boundary_region) {
    if (grid_.numbf <= 0) return;
    int igrid = grid_.headf;
    grid_.headf = grid_.next[igrid - 1];
    if (grid_.headf > 0) grid_.prev[grid_.headf - 1] = 0;
    grid_.numbf--;

    int ind_cell_father = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + ind_grid_father;
    grid_.son[ind_cell_father] = igrid;
    grid_.father[igrid - 1] = ind_cell_father;

    real_t scale = 0.5f / static_cast<real_t>(1 << (ilevel - 1));
    int ix = (icell_pos - 1) & 1;
    int iy = ((icell_pos - 1) & 2) >> 1;
    int iz = ((icell_pos - 1) & 4) >> 2;

    grid_.get_xg(igrid, 1) = grid_.get_xg(ind_grid_father, 1) + (static_cast<real_t>(ix) - 0.5f) * scale;
    if (NDIM > 1) grid_.get_xg(igrid, 2) = grid_.get_xg(ind_grid_father, 2) + (static_cast<real_t>(iy) - 0.5f) * scale;
    if (NDIM > 2) grid_.get_xg(igrid, 3) = grid_.get_xg(ind_grid_father, 3) + (static_cast<real_t>(iz) - 0.5f) * scale;

    int igridn[7]; grid_.get_nbor_grids(ind_grid_father, igridn);
    for (int n = 1; n <= 2*NDIM; ++n) {
        int idim = (n - 1) / 2; int inbor = (n - 1) % 2;
        int ig = constants::iii[idim][inbor][icell_pos - 1];
        int ih = constants::jjj[idim][inbor][icell_pos - 1];
        if (ig == 0) grid_.get_nbor(igrid, n) = grid_.ncoarse + (ih - 1) * grid_.ngridmax + ind_grid_father;
        else if (igridn[ig] > 0) grid_.get_nbor(igrid, n) = grid_.ncoarse + (ih - 1) * grid_.ngridmax + igridn[ig];
        else grid_.get_nbor(igrid, n) = igridn[ig];
    }

    int icpu = grid_.cpu_map[ind_cell_father];
    for (int j = 1; j <= constants::twotondim; ++j) {
        int cell_idx = grid_.ncoarse + (j - 1) * grid_.ngridmax + igrid;
        grid_.cpu_map[cell_idx] = icpu;
        for (int iv = 1; iv <= grid_.nvar; ++iv) {
            grid_.uold(cell_idx, iv) = grid_.uold(ind_cell_father, iv);
            grid_.unew(cell_idx, iv) = grid_.uold(ind_cell_father, iv);
        }
    }

    int n_ilevel = ilevel + 1;
    if (!boundary_region) {
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

void TreeUpdater::mark_cells(int ilevel) {
    real_t err_grad_d = grid_.err_grad_d, err_grad_p = grid_.err_grad_p, gamma = grid_.gamma;
    int myid = MpiManager::instance().rank() + 1;

    if (ilevel == 0) {
        for (int i = 1; i <= grid_.ncoarse; ++i) { grid_.flag1[i] = 1; grid_.flag2[i] = 1; }
        return;
    }

    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        int igridn[7]; grid_.get_nbor_grids(igrid, igridn);
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son[idc] != 0) continue;
            int icelln[6]; grid_.get_nbor_cells(igridn, ic, icelln);
            real_t d = grid_.uold(idc, 1), u = grid_.uold(idc, 2)/d, v = (NDIM>1?grid_.uold(idc,3)/d:0.0), w = (NDIM>2?grid_.uold(idc,4)/d:0.0);
            real_t p = std::max(d*1e-10, (grid_.uold(idc, 5) - 0.5*d*(u*u+v*v+w*w)) * (gamma - 1.0));
            bool refine = false;
            for (int in = 0; in < 2*NDIM; ++in) {
                if (icelln[in] > 0 && grid_.cpu_map[icelln[in]] == myid) {
                    real_t dn = grid_.uold(icelln[in], 1), un = grid_.uold(icelln[in], 2)/dn, vn = (NDIM>1?grid_.uold(icelln[in],3)/dn:0.0), wn = (NDIM>2?grid_.uold(icelln[in],4)/dn:0.0);
                    real_t pn = std::max(dn*1e-10, (grid_.uold(icelln[in], 5) - 0.5*dn*(un*un+vn*vn+wn*wn)) * (gamma - 1.0));
                    if (std::abs(d - dn) / std::max(d, dn) > err_grad_d) refine = true;
                    if (std::abs(p - pn) / std::max(p, pn) > err_grad_p) refine = true;
                }
            }
            if (refine) { grid_.flag1[idc] = 1; grid_.flag2[idc] = 1; }
            else { grid_.flag1[idc] = 0; grid_.flag2[idc] = 0; }
        }
        igrid = grid_.next[igrid - 1];
    }
}

void TreeUpdater::mark_all(int ilevel) {
    int myid = MpiManager::instance().rank() + 1;
    if (ilevel == 0) {
        for (int i = 1; i <= grid_.ncoarse; ++i) { grid_.flag1[i] = 1; grid_.flag2[i] = 1; }
        return;
    }
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ind = 1; ind <= constants::twotondim; ++ind) {
            int idc = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;
            grid_.flag1[idc] = 1; grid_.flag2[idc] = 1;
        }
        igrid = grid_.next[igrid - 1];
    }
}

} // namespace ramses
