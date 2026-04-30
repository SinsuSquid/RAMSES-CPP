#include "ramses/TreeUpdater.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <vector>

namespace ramses {

static inline void check(const std::vector<int>& v, size_t i, const std::string& name) {
    if (i >= v.size()) {
        std::cerr << "OUT OF BOUNDS: " << name << "[" << i << "] size=" << v.size() << std::endl;
        throw std::out_of_range(name);
    }
}

void TreeUpdater::refine_coarse() {
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
    if (grid_.numbf <= 0) return;
    int igrid = grid_.headf;
    check(grid_.next, igrid - 1, "next");
    grid_.headf = grid_.next[igrid - 1];
    if (grid_.headf > 0) {
        check(grid_.prev, grid_.headf - 1, "prev");
        grid_.prev[grid_.headf - 1] = 0;
    }
    grid_.numbf--;

    int nxny = params::nx * params::ny;
    int iz = (ind_cell - 1) / nxny;
    int iy = (ind_cell - 1 - iz * nxny) / params::nx;
    int ix = (ind_cell - 1 - iz * nxny - iy * params::nx);
    
    if (NDIM > 0) grid_.get_xg(igrid, 1) = static_cast<real_t>(ix) + 0.5f;
    if (NDIM > 1) grid_.get_xg(igrid, 2) = static_cast<real_t>(iy) + 0.5f;
    if (NDIM > 2) grid_.get_xg(igrid, 3) = static_cast<real_t>(iz) + 0.5f;

    check(grid_.son, ind_cell, "son");
    grid_.son[ind_cell] = igrid;
    check(grid_.father, igrid - 1, "father");
    grid_.father[igrid - 1] = ind_cell;

    for (int n = 1; n <= constants::twondim; ++n) {
        int idim = (n - 1) / 2; int side = (n - 1) % 2;
        int nix = ix, niy = iy, niz = iz;
        if (idim == 0) nix = (side == 0) ? (ix > 0 ? ix - 1 : params::nx - 1) : (ix < params::nx - 1 ? ix + 1 : 0);
        if (idim == 1) niy = (side == 0) ? (iy > 0 ? iy - 1 : params::ny - 1) : (iy < params::ny - 1 ? iy + 1 : 0);
        if (idim == 2) niz = (side == 0) ? (iz > 0 ? iz - 1 : params::nz - 1) : (iz < params::nz - 1 ? iz + 1 : 0);
        
        int n_idx = (n - 1) * grid_.ngridmax + (igrid - 1);
        check(grid_.nbor, n_idx, "nbor");
        grid_.nbor[n_idx] = 1 + nix + niy * params::nx + niz * nxny;
    }

    int icpu = grid_.cpu_map[ind_cell];
    for (int j = 1; j <= constants::twotondim; ++j) {
        int cell_idx = grid_.ncoarse + (j - 1) * grid_.ngridmax + igrid;
        check(grid_.cpu_map, cell_idx, "cpu_map");
        grid_.cpu_map[cell_idx] = icpu;
        // Copy physical data from father (level 0) to level 1 children
        for (int iv = 1; iv <= grid_.nvar; ++iv) {
            grid_.uold(cell_idx, iv) = grid_.uold(ind_cell, iv);
            grid_.unew(cell_idx, iv) = grid_.uold(ind_cell, iv);
        }
    }

    if (!boundary_region) {
        if (grid_.numbl(icpu, 1) > 0) {
            int tail = grid_.taill(icpu, 1);
            grid_.next[igrid - 1] = 0;
            grid_.prev[igrid - 1] = tail;
            grid_.next[tail - 1] = igrid;
            grid_.taill(icpu, 1) = igrid;
            grid_.numbl(icpu, 1)++;
        } else {
            grid_.next[igrid - 1] = 0; grid_.prev[igrid - 1] = 0;
            grid_.headl(icpu, 1) = igrid; grid_.taill(icpu, 1) = igrid;
            grid_.numbl(icpu, 1) = 1;
        }
    }
}

void TreeUpdater::refine_fine(int ilevel) {
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
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
}

void TreeUpdater::make_grid_fine(int ind_grid_father, int icell_pos, int ilevel, int ibound, bool boundary_region) {
    if (grid_.numbf <= 0) return;
    int igrid = grid_.headf;
    check(grid_.next, igrid - 1, "next");
    grid_.headf = grid_.next[igrid - 1];
    if (grid_.headf > 0) {
        check(grid_.prev, grid_.headf - 1, "prev");
        grid_.prev[grid_.headf - 1] = 0;
    }
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

    int igridn[7];
    grid_.get_nbor_grids(ind_grid_father, igridn);
    for (int n = 1; n <= constants::twondim; ++n) {
        int idim = (n - 1) / 2; int side = (n - 1) % 2;
        int ig = constants::iii[idim][side][icell_pos - 1];
        int ih = constants::jjj[idim][side][icell_pos - 1];
        if (ig == 0) grid_.get_nbor(igrid, n) = grid_.ncoarse + (ih - 1) * grid_.ngridmax + igrid;
        else if (igridn[ig] > 0) grid_.get_nbor(igrid, n) = grid_.ncoarse + (ih - 1) * grid_.ngridmax + igridn[ig];
        else grid_.get_nbor(igrid, n) = igridn[ig];
    }

    int icpu = grid_.cpu_map[ind_cell_father];
    for (int j = 1; j <= constants::twotondim; ++j) {
        int cell_idx = grid_.ncoarse + (j - 1) * grid_.ngridmax + igrid;
        check(grid_.cpu_map, cell_idx, "cpu_map");
        grid_.cpu_map[cell_idx] = icpu;
        // Copy physical data from father to children (piecewise constant)
        for (int iv = 1; iv <= grid_.nvar; ++iv) {
            grid_.uold(cell_idx, iv) = grid_.uold(ind_cell_father, iv);
            grid_.unew(cell_idx, iv) = grid_.uold(ind_cell_father, iv);
        }
    }

    int n_ilevel = ilevel + 1;
    if (!boundary_region) {
        if (grid_.numbl(icpu, n_ilevel) > 0) {
            int tail = grid_.taill(icpu, n_ilevel);
            grid_.next[igrid - 1] = 0;
            grid_.prev[igrid - 1] = tail;
            grid_.next[tail - 1] = igrid;
            grid_.taill(icpu, n_ilevel) = igrid;
            grid_.numbl(icpu, n_ilevel)++;
        } else {
            grid_.next[igrid - 1] = 0; grid_.prev[igrid - 1] = 0;
            grid_.headl(icpu, n_ilevel) = igrid; grid_.taill(icpu, n_ilevel) = igrid;
            grid_.numbl(icpu, n_ilevel) = 1;
        }
    }
}

void TreeUpdater::mark_cells(int ilevel) {
    const real_t gamma = 1.4;
    if (ilevel == 1) {
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            real_t d = std::max(grid_.uold(i, 1), 1e-10);
            real_t e_kin = 0.5 * (grid_.uold(i, 2)*grid_.uold(i, 2) + grid_.uold(i, 3)*grid_.uold(i, 3) + grid_.uold(i, 4)*grid_.uold(i, 4)) / d;
            real_t p = (grid_.uold(i, 5) - e_kin) * (gamma - 1.0);
            if (p > 1e-3) { grid_.flag1[i] = 1; grid_.flag2[i] = 1; }
            else { grid_.flag1[i] = 0; grid_.flag2[i] = 0; }
        }
        return;
    }
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel - 1);
        while (igrid > 0) {
            for (int ind = 1; ind <= constants::twotondim; ++ind) {
                int ind_cell = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;
                real_t d = std::max(grid_.uold(ind_cell, 1), 1e-10);
                real_t e_kin = 0.5 * (grid_.uold(ind_cell, 2)*grid_.uold(ind_cell, 2) + grid_.uold(ind_cell, 3)*grid_.uold(ind_cell, 3) + grid_.uold(ind_cell, 4)*grid_.uold(ind_cell, 4)) / d;
                real_t p = (grid_.uold(ind_cell, 5) - e_kin) * (gamma - 1.0);
                if (p > 1e-3) { grid_.flag1[ind_cell] = 1; grid_.flag2[ind_cell] = 1; }
                else { grid_.flag1[ind_cell] = 0; grid_.flag2[ind_cell] = 0; }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void TreeUpdater::mark_all(int ilevel) {
    if (ilevel == 1) {
        for (int i = 1; i <= grid_.ncoarse; ++i) { 
            check(grid_.flag1, i, "flag1");
            grid_.flag1[i] = 1; grid_.flag2[i] = 1; 
        }
        return;
    }
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel - 1);
        while (igrid > 0) {
            for (int ind = 1; ind <= constants::twotondim; ++ind) {
                int ind_cell = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;
                check(grid_.flag1, ind_cell, "flag1");
                grid_.flag1[ind_cell] = 1; grid_.flag2[ind_cell] = 1;
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

} // namespace ramses
