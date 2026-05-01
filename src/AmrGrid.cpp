#include "ramses/AmrGrid.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace ramses {

void AmrGrid::allocate(int nx_val, int ny_val, int nz_val, int ngridmax_val, int nvar_val, int ncpu_val, int nlevelmax_val) {
    ncoarse = nx_val * ny_val * nz_val;
    ngridmax = ngridmax_val;
    nvar = nvar_val;
    ncpu = ncpu_val;
    nlevelmax = nlevelmax_val;
    ndim = NDIM;
    ncell = ncoarse + constants::twotondim * ngridmax;

    xg.assign(static_cast<size_t>(ngridmax) * 3, 0.0);
    father.assign(ngridmax, 0);
    nbor.assign(static_cast<size_t>(ngridmax) * 6, 0);
    next.assign(ngridmax, 0);
    prev.assign(ngridmax, 0);

    son.assign(ncell, 0);
    flag1.assign(ncell, 0);
    flag2.assign(ncell, 0);
    cpu_map.assign(ncell, 1);
    hilbert_keys.assign(ncell, 0);
    
    uold.allocate(ncell, nvar);
    unew.allocate(ncell, nvar);
    phi.assign(ncell, 0.0);
    f.allocate(ncell, 3);

    headl.allocate(ncpu, nlevelmax);
    taill.allocate(ncpu, nlevelmax);
    numbl.allocate(ncpu, nlevelmax);

    headf = 1; tailf = ngridmax; numbf = ngridmax;
    for (int i = 1; i <= ngridmax; ++i) {
        next[i - 1] = (i < ngridmax) ? (i + 1) : 0;
        prev[i - 1] = (i > 1) ? (i - 1) : 0;
    }

    std::cout << "[AmrGrid] Allocated: ncoarse=" << ncoarse << " ngridmax=" << ngridmax << " ncell=" << ncell << std::endl;
}

void AmrGrid::get_nbor_grids(int igrid, int igridn[7]) const {
    igridn[0] = igrid;
    if (igrid <= 0 || igrid > ngridmax) return;
    for (int i = 1; i <= 6; ++i) {
        int inb = nbor[(i - 1) * ngridmax + (igrid - 1)];
        igridn[i] = (inb > 0) ? son[inb] : 0;
    }
}

void AmrGrid::get_nbor_cells(const int igridn[7], int icell_pos, int icelln[6], int igrid) const {
    for (int i = 0; i < 2 * NDIM; ++i) {
        int idim = i / 2;
        int inbor = i % 2;
        int ig = constants::iii[idim][inbor][icell_pos - 1];
        int ih = constants::jjj[idim][inbor][icell_pos - 1];
        int ig_idx_n = igridn[ig];
        if (ig_idx_n > 0) {
            icelln[i] = ncoarse + (ih - 1) * ngridmax + ig_idx_n;
        } else {
            if (ig == 0) {
                if (igrid > 0) {
                    icelln[i] = ncoarse + (ih - 1) * ngridmax + igrid;
                } else {
                    // Coarse cell neighbor lookup
                    real_t x[3]; get_cell_center(icell_pos, x);
                    x[idim] += (inbor == 0 ? -1.0 : 1.0) / static_cast<real_t>(params::nx);
                    icelln[i] = find_cell_by_coords(x, 1);
                }
            } else {
                if (igrid > 0) {
                    icelln[i] = nbor[(ig - 1) * ngridmax + (igrid - 1)];
                } else {
                    // Coarse cell neighbor lookup
                    real_t x[3]; get_cell_center(icell_pos, x);
                    x[idim] += (inbor == 0 ? -1.0 : 1.0) / static_cast<real_t>(params::nx);
                    icelln[i] = find_cell_by_coords(x, 1);
                }
            }
        }
    }
}

void AmrGrid::setup_root_periodicity() {
    if (ncoarse == 1) {
        // Coarse grid 1 is the only grid. Neighbors are cell 1 itself.
        // Wait, nbor array is for GRIDS, not cell-index in Fortran?
        // Actually nbor array in Fortran stores cell indices.
        // But since we are level 1, there's no grid index, just cell indices.
        // Wait, nbor is only for grids? No.
    }
}

void AmrGrid::get_27_cell_neighbors(int icell, int nbors[27]) const {
    for (int i = 0; i < 27; ++i) nbors[i] = 0;
    real_t x[3]; get_cell_center(icell, x);
    
    int ilevel = 1;
    if (icell > ncoarse) {
        int igrid = ((icell - ncoarse - 1) % ngridmax) + 1;
        int ifath = father[igrid - 1];
        ilevel = 2;
        while (ifath > ncoarse) {
            ilevel++;
            int ig_f = ((ifath - ncoarse - 1) % ngridmax) + 1;
            ifath = father[ig_f - 1];
            if (ilevel > 50) break;
        }
    }
    real_t dx = 1.0 / static_cast<real_t>(params::nx * (1 << (ilevel - 1)));

    for (int iz = (NDIM > 2 ? -1 : 0); iz <= (NDIM > 2 ? 1 : 0); ++iz) {
        for (int iy = (NDIM > 1 ? -1 : 0); iy <= (NDIM > 1 ? 1 : 0); ++iy) {
            for (int ix = -1; ix <= 1; ++ix) {
                real_t xn[3] = { x[0] + ix * dx, x[1] + iy * dx, x[2] + iz * dx };
                for(int d=0; d<3; ++d) { 
                    while (xn[d] < 0.0) xn[d] += 1.0; 
                    while (xn[d] >= 1.0) xn[d] -= 1.0; 
                }
                int idx = (iz + 1) * 9 + (iy + 1) * 3 + (ix + 1);
                nbors[idx] = find_cell_by_coords(xn, ilevel);
            }
        }
    }
}

int AmrGrid::find_cell_by_coords(const real_t x[3], int ilevel_max) const {
    if (ilevel_max < 0) ilevel_max = nlevelmax;
    
    real_t xp[3] = {x[0], x[1], x[2]};
    for(int d=0; d<3; ++d) {
        while(xp[d] < 0.0) xp[d] += 1.0;
        while(xp[d] >= 1.0) xp[d] -= 1.0;
    }

    int nx = params::nx, ny = params::ny, nz = params::nz;
    int ix = std::min(nx - 1, std::max(0, static_cast<int>(std::floor(xp[0] * nx))));
    int iy = std::min(ny - 1, std::max(0, static_cast<int>(std::floor(xp[1] * ny))));
    int iz = std::min(nz - 1, std::max(0, static_cast<int>(std::floor(xp[2] * nz))));
    int icell = 1 + ix + iy * nx + iz * nx * ny;

    int il = 1;
    while (son[icell] > 0 && il < ilevel_max) {
        int igrid = son[icell];
        int cx = (xp[0] > get_xg(igrid, 1)) ? 1 : 0;
        int cy = (NDIM > 1 && xp[1] > get_xg(igrid, 2)) ? 1 : 0;
        int cz = (NDIM > 2 && xp[2] > get_xg(igrid, 3)) ? 1 : 0;
        int ic = 1 + cx + cy * 2 + cz * 4;
        icell = ncoarse + (ic - 1) * ngridmax + igrid;
        il++;
    }
    return icell;
}

void AmrGrid::get_cell_center(int icell, real_t x[3]) const {
    x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
    if (icell <= ncoarse) {
        int nx = params::nx, ny = params::ny, nz = params::nz;
        int iz = (icell - 1) / (nx * ny);
        int iy = (icell - 1 - iz * nx * ny) / nx;
        int ix = (icell - 1 - iy * nx - iz * nx * ny);
        real_t dx = 1.0 / static_cast<real_t>(nx);
        x[0] = (static_cast<real_t>(ix) + 0.5f) * dx;
        if (NDIM > 1) x[1] = (static_cast<real_t>(iy) + 0.5f) * dx;
        if (NDIM > 2) x[2] = (static_cast<real_t>(iz) + 0.5f) * dx;
    } else {
        int igrid = ((icell - ncoarse - 1) % ngridmax) + 1;
        int ic = (icell - ncoarse - 1) / ngridmax;
        int cz = ic / 4; int cy = (ic % 4) / 2; int cx = ic % 2;
        
        int ifath = father[igrid - 1];
        int ilevel = 2;
        while (ifath > ncoarse) {
            ilevel++;
            int ig_f = ((ifath - ncoarse - 1) % ngridmax) + 1;
            ifath = father[ig_f - 1];
            if (ilevel > 50) break;
        }
        real_t dx = 1.0 / static_cast<real_t>(params::nx * (1 << (ilevel - 1)));
        x[0] = get_xg(igrid, 1) + (static_cast<real_t>(cx) - 0.5f) * dx;
        if (NDIM > 1) x[1] = get_xg(igrid, 2) + (static_cast<real_t>(cy) - 0.5f) * dx;
        if (NDIM > 2) x[2] = get_xg(igrid, 3) + (static_cast<real_t>(cz) - 0.5f) * dx;
    }
}

void AmrGrid::restrict_coarse() {
    for (int ic = 1; ic <= ncoarse; ++ic) {
        if (son[ic] > 0) {
            int ison = son[ic];
            for (int iv = 1; iv <= nvar; ++iv) {
                real_t sum = 0;
                for (int jc = 1; jc <= constants::twotondim; ++jc) {
                    int idc = ncoarse + (jc - 1) * ngridmax + ison;
                    sum += uold(idc, iv);
                }
                uold(ic, iv) = sum / static_cast<real_t>(constants::twotondim);
                unew(ic, iv) = uold(ic, iv);
            }
        }
    }
}

void AmrGrid::restrict_fine(int ilevel) {
    if (ilevel < 2 || ilevel > nlevelmax) return;
    int myid = MpiManager::instance().rank() + 1;
    int igrid = get_headl(myid, ilevel - 1);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int id_father = ncoarse + (ic - 1) * ngridmax + igrid;
            if (son[id_father] > 0) {
                int igrid_son = son[id_father];
                for (int iv = 1; iv <= nvar; ++iv) {
                    real_t sum = 0;
                    for (int jc = 1; jc <= constants::twotondim; ++jc) {
                        int id_son = ncoarse + (jc - 1) * ngridmax + igrid_son;
                        sum += uold(id_son, iv);
                    }
                    uold(id_father, iv) = sum / static_cast<real_t>(constants::twotondim);
                    unew(id_father, iv) = uold(id_father, iv);
                }
            }
        }
        igrid = next[igrid - 1];
    }
}

} // namespace ramses
