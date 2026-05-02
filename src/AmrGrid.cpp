#include "ramses/AmrGrid.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace ramses {

void AmrGrid::allocate(int nx, int ny, int nz, int ngridmax_in, int nvar_in, int ncpu_in, int nlevelmax_in) {
    ngridmax = ngridmax_in;
    nvar = nvar_in;
    ncpu = ncpu_in;
    nlevelmax = nlevelmax_in;
    ndim = NDIM;
    ncoarse = nx * ny * nz;
    ncell = calculate_ncell(nx, ny, nz, ngridmax);

    xg.assign(ngridmax * 3, 0.0);
    father.assign(ngridmax, 0);
    nbor.assign(ngridmax * 6, 0);
    next.assign(ngridmax, 0);
    prev.assign(ngridmax, 0);

    son.assign(ncell + 1, 0);
    flag1.assign(ncell + 1, 0);
    flag2.assign(ncell + 1, 0);
    cpu_map.assign(ncell + 1, 0);
    hilbert_keys.assign(ncell + 1, 0.0);

    uold.allocate(ncell, nvar);
    unew.allocate(ncell, nvar);
    divu.assign(ncell + 1, 0.0);
    phi.assign(ncell + 1, 0.0);
    f.allocate(ncell, 3);
    rho.assign(ncell + 1, 0.0);

    headl.allocate(ncpu, nlevelmax);
    taill.allocate(ncpu, nlevelmax);
    numbl.allocate(ncpu, nlevelmax);

    headf = 1;
    tailf = ngridmax;
    numbf = ngridmax;
    for (int i = 1; i <= ngridmax; ++i) {
        next[i - 1] = (i < ngridmax) ? i + 1 : 0;
        prev[i - 1] = (i > 1) ? i - 1 : 0;
    }
}

void AmrGrid::get_nbor_grids(int igrid, int igridn[7]) const {
    igridn[0] = igrid;
    if (igrid <= 0 || igrid > ngridmax) return;
    for (int i = 1; i <= 2 * NDIM; ++i) {
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
                icelln[i] = ncoarse + (ih - 1) * ngridmax + igrid;
            } else {
                icelln[i] = nbor[(ig - 1) * ngridmax + (igrid - 1)];
            }
        }
    }
}

void AmrGrid::get_grids_of_nbor_cells(int igrid, int icell_pos, int nbors_grids[8]) const {
    int iii_loc[8] = {1,2,1,2,1,2,1,2};
    int jjj_loc[8] = {3,3,4,4,3,3,4,4};
    int kkk_loc[8] = {5,5,5,5,6,6,6,6};

    for(int i=0; i<8; ++i) nbors_grids[i] = 0;
    
    int iimin=0, iimax=1, jjmin=0, jjmax=0, kkmin=0, kkmax=0;
    if(NDIM > 1) jjmax=1;
    if(NDIM > 2) kkmax=1;

    for(int kk=kkmin; kk<=kkmax; ++kk) {
        int ig1 = igrid;
        if(kk > 0) ig1 = (igrid > 0) ? son[get_nbor(igrid, kkk_loc[icell_pos-1])] : 0;
        for(int jj=jjmin; jj<=jjmax; ++jj) {
            int ig2 = ig1;
            if(jj > 0) ig2 = (ig1 > 0) ? son[get_nbor(ig1, jjj_loc[icell_pos-1])] : 0;
            for(int ii=iimin; ii<=iimax; ++ii) {
                int ig3 = ig2;
                if(ii > 0) ig3 = (ig2 > 0) ? son[get_nbor(ig2, iii_loc[icell_pos-1])] : 0;
                int inbor = 1 + ii + 2*jj + 4*kk;
                nbors_grids[inbor-1] = ig3;
            }
        }
    }
}

void AmrGrid::get_27_cell_neighbors(int icell, int nbors[27]) const {
    for (int i = 0; i < 27; ++i) nbors[i] = 0;
    if (icell <= ncoarse) {
        real_t x[3]; get_cell_center(icell, x);
        int nx = params::nx, ny = params::ny, nz = params::nz;
        real_t dx = 1.0 / static_cast<real_t>(nx);
        for (int iz = (NDIM > 2 ? -1 : 0); iz <= (NDIM > 2 ? 1 : 0); ++iz) {
            for (int iy = (NDIM > 1 ? -1 : 0); iy <= (NDIM > 1 ? 1 : 0); ++iy) {
                for (int ix = -1; ix <= 1; ++ix) {
                    real_t xn[3] = { x[0] + ix * dx, x[1] + iy * dx, x[2] + iz * dx };
                    int idx = (iz + 1) * 9 + (iy + 1) * 3 + (ix + 1);
                    nbors[idx] = find_cell_by_coords(xn, 1, nboundary == 0);
                }
            }
        }
        return;
    }

    int igrid = ((icell - ncoarse - 1) % ngridmax) + 1;
    int icell_pos = (icell - ncoarse - 1) / ngridmax + 1;
    
    int nbors_grids[8];
    get_grids_of_nbor_cells(igrid, icell_pos, nbors_grids);

    for (int j = 1; j <= constants::threetondim; ++j) {
        int ig = constants::lll[icell_pos - 1][j - 1];
        int ic = constants::mmm[icell_pos - 1][j - 1];
        int ig_idx = nbors_grids[ig - 1];
        if (ig_idx > 0) {
            nbors[j - 1] = ncoarse + (ic - 1) * ngridmax + ig_idx;
        } else {
            // Coarse neighbor fallback
            // In RAMSES, this case is handled by get3cubefather.
            // For now, use coordinate lookup as a slow fallback.
            real_t x[3], xn[3]; get_cell_center(icell, x);
            // This is still slow but happens only at level boundaries.
            // Wait, we need to know the offset for j.
            int ix = (j-1)%3 - 1; int iy = ((j-1)/3)%3 - 1; int iz = (j-1)/9 - 1;
            // ... actually, let's just use coordinate lookup for all 27 if any grid is missing.
            // NO! Let's be smart.
        }
    }
    // Revert to coordinate lookup for the whole cube if any neighbor grid is missing
    // to ensure correctness until get3cubefather is fully ported.
    real_t x[3]; get_cell_center(icell, x);
    int ifath = father[igrid - 1];
    int ilevel = 2;
    while (ifath > ncoarse) {
        ilevel++;
        int ig_f = ((ifath - ncoarse - 1) % ngridmax) + 1;
        ifath = father[ig_f - 1];
    }
    real_t dx = 1.0 / static_cast<real_t>(params::nx * (1 << (ilevel - 1)));
    for (int iz = (NDIM > 2 ? -1 : 0); iz <= (NDIM > 2 ? 1 : 0); ++iz) {
        for (int iy = (NDIM > 1 ? -1 : 0); iy <= (NDIM > 1 ? 1 : 0); ++iy) {
            for (int ix = -1; ix <= 1; ++ix) {
                real_t xn[3] = { x[0] + ix * dx, x[1] + iy * dx, x[2] + iz * dx };
                int idx = (iz + 1) * 9 + (iy + 1) * 3 + (ix + 1);
                nbors[idx] = find_cell_by_coords(xn, ilevel, nboundary == 0);
            }
        }
    }
}

void AmrGrid::setup_root_periodicity() {}

int AmrGrid::find_cell_by_coords(const real_t x[3], int ilevel_max, bool periodic) const {
    if (ilevel_max < 0) ilevel_max = nlevelmax;
    
    real_t xp[3] = {x[0], x[1], x[2]};
    for(int d=0; d<3; ++d) {
        if (periodic) {
            while(xp[d] < 0.0) xp[d] += 1.0;
            while(xp[d] >= 1.0) xp[d] -= 1.0;
        } else {
            if (xp[d] < 0.0 || xp[d] >= 1.0) return 0;
        }
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
