#include "ramses/AmrGrid.hpp"
#include "ramses/Parameters.hpp"

namespace ramses {

void AmrGrid::allocate(int nx_val, int ny_val, int nz_val, int ngridmax_val, int nvar_val, int ncpu_val, int nlevelmax_val) {
    ncoarse = nx_val * ny_val * nz_val;
    ngridmax = ngridmax_val;
    nvar = nvar_val;
    ncpu = ncpu_val;
    nlevelmax = nlevelmax_val;
    ncell = calculate_ncell(nx_val, ny_val, nz_val, ngridmax_val);

    // Allocate tree arrays (oct-based)
    xg.assign(ngridmax * NDIM, 0.0);
    father.assign(ngridmax, 0);
    nbor.assign(ngridmax * 2 * NDIM, 0);
    next.assign(ngridmax, 0);
    prev.assign(ngridmax, 0);

    // Allocate cell-based arrays
    son.assign(ncell + 1, 0);
    flag1.assign(ncell + 1, 0);
    flag2.assign(ncell + 1, 0);
    cpu_map.assign(ncell + 1, 0);

    // Allocate physical fields
    uold.allocate(ncell, nvar);
    unew.allocate(ncell, nvar);

    // Allocate linked list pointers
    headl.allocate(ncpu, nlevelmax);
    taill.allocate(ncpu, nlevelmax);
    numbl.allocate(ncpu, nlevelmax);

    headb.allocate(constants::MAXBOUND, nlevelmax);
    tailb.allocate(constants::MAXBOUND, nlevelmax);
    numbb.allocate(constants::MAXBOUND, nlevelmax);

    // Initialize free list
    for (int i = 1; i < ngridmax; ++i) {
        next[i - 1] = i + 1;
    }
    for (int i = 2; i <= ngridmax; ++i) {
        prev[i - 1] = i - 1;
    }
    
    headf = 1;
    tailf = ngridmax;
    numbf = ngridmax;
    
    prev[headf - 1] = 0;
    next[tailf - 1] = 0;
}

void AmrGrid::get_nbor_grids(int igrid, int igridn[7]) const {
    igridn[0] = igrid;
    for (int j = 1; j <= constants::twondim; ++j) {
        int nbor_cell = nbor[(j - 1) * ngridmax + (igrid - 1)];
        if (nbor_cell > 0) {
            igridn[j] = son[nbor_cell];
        } else {
            igridn[j] = 0;
        }
    }
}

void AmrGrid::get_nbor_cells(const int igridn[7], int icell_pos, int icelln[6]) const {
    for (int inbor = 0; inbor < 2; ++inbor) {
        for (int idim = 0; idim < NDIM; ++idim) {
            int nbor_idx = idim * 2 + inbor;
            int ig = constants::iii[idim][inbor][icell_pos - 1];
            int ih = constants::jjj[idim][inbor][icell_pos - 1];
            if (igridn[ig] > 0) {
                icelln[nbor_idx] = ncoarse + (ih - 1) * ngridmax + igridn[ig];
            } else {
                icelln[nbor_idx] = 0;
            }
        }
    }
}

void AmrGrid::get_3x3x3_father(int igrid, int nbors_father[27]) const {
    int ifather = father[igrid - 1];
    
    // Level 1 (Coarse) logic
    if (ifather <= ncoarse) {
        int nxny = ncpu > 1 ? 0 : (params::nx * params::ny); // Just use global params for now
        if (nxny == 0) nxny = params::nx * params::ny;

        int iz = (ifather - 1) / nxny;
        int iy = (ifather - 1 - iz * nxny) / params::nx;
        int ix = (ifather - 1 - iy * params::nx - iz * nxny);

        for (int k1 = 0; k1 < 3; ++k1) {
            int iiz = iz + k1 - 1;
            if (iiz < 0) iiz = params::nz - 1;
            if (iiz > params::nz - 1) iiz = 0;
            for (int j1 = 0; j1 < 3; ++j1) {
                int iiy = iy + j1 - 1;
                if (iiy < 0) iiy = params::ny - 1;
                if (iiy > params::ny - 1) iiy = 0;
                for (int i1 = 0; i1 < 3; ++i1) {
                    int iix = ix + i1 - 1;
                    if (iix < 0) iix = params::nx - 1;
                    if (iix > params::nx - 1) iix = 0;
                    
                    nbors_father[i1 + 3*j1 + 9*k1] = 1 + iix + iiy * params::nx + iiz * nxny;
                }
            }
        }
        return;
    }

    // Refined level logic (get3cubepos)
    int pos = (ifather - ncoarse - 1) / ngridmax + 1;
    int igrid_father = ifather - ncoarse - (pos - 1) * ngridmax;

    // kkk, jjj, iii from get3cubepos
    static const int kkk[8] = {5,5,5,5,6,6,6,6};
    static const int jjj[8] = {3,3,4,4,3,3,4,4};
    static const int iii[8] = {1,2,1,2,1,2,1,2};

    int nbors_grids[8]; // 2x2x2 cube of grids
    for (int kk = 0; kk < 2; ++kk) {
        int ig1 = igrid_father;
        if (kk > 0) {
            int n_cell = nbor[(kkk[pos-1] - 1) * ngridmax + (igrid_father - 1)];
            ig1 = (n_cell > 0) ? son[n_cell] : 0;
        }
        for (int jj = 0; jj < 2; ++jj) {
            int ig2 = ig1;
            if (jj > 0 && ig1 > 0) {
                int n_cell = nbor[(jjj[pos-1] - 1) * ngridmax + (ig1 - 1)];
                ig2 = (n_cell > 0) ? son[n_cell] : 0;
            }
            for (int ii = 0; ii < 2; ++ii) {
                int ig3 = ig2;
                if (ii > 0 && ig2 > 0) {
                    int n_cell = nbor[(iii[pos-1] - 1) * ngridmax + (ig2 - 1)];
                    ig3 = (n_cell > 0) ? son[n_cell] : 0;
                }
                nbors_grids[ii + 2*jj + 4*kk] = ig3;
            }
        }
    }

    // Map 2x2x2 grids to 3x3x3 cells using lll and mmm
    for (int j = 0; j < 27; ++j) {
        int ig = constants::lll[pos-1][j];
        int ic = constants::mmm[pos-1][j];
        int ig_idx = nbors_grids[ig - 1];
        if (ig_idx > 0) {
            nbors_father[j] = ncoarse + (ic - 1) * ngridmax + ig_idx;
        } else {
            nbors_father[j] = 0; // Should not happen in valid RAMSES tree
        }
    }
}




} // namespace ramses
