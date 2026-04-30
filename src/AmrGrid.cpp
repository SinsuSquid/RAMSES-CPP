#include "ramses/AmrGrid.hpp"
#include "ramses/Parameters.hpp"
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

    std::cout << "[AmrGrid] Allocating Grid: ncoarse=" << ncoarse 
              << " ngridmax=" << ngridmax 
              << " ncell=" << ncell << " ncpu=" << ncpu << std::endl;

    xg.assign(static_cast<size_t>(ngridmax) * NDIM, 0.0);
    father.assign(static_cast<size_t>(ngridmax), 0);
    nbor.assign(static_cast<size_t>(ngridmax) * 2 * NDIM, 0);
    next.assign(static_cast<size_t>(ngridmax), 0);
    prev.assign(static_cast<size_t>(ngridmax), 0);

    son.assign(static_cast<size_t>(ncell) + 1, 0);
    flag1.assign(static_cast<size_t>(ncell) + 1, 0);
    flag2.assign(static_cast<size_t>(ncell) + 1, 0);
    cpu_map.assign(static_cast<size_t>(ncell) + 1, 1);
    hilbert_keys.assign(static_cast<size_t>(ncell) + 1, 0.0);
    divu.assign(static_cast<size_t>(ncell) + 1, 0.0);
    phi.assign(static_cast<size_t>(ncell) + 1, 0.0);
    f.allocate(ncell, NDIM);
    rho.assign(static_cast<size_t>(ncell) + 1, 0.0);

    uold.allocate(ncell, nvar);
    unew.allocate(ncell, nvar);

    headl.allocate(ncpu, nlevelmax);
    taill.allocate(ncpu, nlevelmax);
    numbl.allocate(ncpu, nlevelmax);

    headb.allocate(constants::MAXBOUND, nlevelmax);
    tailb.allocate(constants::MAXBOUND, nlevelmax);
    numbb.allocate(constants::MAXBOUND, nlevelmax);

    for (int i = 1; i < ngridmax; ++i) next[i - 1] = i + 1;
    next[ngridmax - 1] = 0;
    headf = 1;
    tailf = ngridmax;
    numbf = ngridmax;

    for (int i = 2; i <= ngridmax; ++i) prev[i - 1] = i - 1;
    prev[0] = 0;
}

void AmrGrid::get_nbor_grids(int igrid, int igridn[7]) const {
    igridn[0] = igrid;
    for (int j = 1; j <= 2 * NDIM; ++j) {
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

void AmrGrid::get_cell_center(int icell, real_t x[3]) const {
    if (icell <= ncoarse) {
        int nxny = params::nx * params::ny;
        int iz = (icell - 1) / nxny;
        int iy = (icell - 1 - iz * nxny) / params::nx;
        int ix = (icell - 1 - iy * params::nx - iz * nxny);

        x[0] = (static_cast<real_t>(ix) + 0.5f) / static_cast<real_t>(params::nx);
        if (ndim > 1) x[1] = (static_cast<real_t>(iy) + 0.5f) / static_cast<real_t>(params::ny);
        else x[1] = 0.5f;
        if (ndim > 2) x[2] = (static_cast<real_t>(iz) + 0.5f) / static_cast<real_t>(params::nz);
        else x[2] = 0.5f;
        return;
    }

    int pos = (icell - ncoarse - 1) / ngridmax + 1;
    int igrid = icell - ncoarse - (pos - 1) * ngridmax;
    
    int ifather = father[igrid - 1];
    int ilevel = 2;
    int tmp_if = ifather;
    while (tmp_if > ncoarse) {
        ilevel++;
        int igp = (tmp_if - ncoarse - 1) / ngridmax + 1;
        int ig_idx = tmp_if - ncoarse - (igp - 1) * ngridmax;
        tmp_if = father[ig_idx - 1];
    }

    real_t dx = 1.0f / static_cast<real_t>(1 << (ilevel - 1));

    for (int idim = 0; idim < ndim; ++idim) {
        x[idim] = xg[idim * ngridmax + igrid - 1];
    }
    for (int idim = ndim; idim < 3; ++idim) x[idim] = 0.5f;

    int ix = (pos - 1) & 1;
    int iy = ((pos - 1) & 2) >> 1;
    int iz = ((pos - 1) & 4) >> 2;

    x[0] += (static_cast<real_t>(ix) - 0.5f) * dx;
    if (ndim > 1) x[1] += (static_cast<real_t>(iy) - 0.5f) * dx;
    if (ndim > 2) x[2] += (static_cast<real_t>(iz) - 0.5f) * dx;
}

int AmrGrid::find_cell_by_coords(const real_t x[3], int ilevel_max) const {
    int nxny = params::nx * params::ny;
    int ix = std::min(params::nx - 1, std::max(0, static_cast<int>(std::floor(x[0] * params::nx))));
    int iy = (ndim > 1) ? std::min(params::ny - 1, std::max(0, static_cast<int>(std::floor(x[1] * params::ny)))) : 0;
    int iz = (ndim > 2) ? std::min(params::nz - 1, std::max(0, static_cast<int>(std::floor(x[2] * params::nz)))) : 0;
    int icell = 1 + ix + iy * params::nx + iz * nxny;

    int ilevel = 1;
    while (icell > 0 && icell <= ncell && son[icell] > 0 && (ilevel_max < 0 || ilevel < ilevel_max)) {
        int igrid = son[icell];
        real_t xg0 = get_xg(igrid, 1), xg1 = (ndim > 1 ? get_xg(igrid, 2) : 0.5), xg2 = (ndim > 2 ? get_xg(igrid, 3) : 0.5);
        int iix = (x[0] > xg0) ? 1 : 0;
        int iiy = (ndim > 1 && x[1] > xg1) ? 1 : 0;
        int iiz = (ndim > 2 && x[2] > xg2) ? 1 : 0;
        int pos = 1 + iix + 2 * iiy + 4 * iiz;
        icell = ncoarse + (pos - 1) * ngridmax + igrid;
        ilevel++;
    }
    if (icell > ncell || icell <= 0) {
        std::cerr << "[AmrGrid] Error: find_cell_by_coords out of bounds! icell=" << icell << " ncell=" << ncell << " ilevel=" << ilevel << std::endl;
    }
    return icell;
}

void AmrGrid::get_3x3x3_father(int igrid, int nbors_father[27]) const {
    int ifather = father[igrid - 1];
    if (ifather <= ncoarse) {
        int nxny = params::nx * params::ny;
        int iz = (ifather - 1) / nxny;
        int iy = (ifather - 1 - iz * nxny) / params::nx;
        int ix = (ifather - 1 - iy * params::nx - iz * nxny);
        for (int k1 = 0; k1 < 3; ++k1) for (int j1 = 0; j1 < 3; ++j1) for (int i1 = 0; i1 < 3; ++i1) {
            int iix = ix + i1 - 1, iiy = iy + j1 - 1, iiz = iz + k1 - 1;
            if (iix < 0) iix = params::nx - 1; if (iix >= params::nx) iix = 0;
            if (iiy < 0) iiy = params::ny - 1; if (iiy >= params::ny) iiy = 0;
            if (iiz < 0) iiz = params::nz - 1; if (iiz >= params::nz) iiz = 0;
            nbors_father[i1 + 3*j1 + 9*k1] = 1 + iix + iiy * params::nx + iiz * nxny;
        }
        return;
    }

    int pos = (ifather - ncoarse - 1) / ngridmax + 1;
    int ig_idx = ifather - ncoarse - (pos - 1) * ngridmax;
    static const int kkk[8] = {5,5,5,5,6,6,6,6}, jjj[8] = {3,3,4,4,3,3,4,4}, iii[8] = {1,2,1,2,1,2,1,2};
    int nbors_grids[8] = {0}; 
    for(int kk=0; kk<2; ++kk) for(int jj=0; jj<2; ++jj) for(int ii=0; ii<2; ++ii) {
        int ig1 = ig_idx;
        if (kk == 1) { int n_cell = nbor[(kkk[pos-1] - 1) * ngridmax + (ig1 - 1)]; ig1 = (n_cell > 0) ? son[n_cell] : 0; }
        if (ig1 > 0) {
            int ig2 = ig1;
            if (jj == 1) { int n_cell = nbor[(jjj[pos-1] - 1) * ngridmax + (ig2 - 1)]; ig2 = (n_cell > 0) ? son[n_cell] : 0; }
            if (ig2 > 0) {
                int ig3 = ig2;
                if (ii == 1) { int n_cell = nbor[(iii[pos-1] - 1) * ngridmax + (ig3 - 1)]; ig3 = (n_cell > 0) ? son[n_cell] : 0; }
                nbors_grids[ii + 2*jj + 4*kk] = ig3;
            }
        }
    }
    for (int j = 0; j < 27; ++j) {
        int ig = constants::lll[pos-1][j], ic = constants::mmm[pos-1][j];
        int ig_idx_n = nbors_grids[ig - 1];
        nbors_father[j] = (ig_idx_n > 0) ? (ncoarse + (ic - 1) * ngridmax + ig_idx_n) : 0;
    }
}

} // namespace ramses
