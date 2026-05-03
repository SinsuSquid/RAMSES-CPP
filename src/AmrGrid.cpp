#include "ramses/AmrGrid.hpp"
#include "ramses/Constants.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace ramses {

void AmrGrid::allocate(int nx_val, int ny_val, int nz_val, int ngridmax_val, int nvar_val, int ncpu_val, int nlevelmax_val) {
    nx = std::max(nx_val, 1); ny = std::max(ny_val, 1); nz = std::max(nz_val, 1);
    ngridmax = ngridmax_val; nvar = nvar_val; ncpu = ncpu_val; nlevelmax = nlevelmax_val;
    ncoarse = nx * ny * nz;
    ncell = ncoarse + constants::twotondim * ngridmax;

    headl_vec.assign(ncpu * nlevelmax, 0);
    taill_vec.assign(ncpu * nlevelmax, 0);
    numbl_vec.assign(ncpu * nlevelmax, 0);

    next.assign(ngridmax, 0);
    prev.assign(ngridmax, 0);
    father.assign(ngridmax, 0);

    son.assign(ncell, 0);
    nbor.assign(ngridmax * 6, 0);
    flag1.assign(ncell, 0);
    cpu_map.assign(ncell, 0);

    xg.assign(NDIM * ngridmax, 0.0);
    uold_vec.assign(ncell * nvar, 0.0);
    unew_vec.assign(ncell * nvar, 0.0);
    
    phi.assign(ncell, 0.0);
    rho.assign(ncell, 0.0);
    f_vec.assign(ncell * 3, 0.0);

    // Initialize free list
    headf = 1; tailf = ngridmax; numbf = ngridmax;
    for (int i = 1; i < ngridmax; ++i) next[i - 1] = i + 1;
    for (int i = 2; i <= ngridmax; ++i) prev[i - 1] = i - 1;

    // Initialize coarse grid
    for (int i = 1; i <= ncoarse; ++i) {
        cpu_map[i - 1] = 1; // Default to CPU 1
    }
}

int AmrGrid::get_free_grid() {
    if (numbf == 0) return 0;
    int igrid = headf;
    headf = next[igrid - 1];
    if (headf > 0) prev[headf - 1] = 0;
    else tailf = 0;
    next[igrid - 1] = 0;
    numbf--;
    return igrid;
}

void AmrGrid::add_to_level_list(int igrid, int ilevel) {
    int myid = 1; // Simplification for serial
    int head = get_headl(myid, ilevel);
    if (head == 0) {
        headl_vec[(ilevel - 1) * ncpu + (myid - 1)] = igrid;
        taill_vec[(ilevel - 1) * ncpu + (myid - 1)] = igrid;
        prev[igrid - 1] = 0;
        next[igrid - 1] = 0;
    } else {
        int tail = taill(myid, ilevel);
        next[tail - 1] = igrid;
        prev[igrid - 1] = tail;
        next[igrid - 1] = 0;
        taill_vec[(ilevel - 1) * ncpu + (myid - 1)] = igrid;
    }
    numbl_vec[(ilevel - 1) * ncpu + (myid - 1)]++;
}

int AmrGrid::count_grids_at_level(int ilevel) const {
    int total = 0;
    for (int i = 1; i <= ncpu; ++i) total += numbl(i, ilevel);
    return total;
}

void AmrGrid::get_nbor_grids(int igrid, int ign[7]) const {
    ign[0] = igrid;
    for (int i = 1; i <= 6; ++i) {
        int ic_nbor = nbor[(i - 1) * ngridmax + igrid - 1];
        if (ic_nbor > 0) ign[i] = son[ic_nbor - 1];
        else ign[i] = 0;
    }
}

void AmrGrid::get_nbor_cells(const int ign[7], int ic, int icn[6], int igrid) const {
    for (int i = 0; i < 6; ++i) icn[i] = 0;
    for (int idim = 0; idim < NDIM; ++idim) {
        for (int inbor = 0; inbor < 2; ++inbor) {
            int ig_idx = constants::iii[idim][inbor][ic - 1];
            int ic_pos = constants::jjj[idim][inbor][ic - 1];
            int ig = ign[ig_idx];
            if (ig > 0) {
                icn[idim * 2 + inbor] = ncoarse + (ic_pos - 1) * ngridmax + ig;
            } else {
                icn[idim * 2 + inbor] = nbor[(idim * 2 + inbor) * ngridmax + igrid - 1];
            }
        }
    }
}

} // namespace ramses
