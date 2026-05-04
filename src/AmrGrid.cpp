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

void AmrGrid::free_grid(int igrid) {
    if (igrid <= 0) return;
    if (tailf == 0) {
        headf = igrid;
        tailf = igrid;
        prev[igrid - 1] = 0;
        next[igrid - 1] = 0;
    } else {
        next[tailf - 1] = igrid;
        prev[igrid - 1] = tailf;
        next[igrid - 1] = 0;
        tailf = igrid;
    }
    numbf++;
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

void AmrGrid::get_27_cell_neighbors(int icell, int nbors[27]) const {
    // 1. Identify current cell level
    int ilevel = 0;
    if (icell > ncoarse) {
        int curr_ig = ((icell - ncoarse - 1) % ngridmax) + 1;
        int ind_f = father[curr_ig - 1];
        while (ind_f > 0) {
            ilevel++;
            if (ind_f <= ncoarse) break;
            curr_ig = ((ind_f - ncoarse - 1) % ngridmax) + 1;
            ind_f = father[curr_ig - 1];
        }
    }

    // 2. Find ancestors up to coarse grid
    std::vector<int> ancestors; int curr_c = icell;
    for(int l=ilevel; l>0; --l) {
        ancestors.push_back(curr_c);
        int ig = ((curr_c - ncoarse - 1) % ngridmax) + 1;
        curr_c = father[ig-1];
    }
    std::reverse(ancestors.begin(), ancestors.end());

    // 3. Find 27 neighbors at coarse level
    int idx = curr_c - 1;
    int iz = idx / (nx * ny); int iy = (idx % (nx * ny)) / nx; int ix = idx % nx;
    int n_curr[27];
    for(int k=-1; k<=1; ++k) for(int j=-1; j<=1; ++j) for(int i=-1; i<=1; ++i) {
        int ni = ix + i, nj = iy + j, nk = iz + k;
        // Periodic wrap
        ni = (ni + nx) % nx; nj = (nj + ny) % ny; nk = (nk + nz) % nz;
        n_curr[(k+1)*9 + (j+1)*3 + (i+1)] = 1 + ni + nj * nx + nk * nx * ny;
    }

    // 4. Descend back to ilevel
    for(int l=1; l<=ilevel; ++l) {
        int target_c = ancestors[l-1];
        int pos = ((target_c - ncoarse - 1) / ngridmax) + 1;
        int n_next[27];
        for(int j=0; j<27; ++j) {
            int f_idx = constants::lll[pos-1][j];
            int c_pos = constants::mmm[pos-1][j];
            int ifather_n = n_curr[f_idx-1];
            if (ifather_n > 0 && son[ifather_n-1] > 0) {
                n_next[j] = ncoarse + (c_pos - 1) * ngridmax + son[ifather_n-1];
            } else {
                n_next[j] = ifather_n;
            }
        }
        for(int i=0; i<27; ++i) n_curr[i] = n_next[i];
    }
    for(int i=0; i<27; ++i) nbors[i] = n_curr[i];
}

} // namespace ramses
