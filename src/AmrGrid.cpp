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
    nbor.assign(6 * ngridmax, 0);
    flag1.assign(ncell, 0);
    cpu_map.assign(ncell, 0);

    xg.assign(3 * ngridmax, 0.0);
    rho.assign(ncell, 0.0);
    phi.assign(ncell, 0.0);
    f_vec.assign(3 * ncell, 0.0);
    uold_vec.assign(nvar * ncell, 0.0);
    unew_vec.assign(nvar * ncell, 0.0);

    // Particles (default allocation, can be resized)
    npart = 0; npartmax = ngridmax * 10; // Rule of thumb
    xp.assign(3 * npartmax, 0.0);
    vp.assign(3 * npartmax, 0.0);
    mp.assign(npartmax, 0.0);
    idp.assign(npartmax, 0);
    levelp.assign(npartmax, 0);
    headp.assign(ngridmax, 0);
    tailp.assign(ngridmax, 0);
    numbp.assign(ngridmax, 0);
    nextp.assign(npartmax, 0);
    prevp.assign(npartmax, 0);

    headp_free = 1; tailp_free = npartmax;
    numbp_free = npartmax;
    for (int i = 1; i <= npartmax; ++i) {
        nextp[i - 1] = (i < npartmax) ? i + 1 : 0;
        prevp[i - 1] = (i > 1) ? i - 1 : 0;
    }

    headf = 1; tailf = ngridmax;
    for (int i = 1; i <= ngridmax; ++i) {
        next[i - 1] = (i < ngridmax) ? i + 1 : 0;
        prev[i - 1] = (i > 1) ? i - 1 : 0;
    }
    numbf = ngridmax;
}

void AmrGrid::resize_particles(int new_npartmax) {
    if (new_npartmax <= npartmax) return;
    int old_npartmax = npartmax;
    std::vector<real_t> new_xp(3 * new_npartmax, 0.0);
    std::vector<real_t> new_vp(3 * new_npartmax, 0.0);
    for (int d = 0; d < 3; ++d) {
        for (int i = 0; i < old_npartmax; ++i) {
            new_xp[d * new_npartmax + i] = xp[d * old_npartmax + i];
            new_vp[d * new_npartmax + i] = vp[d * old_npartmax + i];
        }
    }
    xp = std::move(new_xp);
    vp = std::move(new_vp);
    mp.resize(new_npartmax, 0.0);
    idp.resize(new_npartmax, 0);
    levelp.resize(new_npartmax, 0);
    
    std::vector<int> new_nextp(new_npartmax, 0);
    std::vector<int> new_prevp(new_npartmax, 0);
    for(int i=0; i<old_npartmax; ++i) { new_nextp[i] = nextp[i]; new_prevp[i] = prevp[i]; }
    
    // Add new space to free list
    for (int i = old_npartmax + 1; i <= new_npartmax; ++i) {
        new_nextp[i - 1] = (i < new_npartmax) ? i + 1 : 0;
        new_prevp[i - 1] = (i > old_npartmax + 1) ? i - 1 : 0;
    }
    
    if (numbp_free == 0) {
        headp_free = old_npartmax + 1;
        tailp_free = new_npartmax;
    } else {
        new_nextp[tailp_free - 1] = old_npartmax + 1;
        new_prevp[old_npartmax] = tailp_free;
        tailp_free = new_npartmax;
    }
    numbp_free += (new_npartmax - old_npartmax);
    
    nextp = std::move(new_nextp);
    prevp = std::move(new_prevp);
    npartmax = new_npartmax;
}

int AmrGrid::get_free_grid() {
    if (headf == 0) return 0;
    int igrid = headf;
    headf = next[igrid - 1];
    if (headf == 0) tailf = 0; else prev[headf - 1] = 0;
    next[igrid - 1] = 0; prev[igrid - 1] = 0;
    numbf--;
    return igrid;
}

void AmrGrid::free_grid(int igrid) {
    if (igrid <= 0) return;
    for (int i = 0; i < 6; ++i) nbor[i * ngridmax + igrid - 1] = 0;
    father[igrid - 1] = 0;
    for (int ic = 1; ic <= constants::twotondim; ++ic) son[ncoarse + (ic - 1) * ngridmax + igrid - 1] = 0;
    if (tailf == 0) { headf = igrid; tailf = igrid; prev[igrid - 1] = 0; next[igrid - 1] = 0; }
    else { next[tailf - 1] = igrid; prev[igrid - 1] = tailf; next[igrid - 1] = 0; tailf = igrid; }
    numbf++;
}

void AmrGrid::add_to_level_list(int igrid, int ilevel) {
    int myid = 1;
    int head = get_headl(myid, ilevel);
    if (head == 0) {
        headl_vec[(ilevel - 1) * ncpu + (myid - 1)] = igrid;
        taill_vec[(ilevel - 1) * ncpu + (myid - 1)] = igrid;
        prev[igrid - 1] = 0; next[igrid - 1] = 0;
    } else {
        int tail = taill(myid, ilevel);
        next[tail - 1] = igrid; prev[igrid - 1] = tail; next[igrid - 1] = 0;
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
    for (int i = 1; i <= 6; ++i) ign[i] = nbor[(i - 1) * ngridmax + igrid - 1];
}

void AmrGrid::get_nbor_cells(const int ign[7], int ic, int icn[6], int igrid) const {
    for (int i = 0; i < 6; ++i) icn[i] = 0;
    if (igrid <= 0) return;
    for (int idim = 0; idim < NDIM; ++idim) {
        for (int inbor = 0; inbor < 2; ++inbor) {
            int ig_idx = constants::iii[idim][inbor][ic - 1];
            int ic_pos = constants::jjj[idim][inbor][ic - 1];
            int ig = ign[ig_idx];
            if (ig > 0) icn[idim * 2 + inbor] = ncoarse + (ic_pos - 1) * ngridmax + ig - 1 + 1;
            else icn[idim * 2 + inbor] = - (idim * 2 + inbor + 1);
        }
    }
}

void AmrGrid::get_27_cell_neighbors(int icell, int nbors[27]) const {
    for(int i=0; i<27; ++i) nbors[i] = 0;
    if (icell <= 0) return;
    
    nbors[13] = icell;
    int ig = (icell > ncoarse) ? ((icell - ncoarse - 1) % ngridmax) + 1 : 0;
    int ic = (icell > ncoarse) ? ((icell - ncoarse - 1) / ngridmax) + 1 : icell;
    
    int ign[7];
    if (ig > 0) get_nbor_grids(ig, ign);
    else {
        ign[0] = 0; // Coarse level
        for(int i=1; i<=6; ++i) ign[i] = 0; // Needs coarse neighbor logic
    }
    
    // Simplification: only fill 6 direct neighbors for now
    int icn[6];
    if (ig > 0) {
        get_nbor_cells(ign, ic, icn, ig);
        nbors[13+1] = icn[1]; nbors[13-1] = icn[0];
        if (NDIM > 1) { nbors[13+3] = icn[3]; nbors[13-3] = icn[2]; }
        if (NDIM > 2) { nbors[13+9] = icn[5]; nbors[13-9] = icn[4]; }
    }
}

} // namespace ramses
