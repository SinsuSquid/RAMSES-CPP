#include "ramses/AmrGrid.hpp"
#include "ramses/Constants.hpp"
#include "ramses/MpiManager.hpp"
#ifdef RAMSES_USE_MPI
#include <mpi.h>
#endif
#include <algorithm>
#include <cmath>
#include <iostream>

namespace ramses {

void AmrGrid::allocate(int nx_val, int ny_val, int nz_val, int ngridmax_val, int nvar_val, int ncpu_val, int nlevelmax_val) {
    nx = std::max(nx_val, 1); ny = std::max(ny_val, 1); nz = std::max(nz_val, 1);
    ngridmax = std::max(ngridmax_val, 1); nvar = std::max(nvar_val, 1); ncpu = std::max(ncpu_val, 1); nlevelmax = std::max(nlevelmax_val, 1);
    ncoarse = (long long)nx * ny * nz;
    ncell = ncoarse + (long long)constants::twotondim * ngridmax;
    std::cout << "[AmrGrid] Allocating ngridmax=" << ngridmax << " ncell=" << ncell << " nvar=" << nvar << std::endl;

    headl_vec.assign(ncpu * (nlevelmax + 1), 0);
    taill_vec.assign(ncpu * (nlevelmax + 1), 0);
    numbl_vec.assign(ncpu * (nlevelmax + 1), 0);

    next.assign(ngridmax, 0);
    prev.assign(ngridmax, 0);
    father.assign(ngridmax, 0);
    son.assign(ncell, 0);
    nbor.assign(6 * ngridmax, 0);
    flag1.assign(ncell, 0);
    flag2.assign(ncell, 0);
    cpu_map.assign(ncell, MpiManager::instance().rank() + 1);

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
    tp.assign(npartmax, 0.0);
    zp.assign(npartmax, 0.0);
    family.assign(npartmax, 0);
    tag.assign(npartmax, 0);
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
    tp.resize(new_npartmax, 0.0);
    zp.resize(new_npartmax, 0.0);
    family.resize(new_npartmax, 0);
    tag.resize(new_npartmax, 0);
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

int AmrGrid::get_free_particle() {
    if (headp_free == 0) return 0;
    int ip = headp_free;
    headp_free = nextp[ip - 1];
    if (headp_free == 0) tailp_free = 0; else prevp[headp_free - 1] = 0;
    nextp[ip - 1] = 0; prevp[ip - 1] = 0;
    numbp_free--;
    npart++;
    return ip;
}

void AmrGrid::free_particle(int ip) {
    if (ip <= 0) return;
    if (tailp_free == 0) { headp_free = ip; tailp_free = ip; prevp[ip - 1] = 0; nextp[ip - 1] = 0; }
    else { nextp[tailp_free - 1] = ip; prevp[ip - 1] = tailp_free; nextp[ip - 1] = 0; tailp_free = ip; }
    numbp_free++;
    npart--;
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
    int myid = MpiManager::instance().rank() + 1;
    int head = get_headl(myid, ilevel);
    if (head == 0) {
        headl_vec[ilevel * ncpu + (myid - 1)] = igrid;
        taill_vec[ilevel * ncpu + (myid - 1)] = igrid;
        prev[igrid - 1] = 0; next[igrid - 1] = 0;
    } else {
        int tail = taill(myid, ilevel);
        next[tail - 1] = igrid; prev[igrid - 1] = tail; next[igrid - 1] = 0;
        taill_vec[ilevel * ncpu + (myid - 1)] = igrid;
    }
    numbl_vec[ilevel * ncpu + (myid - 1)]++;
}

int AmrGrid::count_grids_at_level(int ilevel) const {
    int total = 0;
    for (int i = 1; i <= ncpu; ++i) total += numbl(i, ilevel);
    return total;
}

void AmrGrid::synchronize_level_counts() {
#ifdef RAMSES_USE_MPI
    if (ncpu > 1) {
        std::vector<int> global_numbl(ncpu * (nlevelmax + 1));
        MPI_Allreduce(numbl_vec.data(), global_numbl.data(), ncpu * (nlevelmax + 1), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        numbl_vec = std::move(global_numbl);
    }
#endif
}

void AmrGrid::get_cell_center(int icell, real_t xc[3]) const {
    // Default: center of box for unused dimensions
    xc[0] = 0.5 * boxlen; xc[1] = 0.5 * boxlen; xc[2] = 0.5 * boxlen;
    if (icell <= ncoarse) {
        int idx = icell - 1;
        int iz = idx / (nx * ny); idx %= (nx * ny);
        int iy = idx / nx; int ix = idx % nx;
        real_t dx_coarse = boxlen / std::max({nx, ny, nz});
        xc[0] = (ix + 0.5) * dx_coarse;
        if (NDIM > 1) xc[1] = (iy + 0.5) * dx_coarse;
        if (NDIM > 2) xc[2] = (iz + 0.5) * dx_coarse;
    } else {
        int ig = ((icell - ncoarse - 1) % ngridmax) + 1;
        int ic = ((icell - ncoarse - 1) / ngridmax) + 1;
        int idc = (ic - 1);
        int ix = idc & 1, iy = (idc & 2) >> 1, iz = (idc & 4) >> 2;
        int ixyz[3] = {ix, iy, iz};
        // We need level to compute dx. But we don't store level in cell!
        // This is a problem. RAMSES usually knows the level from context.
        // For now, let's assume level is 2 (first refined) if we don't know.
        // Better: find level by traversing upward to coarse cell.
        int level = 1; int curr_ig = ig;
        while(curr_ig > 0) {
            int father_cell = father[curr_ig-1];
            if (father_cell <= ncoarse) break;
            curr_ig = ((father_cell - ncoarse - 1) % ngridmax) + 1;
            level++;
        }
        real_t dx = boxlen / (real_t)(nx * (1 << level));
        real_t scale = boxlen / (real_t)nx;
        // Only compute coordinates for active dimensions; unused dims default to 0.5*boxlen
        xc[0] = 0.5 * boxlen; xc[1] = 0.5 * boxlen; xc[2] = 0.5 * boxlen;
        for(int d=0; d<NDIM; ++d) xc[d] = xg[d * ngridmax + ig - 1] * scale + (real_t)(ixyz[d] - 0.5) * dx;
    }
}

int AmrGrid::find_cell_by_coords(const real_t x[3], int level) const {
    real_t dx_coarse = boxlen / std::max({nx, ny, nz});
    int ix = std::clamp((int)(x[0] / dx_coarse), 0, nx - 1);
    int iy = std::clamp((int)(x[1] / dx_coarse), 0, ny - 1);
    int iz = std::clamp((int)(x[2] / dx_coarse), 0, nz - 1);
    int curr_cell = iz * nx * ny + iy * nx + ix + 1;

    int curr_level = 1;
    while (curr_level < level && son[curr_cell - 1] > 0) {
        int ig = son[curr_cell - 1];
        curr_level++;
        real_t dx = boxlen / (real_t)(nx * (1 << (curr_level - 1)));
        int ixc = (x[0] > xg[0 * ngridmax + ig - 1]) ? 1 : 0;
        int iyc = (x[1] > xg[1 * ngridmax + ig - 1]) ? 1 : 0;
        int izc = (x[2] > xg[2 * ngridmax + ig - 1]) ? 1 : 0;
        int ic = izc * 4 + iyc * 2 + ixc + 1;
        curr_cell = ncoarse + (ic - 1) * ngridmax + ig;
    }
    return curr_cell;
}

void AmrGrid::get_nbor_grids(int igrid, int ign[7]) const {
    ign[0] = igrid;
    for (int i = 1; i <= 6; ++i) {
        int idx = (i - 1) * ngridmax + igrid - 1;
        int val = nbor.at(idx);
        if (val > ncell) {
            ign[i] = 0; // Nullify corrupted pointer
        } else {
            ign[i] = val;
        }
    }
}

void AmrGrid::get_nbor_cells(const int ign[7], int ic, int icn[6], int igrid) const {
    for (int i = 0; i < 6; ++i) icn[i] = 0;
    if (igrid <= 0) return;
    int ifather = father[igrid - 1];

    for (int idim = 0; idim < NDIM; ++idim) {
        for (int inbor = 0; inbor < 2; ++inbor) {
            int ig_idx = constants::iii[idim][inbor][ic - 1];
            int ic_pos = constants::jjj[idim][inbor][ic - 1];
            int ig = ign[ig_idx];
            if (ig > 0) {
                // Neighbor grid exists at same level
                icn[idim * 2 + inbor] = ncoarse + (ic_pos - 1) * ngridmax + ig;
            } else if (ig < 0) {
                // Boundary
                icn[idim * 2 + inbor] = ig;
            } else {
                // Coarse neighbor
                int ifn[6];
                if (ifather <= ncoarse) {
                    get_nbor_cells_coarse(ifather, ifn);
                } else {
                    int p_ig = ((ifather - ncoarse - 1) % ngridmax) + 1;
                    int p_ign[7]; get_nbor_grids(p_ig, p_ign);
                    int p_ic = ((ifather - ncoarse - 1) / ngridmax) + 1;
                    get_nbor_cells(p_ign, p_ic, ifn, p_ig);
                }
                icn[idim * 2 + inbor] = ifn[idim * 2 + inbor];
            }
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

void AmrGrid::get_nbor_cells_coarse(int icell, int icn[6]) const {
    for (int i = 0; i < 6; ++i) icn[i] = 0;
    if (icell <= 0 || icell > ncoarse) return;
    
    int idx = icell - 1;
    int iz = idx / (nx * ny); int rem = idx % (nx * ny);
    int iy = rem / nx; int ix = rem % nx;
    int ixyz[3] = {ix, iy, iz};
    int n[3] = {nx, ny, nz};

    for (int idim = 0; idim < 3; ++idim) {
        for (int side = 0; side < 2; ++side) {
            int coord[3] = {ix, iy, iz};
            coord[idim] += (side == 0) ? -1 : 1;
            
            if (coord[idim] < 0 || coord[idim] >= n[idim]) {
                // Check if boundary exists
                int b_idx = -1;
                for (int ib = 0; ib < nboundary; ++ib) {
                    if (idim == 0 && ((side == 0 && ibound_min[ib] == -1) || (side == 1 && ibound_max[ib] == 1))) b_idx = ib;
                    if (idim == 1 && ((side == 0 && jbound_min[ib] == -1) || (side == 1 && jbound_max[ib] == 1))) b_idx = ib;
                    if (idim == 2 && ((side == 0 && kbound_min[ib] == -1) || (side == 1 && kbound_max[ib] == 1))) b_idx = ib;
                }
                if (b_idx != -1) icn[idim * 2 + side] = -(b_idx + 1);
                else icn[idim * 2 + side] = ((coord[idim] + n[idim]) % n[idim]) + 1; // Fallback to periodic? Or 0?
                // Actually, if no boundary is defined, it should probably be periodic or 0. 
                // Legacy RAMSES defaults coarse level to periodic if no boundary patches are there.
                if (b_idx == -1) {
                    int p_ixyz[3] = {ix, iy, iz};
                    p_ixyz[idim] = (p_ixyz[idim] + (side == 0 ? -1 : 1) + n[idim]) % n[idim];
                    icn[idim * 2 + side] = p_ixyz[2] * nx * ny + p_ixyz[1] * nx + p_ixyz[0] + 1;
                }
            } else {
                icn[idim * 2 + side] = coord[2] * nx * ny + coord[1] * nx + coord[0] + 1;
            }
        }
    }
}

} // namespace ramses
