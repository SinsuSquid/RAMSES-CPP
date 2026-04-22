#include "ramses/AmrGrid.hpp"

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
    // icell_pos: 1-8
    // neighbors: 1-6
    for (int inbor = 0; inbor < 2; ++inbor) { // left, right
        for (int idim = 0; idim < NDIM; ++idim) { // x, y, z
            int nbor_idx = idim * 2 + inbor; // 0-5
            
            // Map to Constants.cpp layout: [ndim][left/right][pos]
            int ig = constants::iii[idim][inbor][icell_pos - 1];
            int ih = constants::jjj[idim][inbor][icell_pos - 1];
            
            if (igridn[ig] > 0) {
                // Absolute cell index = ncoarse + (cell_pos_in_oct - 1) * ngridmax + igrid
                icelln[nbor_idx] = ncoarse + (ih - 1) * ngridmax + igridn[ig];
            } else {
                icelln[nbor_idx] = 0;
            }
        }
    }
}


} // namespace ramses
