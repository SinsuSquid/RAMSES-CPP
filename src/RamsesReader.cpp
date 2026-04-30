#include "ramses/RamsesReader.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Parameters.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

namespace ramses {

bool RamsesReader::load_amr(AmrGrid& grid) {
    if (!file_.is_open()) return false;

    int ncpu = read_single<int>();
    int ndim = read_single<int>();
    std::vector<int> nxyz; read_record(nxyz);
    int nx = nxyz[0], ny = nxyz[1], nz = nxyz[2];
    int nlevelmax = read_single<int>();
    int ngridmax = read_single<int>();
    int nboundary = read_single<int>();
    int ngrid_current = read_single<int>();
    real_t boxlen = read_single<real_t>();

    grid.allocate(nx, ny, nz, ngridmax, grid.nvar > 0 ? grid.nvar : 5, ncpu, nlevelmax);
    grid.ndim = ndim;

    for(int i=0; i<3; ++i) skip_record(); // nout, tout, aout
    real_t t = read_single<real_t>(); // t
    for(int i=0; i<7; ++i) skip_record(); // dtold, dtnew, nstep, einit, cosmo, aexp, msph

    std::vector<int> buf;
    read_record(buf); for(int i=0; i<ncpu*nlevelmax; ++i) grid.headl.data()[i] = buf[i];
    read_record(buf); for(int i=0; i<ncpu*nlevelmax; ++i) grid.taill.data()[i] = buf[i];
    read_record(buf); for(int i=0; i<ncpu*nlevelmax; ++i) grid.numbl.data()[i] = buf[i];
    skip_record(); // numbtot
    if (nboundary > 0) { skip_record(); skip_record(); skip_record(); }
    skip_record(); // headf...
    skip_record(); // order
    skip_record(); // key_size

    read_record(buf); for(int i=0; i<grid.ncoarse; ++i) grid.son[i+1] = buf[i];
    read_record(buf); for(int i=0; i<grid.ncoarse; ++i) grid.flag1[i+1] = buf[i];
    read_record(buf); for(int i=0; i<grid.ncoarse; ++i) grid.cpu_map[i+1] = buf[i];

    for(int il=1; il <= nlevelmax; ++il) for(int ic=1; ic <= ncpu; ++ic) {
        int nca = grid.numbl(ic, il);
        if (nca > 0) {
            std::vector<int> igs; read_record(igs);
            std::vector<int> nxt; read_record(nxt);
            std::vector<int> prv; read_record(prv);
            for(int i=0; i<nca; ++i) { grid.next[igs[i]-1] = nxt[i]; grid.prev[igs[i]-1] = prv[i]; }
            for(int d=1; d<=ndim; ++d) { std::vector<double> xgd; read_record(xgd); for(int i=0; i<nca; ++i) grid.xg[(d-1)*ngridmax + (igs[i]-1)] = xgd[i]; }
            std::vector<int> ftr; read_record(ftr); for(int i=0; i<nca; ++i) grid.father[igs[i]-1] = ftr[i];
            for(int n=1; n<=2*ndim; ++n) { std::vector<int> nbd; read_record(nbd); for(int i=0; i<nca; ++i) grid.nbor[(n-1)*ngridmax + (igs[i]-1)] = nbd[i]; }
            for(int s=1; s<=constants::twotondim; ++s) { read_record(buf); for(int i=0; i<nca; ++i) grid.son[grid.ncoarse + (s-1)*ngridmax + igs[i]] = buf[i]; }
            for(int s=1; s<=constants::twotondim; ++s) { read_record(buf); for(int i=0; i<nca; ++i) grid.flag1[grid.ncoarse + (s-1)*ngridmax + igs[i]] = buf[i]; }
            for(int s=1; s<=constants::twotondim; ++s) { read_record(buf); for(int i=0; i<nca; ++i) grid.cpu_map[grid.ncoarse + (s-1)*ngridmax + igs[i]] = buf[i]; }
        }
    }
    return true;
}

bool RamsesReader::load_hydro(AmrGrid& grid) {
    if (!file_.is_open()) return false;
    
    int ncpu = read_single<int>();
    int nvar_file = read_single<int>();
    int ndim = read_single<int>();
    int nlevelmax = read_single<int>();
    int nboundary = read_single<int>();
    real_t gamma = read_single<real_t>();

    for (int il = 1; il <= nlevelmax; ++il) {
        for (int ib = 1; ib <= nboundary + ncpu; ++ib) {
            int ilevel = read_single<int>();
            int ncache = read_single<int>();
            if (ncache > 0) {
                std::vector<int> active_grids;
                if (ib <= ncpu) {
                    int curr = grid.headl(ib, il);
                    while (curr > 0) { active_grids.push_back(curr); curr = grid.next[curr-1]; }
                }
                int n_to_read = 1 + ndim; // density + velocities
                // We don't know nener, so we'll try to read as many as possible
                // for (int is = 1; is <= (1<<ndim); ++is) {
                //     for (int iv = 1; iv <= nvar_file; ++iv) { ... }
                // }
                // Standard RAMSES backup_hydro writes: 1 + ndim + nener + 1 + npscal
                
                for (int is = 1; is <= (1<<ndim); ++is) {
                    int iv_in = 1;
                    // 1. Density
                    std::vector<double> data; read_record(data);
                    for(int i=0; i<ncache; ++i) grid.uold(grid.ncoarse + (is-1)*grid.ngridmax + active_grids[i], 1) = data[i];
                    iv_in++;

                    // 2. Velocities
                    for(int d=1; d<=ndim; ++d) {
                        read_record(data);
                        for(int i=0; i<ncache; ++i) grid.uold(grid.ncoarse + (is-1)*grid.ngridmax + active_grids[i], 1+d) = data[i];
                        iv_in++;
                    }

                    // 3. Pressure (we store at 5)
                    read_record(data);
                    for(int i=0; i<ncache; ++i) grid.uold(grid.ncoarse + (is-1)*grid.ngridmax + active_grids[i], 5) = data[i];
                    iv_in++;

                    // Skip remaining records for this cell
                    for(; iv_in <= nvar_file; ++iv_in) skip_record();
                }
            }
        }
    }
    return true;
}

} // namespace ramses
