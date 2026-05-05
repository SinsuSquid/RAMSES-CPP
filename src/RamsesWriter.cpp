#include "ramses/RamsesWriter.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Parameters.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <sstream>
#include <cstring>
#include <cmath>

namespace ramses {

RamsesWriter::RamsesWriter(const std::string& filename) : filename_(filename) {}

void RamsesWriter::write_amr(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary);
    auto write_rec = [&](const void* data, size_t size) { int s = size; file.write((char*)&s, 4); file.write((char*)data, size); file.write((char*)&s, 4); };
    
    write_rec(&grid.ncpu, 4);
    int ndim = NDIM; write_rec(&ndim, 4);
    int nxnyz[3] = {params::nx, params::ny, params::nz}; write_rec(nxnyz, 12);
    write_rec(&params::nlevelmax, 4);
    write_rec(&grid.ngridmax, 4);
    write_rec(&grid.nboundary, 4);
    int ngrid_tot = 0; for(int l=1; l<=grid.nlevelmax; ++l) ngrid_tot += grid.numbl(1, l);
    write_rec(&ngrid_tot, 4);
    write_rec(&params::boxlen, 8);
    
    int nout[3] = {info.noutput, info.iout, info.iout}; write_rec(nout, 12);
    std::vector<real_t> tout(info.noutput, info.t); write_rec(tout.data(), tout.size() * 8);
    std::vector<real_t> aout(info.noutput, 1.0); write_rec(aout.data(), aout.size() * 8);
    write_rec(&info.t, 8);
    std::vector<real_t> dt_vec(params::nlevelmax, 1e-3); write_rec(dt_vec.data(), dt_vec.size() * 8);
    write_rec(dt_vec.data(), dt_vec.size() * 8);
    int step[2] = {info.nstep, 0}; write_rec(step, 8);
    real_t cosmo[3] = {0, 1.0, 1.0}; write_rec(cosmo, 24);
    real_t cosmo2[7] = {0.3, 0.7, 0.0, 0.045, 0.7, 1.0, params::boxlen}; write_rec(cosmo2, 56);
    real_t cosmo3[5] = {1.0, 0.0, 1.0, 0.0, 0.0}; write_rec(cosmo3, 40);
    real_t mass_sph = 0; write_rec(&mass_sph, 8);
    
    std::vector<int> headl(grid.ncpu * grid.nlevelmax);
    for (int l = 1; l <= grid.nlevelmax; ++l) headl[(l-1)*grid.ncpu] = grid.get_headl(1, l);
    write_rec(headl.data(), headl.size() * 4);
    std::vector<int> taill(grid.ncpu * grid.nlevelmax, 0); write_rec(taill.data(), taill.size() * 4);
    std::vector<int> numbl(grid.ncpu * grid.nlevelmax);
    for (int l = 1; l <= grid.nlevelmax; ++l) numbl[(l-1)*grid.ncpu] = grid.numbl(1, l);
    write_rec(numbl.data(), numbl.size() * 4);
    
    std::vector<int> numbf(grid.ncpu, 0); write_rec(numbf.data(), numbf.size() * 4);
    std::vector<int> headf(grid.ncpu, 0); write_rec(headf.data(), headf.size() * 4);
    std::vector<int> tailf(grid.ncpu, 0); write_rec(tailf.data(), tailf.size() * 4);
    int ngrid_curr = ngrid_tot; write_rec(&ngrid_curr, 4);
    std::vector<int> cp_map(ngrid_tot, 1); write_rec(cp_map.data(), cp_map.size() * 4);

    for (int il = 1; il <= grid.nlevelmax; ++il) {
        int myid = 1; int ngrid = grid.numbl(myid, il);
        if (ngrid > 0) {
            std::vector<int> g_list; int ig = grid.get_headl(myid, il); while(ig > 0) { g_list.push_back(ig); ig = grid.next[ig-1]; }
            int n = g_list.size(); int s_r = n * 8;
            std::vector<real_t> tmp(n);
            for(int d=0; d<ndim; ++d) { for(int i=0; i<n; ++i) tmp[i] = grid.xg[d*grid.ngridmax + g_list[i]-1]; write_rec(tmp.data(), s_r); }
            std::vector<int> tmp_i(n);
            for(int i=0; i<n; ++i) tmp_i[i] = grid.father[g_list[i]-1]; write_rec(tmp_i.data(), n*4);
            for(int d=0; d<ndim*2; ++d) { for(int i=0; i<n; ++i) tmp_i[i] = grid.nbor[d*grid.ngridmax + g_list[i]-1]; write_rec(tmp_i.data(), n*4); }
            for(int d=0; d<constants::twotondim; ++d) { for(int i=0; i<n; ++i) tmp_i[i] = grid.son[grid.ncoarse + d*grid.ngridmax + g_list[i]-1]; write_rec(tmp_i.data(), n*4); }
            std::vector<int> cp_m(n, 1); write_rec(cp_m.data(), n*4);
        }
    }
}

void RamsesWriter::write_hydro(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary);
    auto write_rec = [&](const void* data, size_t size) { int s = size; file.write((char*)&s, 4); file.write((char*)data, size); file.write((char*)&s, 4); };
    
    write_rec(&grid.ncpu, 4); write_rec(&grid.nvar, 4);
    int ngrid_tot = 0; for(int l=1; l<=grid.nlevelmax; ++l) ngrid_tot += grid.numbl(1, l);
    write_rec(&ngrid_tot, 4);
    write_rec(&grid.gamma, 8);

    for (int il = 1; il <= grid.nlevelmax; ++il) {
        int myid = 1; int ngrid = grid.numbl(myid, il);
        if (ngrid > 0) {
            std::vector<int> g_list; int ig = grid.get_headl(myid, il); while(ig > 0) { g_list.push_back(ig); ig = grid.next[ig-1]; }
            int n = g_list.size();
            for(int iv=1; iv<=grid.nvar; ++iv) {
                for(int ic=0; ic<constants::twotondim; ++ic) {
                    std::vector<real_t> tmp(n);
                    for(int i=0; i<n; ++i) tmp[i] = grid.uold(grid.ncoarse + ic*grid.ngridmax + g_list[i]-1 + 1, iv);
                    write_rec(tmp.data(), n * 8);
                }
            }
        }
    }
}

void RamsesWriter::write_header(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_);
    file << "ncpu        = " << grid.ncpu << std::endl;
    file << "ndim        = " << NDIM << std::endl;
    file << "levelmin    = " << params::levelmin << std::endl;
    file << "levelmax    = " << params::nlevelmax << std::endl;
    file << "ngridmax    = " << grid.ngridmax << std::endl;
    file << "nstep       = " << info.nstep << std::endl;
    file << "nstep_coarse= " << info.nstep << std::endl;
    file << "noutput     = " << info.noutput << std::endl;
    file << "time        = " << std::scientific << std::setprecision(16) << info.t << std::endl;
    file << "boxlen      = " << params::boxlen << std::endl;
    file << "omega_m     = 0.3" << std::endl;
    file << "omega_l     = 0.7" << std::endl;
    file << "omega_k     = 0.0" << std::endl;
    file << "omega_b     = 0.045" << std::endl;
    file << "h0          = 0.7" << std::endl;
    file << "aexp        = 1.0" << std::endl;
    file << "unit_l      = 1.0" << std::endl;
    file << "unit_d      = 1.0" << std::endl;
    file << "unit_t      = 1.0" << std::endl;
    file << "ordering type= hilbert" << std::endl;
}

void RamsesWriter::write_header_file(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_);
    file << "# This is a fake header file for visu_ramses" << std::endl;
}

void RamsesWriter::write_hydro_descriptor(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_);
    file << "# version: 1" << std::endl;
    int ivar = 1;
    file << ivar++ << ", density, double" << std::endl;
    file << ivar++ << ", velocity_x, double" << std::endl;
    if (NDIM > 1) file << ivar++ << ", velocity_y, double" << std::endl; else ivar++;
    if (NDIM > 2) file << ivar++ << ", velocity_z, double" << std::endl; else ivar++;
    file << ivar++ << ", pressure, double" << std::endl;
    for (int ie = 1; ie <= info.nener; ++ie) { std::stringstream ss; ss << "non_thermal_energy_" << std::setfill('0') << std::setw(2) << (ie - 1); file << ivar++ << ", " << ss.str() << ", double" << std::endl; }
    int npassive = grid.nvar - (5 + info.nener);
    for (int ip = 1; ip <= npassive; ++ip) { std::stringstream ss; ss << "scalar_" << std::setfill('0') << std::setw(2) << (ip - 1); file << ivar++ << ", " << ss.str() << ", double" << std::endl; }
}

void RamsesWriter::write_grav(const AmrGrid&, const SnapshotInfo&) {}

} // namespace ramses
