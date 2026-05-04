#include "ramses/RamsesWriter.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Parameters.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <sstream>
#include <cstring>

namespace ramses {

RamsesWriter::RamsesWriter(const std::string& filename) : filename_(filename) {}

template <typename T>
void RamsesWriter::write_record_internal(std::ofstream& file, const T* data, size_t count) {
    int32_t size = static_cast<int32_t>(count * sizeof(T));
    file.write(reinterpret_cast<const char*>(&size), sizeof(int32_t));
    file.write(reinterpret_cast<const char*>(data), size);
    file.write(reinterpret_cast<const char*>(&size), sizeof(int32_t));
}

void RamsesWriter::write_amr(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary | std::ios::trunc); if (!file.is_open()) return;
    int32_t val = grid.ncpu; write_record_internal(file, &val, 1);
    val = NDIM; write_record_internal(file, &val, 1);
    int32_t dims[3] = {grid.nx, grid.ny, grid.nz}; write_record_internal(file, dims, 3);
    val = grid.nlevelmax; write_record_internal(file, &val, 1);
    val = grid.ngridmax; write_record_internal(file, &val, 1);
    val = grid.nboundary; write_record_internal(file, &val, 1);
    val = info.nstep; write_record_internal(file, &val, 1);
    val = info.noutput; write_record_internal(file, &val, 1);
    val = info.nstep_coarse; write_record_internal(file, &val, 1);
    val = info.iout; write_record_internal(file, &val, 1);
    
    double dval = info.t; write_record_internal(file, &dval, 1);
    dval = 0.0; write_record_internal(file, &dval, 1); // dtold_coarse
    write_record_internal(file, &dval, 1); // dtnew_coarse
    std::vector<double> consts(10, 0.0); consts[9] = grid.boxlen; write_record_internal(file, consts.data(), 10);
    
    std::vector<double> dt(grid.nlevelmax, 0.0); write_record_internal(file, dt.data(), dt.size()); // dtold
    write_record_internal(file, dt.data(), dt.size()); // dtnew
    
    std::vector<int32_t> ngridlevel(grid.ncpu * grid.nlevelmax, 0);
    for(int il=1; il<=grid.nlevelmax; ++il) for(int ic=1; ic<=grid.ncpu; ++ic) ngridlevel[(il-1)*grid.ncpu + (ic-1)] = grid.numbl(ic, il);
    write_record_internal(file, ngridlevel.data(), ngridlevel.size());
    
    std::vector<int32_t> ngridboundary(grid.nboundary * grid.nlevelmax, 0);
    if(grid.nboundary > 0) write_record_internal(file, ngridboundary.data(), ngridboundary.size());
    
    char ordering[128] = {0}; std::strncpy(ordering, "hilbert", 127);
    write_record_internal(file, ordering, 128);
    
    for (int il = 1; il <= grid.nlevelmax; ++il) for (int ic = 1; ic <= grid.ncpu; ++ic) {
        int32_t il_32 = il, nca = grid.numbl(ic, il);
        write_record_internal(file, &il_32, 1); write_record_internal(file, &nca, 1);
        if (nca > 0) {
            int cr;
            for (int d = 0; d < NDIM; ++d) { std::vector<double> xg(nca); cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { xg[i] = grid.xg[d * grid.ngridmax + cr - 1]; cr = grid.next[cr - 1]; } write_record_internal(file, xg.data(), nca); }
            std::vector<int32_t> id(nca);
            cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { id[i] = grid.father[cr - 1]; cr = grid.next[cr - 1]; } write_record_internal(file, id.data(), nca);
            cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { id[i] = grid.next[cr - 1]; cr = grid.next[cr - 1]; } write_record_internal(file, id.data(), nca);
            cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { id[i] = grid.prev[cr - 1]; cr = grid.next[cr - 1]; } write_record_internal(file, id.data(), nca);
            for (int d = 0; d < (1 << NDIM); ++d) { cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { id[i] = grid.son[grid.ncoarse + d * grid.ngridmax + cr - 1]; cr = grid.next[cr - 1]; } write_record_internal(file, id.data(), nca); }
            for (int d = 0; d < 6; ++d) { cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { id[i] = grid.nbor[d * grid.ngridmax + cr - 1]; cr = grid.next[cr - 1]; } write_record_internal(file, id.data(), nca); }
            cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { id[i] = grid.cpu_map[cr - 1]; cr = grid.next[cr - 1]; } write_record_internal(file, id.data(), nca);
        }
    }
    file.close();
}

void RamsesWriter::write_hydro(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary | std::ios::trunc); if (!file.is_open()) return;
    int32_t val = grid.ncpu; write_record_internal(file, &val, 1);
    int32_t nv = 5; 
#ifdef MHD
    nv = 11;
#endif
    write_record_internal(file, &nv, 1);
    val = NDIM; write_record_internal(file, &val, 1);
    val = grid.nlevelmax; write_record_internal(file, &val, 1);
    val = grid.nboundary; write_record_internal(file, &val, 1);
    double gamma = info.gamma; write_record_internal(file, &gamma, 1);
    for (int il = 1; il <= grid.nlevelmax; ++il) for (int ic = 1; ic <= grid.ncpu; ++ic) {
        int32_t il_32 = il, nca = grid.numbl(ic, il);
        write_record_internal(file, &il_32, 1); write_record_internal(file, &nca, 1);
        if (nca > 0) for (int isc = 1; isc <= (1 << NDIM); ++isc) {
            std::vector<double> d(nca);
            int cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { d[i] = grid.uold(grid.ncoarse + (isc - 1) * grid.ngridmax + cr, 1); cr = grid.next[cr - 1]; } write_record_internal(file, d.data(), nca);
            for (int k = 1; k <= 3; ++k) { cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { int id = grid.ncoarse + (isc - 1) * grid.ngridmax + cr; d[i] = (k <= NDIM) ? grid.uold(id, 1 + k) / std::max(grid.uold(id, 1), 1e-10) : 0.0; cr = grid.next[cr - 1]; } write_record_internal(file, d.data(), nca); }
            cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) {
                int id = grid.ncoarse + (isc - 1) * grid.ngridmax + cr;
                real_t r = std::max(grid.uold(id, 1), 1e-10), v2 = 0; for (int j = 1; j <= NDIM; ++j) { real_t v = grid.uold(id, 1 + j) / r; v2 += v * v; }
                real_t em = 0; 
#ifdef MHD
                real_t Ax = 0.5 * (grid.uold(id, 6) + grid.uold(id, grid.nvar - 2)), Ay = 0.5 * (grid.uold(id, 7) + grid.uold(id, grid.nvar - 1)), Az = 0.5 * (grid.uold(id, 8) + grid.uold(id, grid.nvar)); em = 0.5 * (Ax * Ax + Ay * Ay + Az * Az);
#endif
                d[i] = (info.gamma - 1.0) * (grid.uold(id, 5) - 0.5 * r * v2 - em); cr = grid.next[cr - 1];
            } write_record_internal(file, d.data(), nca);
#ifdef MHD
            for (int k = 1; k <= 3; ++k) { cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { d[i] = grid.uold(grid.ncoarse + (isc - 1) * grid.ngridmax + cr, 5 + k); cr = grid.next[cr - 1]; } write_record_internal(file, d.data(), nca); }
            for (int k = 1; k <= 3; ++k) { cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { d[i] = grid.uold(grid.ncoarse + (isc - 1) * grid.ngridmax + cr, grid.nvar - 3 + k); cr = grid.next[cr - 1]; } write_record_internal(file, d.data(), nca); }
#endif
        }
    }
}

void RamsesWriter::write_grav(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary | std::ios::trunc); if (!file.is_open()) return;
    int32_t val = grid.ncpu; write_record_internal(file, &val, 1);
    val = NDIM + 1; write_record_internal(file, &val, 1);
    val = grid.nlevelmax; write_record_internal(file, &val, 1);
    val = 0; write_record_internal(file, &val, 1);
    for (int il = 1; il <= grid.nlevelmax; ++il) for (int ic = 1; ic <= grid.ncpu; ++ic) {
        int32_t il_32 = il, nca = grid.numbl(ic, il);
        write_record_internal(file, &il_32, 1); write_record_internal(file, &nca, 1);
        if (nca > 0) for (int is = 1; is <= (1 << NDIM); ++is) {
            std::vector<double> pd(nca); int cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { pd[i] = grid.phi[grid.ncoarse + (is - 1) * grid.ngridmax + cr - 1]; cr = grid.next[cr - 1]; } write_record_internal(file, pd.data(), nca);
            for (int d = 1; d <= NDIM; ++d) { std::vector<double> fd(nca); cr = grid.headl(ic, il); for (int i = 0; i < nca; ++i) { fd[i] = grid.f(grid.ncoarse + (is - 1) * grid.ngridmax + cr - 1, d); cr = grid.next[cr - 1]; } write_record_internal(file, fd.data(), nca); }
        }
    }
}

void RamsesWriter::write_header(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_); if (!file.is_open()) return;
    file << "ncpu        = " << std::setw(15) << grid.ncpu << std::endl << "ndim        = " << std::setw(15) << NDIM << std::endl << "levelmin    = " << std::setw(15) << 1 << std::endl << "levelmax    = " << std::setw(15) << grid.nlevelmax << std::endl << "ngridmax    = " << std::setw(15) << grid.ngridmax << std::endl << "nstep       = " << std::setw(15) << info.nstep << std::endl << "nstep_coarse=" << std::setw(15) << info.nstep_coarse << std::endl << "noutput     =" << std::setw(15) << info.noutput << std::endl << "time        = " << std::setw(15) << std::scientific << info.t << std::endl << "boxlen      = " << std::setw(15) << std::scientific << grid.boxlen << std::endl << "omega_m     = 0.0" << std::endl << "omega_l     = 0.0" << std::endl << "omega_k     = 0.0" << std::endl << "omega_b     = 0.0" << std::endl << "h0          = 0.0" << std::endl << "aexp        = 1.0" << std::endl << "unit_l      = 1.0" << std::endl << "unit_d      = 1.0" << std::endl << "unit_t      = 1.0" << std::endl << std::endl << "ordering type=hilbert" << std::endl << std::endl;
}

void RamsesWriter::write_header_file(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary | std::ios::trunc); if (!file.is_open()) return;
    int32_t val = grid.ncpu; write_record_internal(file, &val, 1);
    val = grid.nlevelmax; write_record_internal(file, &val, 1);
    std::vector<int32_t> counts(grid.ncpu * grid.nlevelmax);
    for (int il = 1; il <= grid.nlevelmax; ++il) for (int ic = 1; ic <= grid.ncpu; ++ic) counts[(il-1)*grid.ncpu + (ic-1)] = grid.numbl(ic, il);
    write_record_internal(file, counts.data(), counts.size());
}

void RamsesWriter::write_hydro_descriptor(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_); if (!file.is_open()) return;
    bool mhd = false;
#ifdef MHD
    mhd = true;
#endif
    int ivar = 1; file << ivar++ << ", density, double" << std::endl;
    file << ivar++ << ", velocity_x, double" << std::endl << ivar++ << ", velocity_y, double" << std::endl << ivar++ << ", velocity_z, double" << std::endl << ivar++ << ", pressure, double" << std::endl;
    if (mhd) { 
        file << ivar++ << ", B_x_left, double" << std::endl << ivar++ << ", B_y_left, double" << std::endl << ivar++ << ", B_z_left, double" << std::endl;
        file << ivar++ << ", B_x_right, double" << std::endl << ivar++ << ", B_y_right, double" << std::endl << ivar++ << ", B_z_right, double" << std::endl;
    }
}

} // namespace ramses
