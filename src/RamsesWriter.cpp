#include "ramses/RamsesWriter.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace ramses {
RamsesWriter::RamsesWriter(const std::string& filename) : filename_(filename) {}
bool RamsesWriter::is_open() const { return true; }
template <typename T> void RamsesWriter::write_record_internal(std::ofstream& file, const T* data, size_t count) {
    int32_t size_bytes = static_cast<int32_t>(count * sizeof(T));
    file.write(reinterpret_cast<const char*>(&size_bytes), sizeof(int32_t));
    if (count > 0) file.write(reinterpret_cast<const char*>(data), size_bytes);
    file.write(reinterpret_cast<const char*>(&size_bytes), sizeof(int32_t));
}
template <typename T> void RamsesWriter::write_record(const T* data, size_t count) {
    std::ofstream file(filename_, std::ios::binary | std::ios::out);
    if (!file.is_open()) return;
    write_record_internal(file, data, count);
    file.close();
}
template void RamsesWriter::write_record<int32_t>(const int32_t*, size_t);
template void RamsesWriter::write_record<double>(const double*, size_t);
template void RamsesWriter::write_record<char>(const char*, size_t);
void RamsesWriter::write_amr(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary | std::ios::trunc);
    if (!file.is_open()) return;
    write_record_internal(file, &grid.ncpu, 1);
    int32_t ndim = NDIM; write_record_internal(file, &ndim, 1);
    int32_t nxnyz[3] = {params::nx, params::ny, params::nz}; write_record_internal(file, nxnyz, 3);
    write_record_internal(file, &grid.nlevelmax, 1);
    write_record_internal(file, &grid.ngridmax, 1);
    write_record_internal(file, &grid.nboundary, 1);
    int32_t ngrid_current = grid.ngridmax - grid.numbf;
    write_record_internal(file, &ngrid_current, 1);
    double boxlen = params::boxlen; write_record_internal(file, &boxlen, 1);
    int32_t rec9[3] = {info.noutput, info.iout, 1}; write_record_internal(file, rec9, 3);
    std::vector<double> tout(info.noutput, 0.0);
    for(int i=0; i<std::min((int)info.tout.size(), info.noutput); ++i) tout[i] = info.tout[i];
    write_record_internal(file, tout.data(), info.noutput);
    std::vector<double> aout(info.noutput, 0.0); write_record_internal(file, aout.data(), info.noutput);
    double t = info.t; write_record_internal(file, &t, 1);
    std::vector<double> dtold(grid.nlevelmax, 0.0); write_record_internal(file, dtold.data(), grid.nlevelmax);
    std::vector<double> dtnew(grid.nlevelmax, 0.0); write_record_internal(file, dtnew.data(), grid.nlevelmax);
    int32_t rec15[2] = {info.nstep, info.nstep_coarse}; write_record_internal(file, rec15, 2);
    double rec16[3]={0,0,0}, rec17[7]={0,0,0,0,0,0,0}, rec18[5]={0,0,0,0,0}, msph=0;
    write_record_internal(file, rec16, 3); write_record_internal(file, rec17, 7); write_record_internal(file, rec18, 5); write_record_internal(file, &msph, 1);
    std::vector<int32_t> hl(grid.ncpu * grid.nlevelmax), tl(grid.ncpu * grid.nlevelmax), nl(grid.ncpu * grid.nlevelmax);
    for(int il=0; il<grid.nlevelmax; ++il) for(int ic=0; ic<grid.ncpu; ++ic) { hl[il*grid.ncpu+ic]=grid.headl(ic+1,il+1); tl[il*grid.ncpu+ic]=grid.taill(ic+1,il+1); nl[il*grid.ncpu+ic]=grid.numbl(ic+1,il+1); }
    write_record_internal(file, hl.data(), hl.size()); write_record_internal(file, tl.data(), tl.size()); write_record_internal(file, nl.data(), nl.size());
    std::vector<int32_t> nt(10 * grid.nlevelmax, 0); for(int il=0; il<grid.nlevelmax; ++il) nt[il * 10] = grid.count_grids_at_level(il+1);
    write_record_internal(file, nt.data(), nt.size());
    if (grid.nboundary > 0) { std::vector<int32_t> b_null(grid.nboundary * grid.nlevelmax, 0); write_record_internal(file, b_null.data(), b_null.size()); write_record_internal(file, b_null.data(), b_null.size()); write_record_internal(file, b_null.data(), b_null.size()); }
    int32_t rec_mem[5] = {grid.headf, grid.tailf, grid.numbf, 0, 0}; write_record_internal(file, rec_mem, 5);
    char ord[128] = {0}; std::strncpy(ord, "hilbert", 127); write_record_internal(file, ord, 128);
    std::vector<double> bound_key(grid.ncpu + 1, 0.0); bound_key[grid.ncpu] = 1.0; write_record_internal(file, bound_key.data(), bound_key.size());
    write_record_internal(file, grid.son.data(), grid.ncoarse); write_record_internal(file, grid.flag1.data(), grid.ncoarse); write_record_internal(file, grid.cpu_map.data(), grid.ncoarse);
    for(int il=1; il <= grid.nlevelmax; ++il) {
        for(int ic=1; ic <= grid.ncpu + grid.nboundary; ++ic) {
            int nca = (ic <= grid.ncpu) ? grid.numbl(ic, il) : 0;
            if (nca > 0) {
                std::vector<int32_t> igs(nca), nxt(nca), prv(nca), ftr(nca), nbd(nca, 0), sd(nca), f1d(nca), cmd(nca); std::vector<double> xgd(nca);
                int curr = grid.headl(ic, il); for(int i=0; i<nca; ++i) { igs[i]=curr; nxt[i]=grid.next[curr-1]; prv[i]=grid.prev[curr-1]; ftr[i]=grid.father[curr-1]; curr=grid.next[curr-1]; }
                write_record_internal(file, igs.data(), nca); write_record_internal(file, nxt.data(), nca); write_record_internal(file, prv.data(), nca);
                for(int d=1; d<=NDIM; ++d) { curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ xgd[i]=grid.xg[(d-1)*grid.ngridmax+(curr-1)]; curr=grid.next[curr-1]; } write_record_internal(file, xgd.data(), nca); }
                write_record_internal(file, ftr.data(), nca); for(int n=1; n<=2*NDIM; ++n) { write_record_internal(file, nbd.data(), nca); }
                for(int s=1; s<=constants::twotondim; ++s) { curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ sd[i]=grid.son[grid.ncoarse+(s-1)*grid.ngridmax+curr-1]; curr=grid.next[curr-1]; } write_record_internal(file, sd.data(), nca); }
                for(int s=1; s<=constants::twotondim; ++s) { curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ cmd[i]=grid.cpu_map[grid.ncoarse+(s-1)*grid.ngridmax+curr-1]; curr=grid.next[curr-1]; } write_record_internal(file, cmd.data(), nca); }
                for(int s=1; s<=constants::twotondim; ++s) { curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ f1d[i]=grid.flag1[grid.ncoarse+(s-1)*grid.ngridmax+curr-1]; curr=grid.next[curr-1]; } write_record_internal(file, f1d.data(), nca); }
            }
        }
    }
    file.close();
}
void RamsesWriter::write_hydro(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary | std::ios::trunc);
    if (!file.is_open()) return;
    bool mhd = false;
#ifdef MHD
    mhd = true;
#endif
    int nener = info.nener; int nvar_f = 1 + NDIM + nener + 1 + ((grid.nvar > (mhd ? 8 : 5) + nener) ? grid.nvar - ((mhd ? 8 : 5) + nener) : 0);
    if (mhd) nvar_f += 6;
    int32_t rec_h[5] = {grid.ncpu, nvar_f, (int32_t)NDIM, grid.nlevelmax, grid.nboundary};
    int32_t s_rec = 5*4 + 8; file.write((char*)&s_rec, 4); file.write((char*)rec_h, 20); file.write((char*)&info.gamma, 8); file.write((char*)&s_rec, 4);
    for(int il=1; il <= grid.nlevelmax; ++il) for(int ic=1; ic <= grid.ncpu + grid.nboundary; ++ic) {
        int nca = (ic <= grid.ncpu) ? grid.numbl(ic, il) : 0;
        write_record_internal(file, &il, 1); write_record_internal(file, &nca, 1);
        if (nca > 0) for(int is=1; is<=constants::twotondim; ++is) {
            std::vector<double> d(nca); int cr = grid.headl(ic, il);
            for(int i=0; i<nca; ++i) { d[i]=grid.uold(grid.ncoarse+(is-1)*grid.ngridmax+cr, 1); cr=grid.next[cr-1]; } write_record_internal(file, d.data(), nca);
            for(int k=1; k<=NDIM; ++k) { cr=grid.headl(ic, il); for(int i=0; i<nca; ++i) { int id=grid.ncoarse+(is-1)*grid.ngridmax+cr; d[i]=grid.uold(id, 1+k)/std::max(grid.uold(id, 1), 1e-10); cr=grid.next[cr-1]; } write_record_internal(file, d.data(), nca); }
            if (mhd) { for(int k=1; k<=3; ++k) { cr=grid.headl(ic, il); for(int i=0; i<nca; ++i) { d[i]=grid.uold(grid.ncoarse+(is-1)*grid.ngridmax+cr, 5+k); cr=grid.next[cr-1]; } write_record_internal(file, d.data(), nca); } for(int k=1; k<=3; ++k) { cr=grid.headl(ic, il); for(int i=0; i<nca; ++i) { d[i]=grid.uold(grid.ncoarse+(is-1)*grid.ngridmax+cr, grid.nvar-3+k); cr=grid.next[cr-1]; } write_record_internal(file, d.data(), nca); } }
            for(int ie=1; ie<=nener; ++ie) { cr=grid.headl(ic, il); for(int i=0; i<nca; ++i) { int id=grid.ncoarse+(is-1)*grid.ngridmax+cr; d[i]=grid.uold(id, (mhd?8:5)+ie)/std::max(grid.uold(id, 1), 1e-10); cr=grid.next[cr-1]; } write_record_internal(file, d.data(), nca); }
            cr = grid.headl(ic, il); for(int i=0; i<nca; ++i) {
                int id=grid.ncoarse+(is-1)*grid.ngridmax+cr; real_t r=std::max(grid.uold(id, 1), 1e-10), v2=0; for(int j=1; j<=3; ++j) { real_t v=grid.uold(id, 1+j)/r; v2+=v*v; }
                real_t em=0, er=0; if(mhd) { real_t Ax=0.5*(grid.uold(id,6)+grid.uold(id,grid.nvar-2)), Ay=0.5*(grid.uold(id,7)+grid.uold(id,grid.nvar-1)), Az=0.5*(grid.uold(id,8)+grid.uold(id,grid.nvar)); em=0.5*(Ax*Ax+Ay*Ay+Az*Az); }
                for(int ie=0; ie<nener; ++ie) er += grid.uold(id, (mhd?8:5)+1+ie);
                d[i]=(info.gamma-1.0)*(grid.uold(id, 5)-0.5*r*v2-em-er); cr=grid.next[cr-1];
            } write_record_internal(file, d.data(), nca);
            int nps = (grid.nvar > (mhd ? 8 : 5) + nener) ? grid.nvar - ((mhd ? 8 : 5) + nener) : 0;
            for(int is=1; is<=nps; ++is) { cr=grid.headl(ic, il); for(int i=0; i<nca; ++i) { int id=grid.ncoarse+(is-1)*grid.ngridmax+cr; d[i]=grid.uold(id, (mhd?8:5)+nener+is)/std::max(grid.uold(id, 1), 1e-10); cr=grid.next[cr-1]; } write_record_internal(file, d.data(), nca); }
        }
    }
    file.close();
}
void RamsesWriter::write_grav(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary | std::ios::trunc);
    if (!file.is_open()) return;
    int32_t val=grid.ncpu; write_record_internal(file, &val, 1); val=NDIM+1; write_record_internal(file, &val, 1); val=grid.nlevelmax; write_record_internal(file, &val, 1); val=0; write_record_internal(file, &val, 1);
    for(int il=1; il <= grid.nlevelmax; ++il) for(int ic=1; ic <= grid.ncpu + grid.nboundary; ++ic) {
        int nca=(ic<=grid.ncpu)?grid.numbl(ic, il):0; write_record_internal(file, &il, 1); write_record_internal(file, &nca, 1);
        if (nca > 0) for(int is=1; is<=constants::twotondim; ++is) {
            std::vector<double> pd(nca); int cr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ pd[i]=grid.phi[grid.ncoarse+(is-1)*grid.ngridmax+cr-1]; cr=grid.next[cr-1]; } write_record_internal(file, pd.data(), nca);
            for(int d=1; d<=NDIM; ++d) { std::vector<double> fd(nca); cr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ fd[i]=grid.f(grid.ncoarse+(is-1)*grid.ngridmax+cr-1, d); cr=grid.next[cr-1]; } write_record_internal(file, fd.data(), nca); }
        }
    }
    file.close();
}
void RamsesWriter::write_header(const SnapshotInfo& info) {
    std::ofstream file(filename_); if (!file.is_open()) return;
    file << "ncpu        = " << std::setw(15) << 1 << std::endl << "ndim        = " << std::setw(15) << NDIM << std::endl << "levelmin    = " << std::setw(15) << params::levelmin << std::endl << "levelmax    = " << std::setw(15) << params::nlevelmax << std::endl << "ngridmax    = " << std::setw(15) << params::ngridmax << std::endl << "nstep       = " << std::setw(15) << info.nstep << std::endl << "nstep_coarse=" << std::setw(15) << info.nstep_coarse << std::endl << "noutput     =" << std::setw(15) << info.noutput << std::endl << "time        = " << std::setw(15) << std::scientific << info.t << std::endl << "boxlen      = " << std::setw(15) << std::scientific << params::boxlen << std::endl << "omega_m     = 0.0" << std::endl << "omega_l     = 0.0" << std::endl << "omega_k     = 0.0" << std::endl << "omega_b     = 0.0" << std::endl << "h0          = 0.0" << std::endl << "aexp        = 1.0" << std::endl << "unit_l      = 1.0" << std::endl << "unit_d      = 1.0" << std::endl << "unit_t      = 1.0" << std::endl << std::endl << "ordering type=hilbert" << std::endl << std::endl;
}
void RamsesWriter::write_hydro_descriptor(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_); if (!file.is_open()) return;
    bool mhd = false;
#ifdef MHD
    mhd = true;
#endif
    int ivar = 1; file << ivar++ << ", density, double" << std::endl;
    for(int d=1; d<=NDIM; ++d) { std::string vn = "velocity_"; if(d==1) vn += "x"; else if(d==2) vn += "y"; else vn += "z"; file << ivar++ << ", " << vn << ", double" << std::endl; }
    if(mhd) { file << ivar++ << ", B_x_left, double" << std::endl << ivar++ << ", B_y_left, double" << std::endl << ivar++ << ", B_z_left, double" << std::endl << ivar++ << ", B_x_right, double" << std::endl << ivar++ << ", B_y_right, double" << std::endl << ivar++ << ", B_z_right, double" << std::endl; }
    for(int ie=1; ie<=info.nener; ++ie) { std::stringstream ss; ss << "non_thermal_pressure_" << std::setfill('0') << std::setw(2) << ie; file << ivar++ << ", " << ss.str() << ", double" << std::endl; }
    file << ivar++ << ", pressure, double" << std::endl;
    int nps = (grid.nvar > (mhd ? 8 : 5) + info.nener) ? grid.nvar - ((mhd ? 8 : 5) + info.nener) : 0;
    for(int is=1; is<=nps; ++is) file << ivar++ << ", scalar_" << is << ", double" << std::endl;
}
void RamsesWriter::write_header_file(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_); if (!file.is_open()) return;
    file << "# Generated by RAMSES-CPP" << std::endl << "ncpu " << grid.ncpu << std::endl << "ndim " << NDIM << std::endl << "nvar " << grid.nvar << std::endl << "nlevelmax " << grid.nlevelmax << std::endl << "nboundary " << grid.nboundary << std::endl << "total 0" << std::endl << "0" << std::endl;
}
}
