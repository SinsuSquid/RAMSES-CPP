#include "ramses/RamsesWriter.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>

namespace ramses {

RamsesWriter::RamsesWriter(const std::string& filename) {
    file_.open(filename, std::ios::binary);
}

bool RamsesWriter::is_open() const {
    return file_.is_open();
}

template <typename T>
void RamsesWriter::write_record(const T* data, size_t count) {
    if (!file_.is_open()) return;
    int32_t size_bytes = static_cast<int32_t>(count * sizeof(T));
    file_.write(reinterpret_cast<const char*>(&size_bytes), sizeof(int32_t));
    file_.write(reinterpret_cast<const char*>(data), size_bytes);
    file_.write(reinterpret_cast<const char*>(&size_bytes), sizeof(int32_t));
}

template void RamsesWriter::write_record<int32_t>(const int32_t*, size_t);
template void RamsesWriter::write_record<double>(const double*, size_t);
template void RamsesWriter::write_record<char>(const char*, size_t);

void RamsesWriter::write_amr(const AmrGrid& grid, const SnapshotInfo& info) {
    int32_t val;
    val = grid.ncpu; write_record(&val, 1);
    val = NDIM; write_record(&val, 1);
    int32_t nxnyz[3] = {params::nx, params::ny, params::nz}; write_record(nxnyz, 3);
    val = grid.nlevelmax; write_record(&val, 1);
    val = grid.ngridmax; write_record(&val, 1);
    val = 0; write_record(&val, 1); // nboundary
    val = 0; write_record(&val, 1); // ngrid_current
    double boxlen = params::boxlen; write_record(&boxlen, 1);
    int32_t nout_iout_ifout[3] = {info.noutput, info.iout, 0}; write_record(nout_iout_ifout, 3);
    std::vector<double> tout(info.noutput, 0.0);
    for(int i=0; i<std::min((int)info.tout.size(), info.noutput); ++i) tout[i] = info.tout[i];
    write_record(tout.data(), info.noutput);
    std::vector<double> aout(info.noutput, 0.0); write_record(aout.data(), info.noutput);
    double t = info.t; write_record(&t, 1);
    std::vector<double> dtold(grid.nlevelmax, 0.0); write_record(dtold.data(), grid.nlevelmax);
    std::vector<double> dtnew(grid.nlevelmax, 0.0); write_record(dtnew.data(), grid.nlevelmax);
    int32_t nstep_coarse[2] = {info.nstep, info.nstep}; write_record(nstep_coarse, 2);
    double einit[3] = {0.0, 0.0, 0.0}; write_record(einit, 3);
    double cosmo[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; write_record(cosmo, 7);
    double aexp_hexp[5] = {1.0, 0.0, 1.0, 0.0, 0.0}; write_record(aexp_hexp, 5);
    double msph = 0.0; write_record(&msph, 1);
    
    std::vector<int32_t> hl(grid.ncpu * grid.nlevelmax), tl(grid.ncpu * grid.nlevelmax), nl(grid.ncpu * grid.nlevelmax);
    for(int il=0; il<grid.nlevelmax; ++il) for(int ic=0; ic<grid.ncpu; ++ic) {
        hl[il * grid.ncpu + ic] = grid.headl(ic+1, il+1);
        tl[il * grid.ncpu + ic] = grid.taill(ic+1, il+1);
        nl[il * grid.ncpu + ic] = grid.numbl(ic+1, il+1);
    }
    write_record(hl.data(), hl.size()); write_record(tl.data(), tl.size()); write_record(nl.data(), nl.size());
    std::vector<int32_t> nt(10 * grid.nlevelmax, 0); 
    for(int il=0; il<grid.nlevelmax; ++il) nt[il] = grid.count_grids_at_level(il+1);
    write_record(nt.data(), nt.size());
    int32_t fm[5] = {grid.headf, grid.tailf, grid.numbf, 0, 0}; write_record(fm, 5);
    char ord[128] = {0}; std::strncpy(ord, "hilbert", 127); write_record(ord, 128);
    int32_t ks = 0; write_record(&ks, 1);
    
    std::vector<int32_t> sc(grid.ncoarse), f1c(grid.ncoarse), cmc(grid.ncoarse);
    for(int i=0; i<grid.ncoarse; ++i) { sc[i]=grid.son[i+1]; f1c[i]=grid.flag1[i+1]; cmc[i]=grid.cpu_map[i+1]; }
    write_record(sc.data(), grid.ncoarse); write_record(f1c.data(), grid.ncoarse); write_record(cmc.data(), grid.ncoarse);

    for(int il=1; il <= grid.nlevelmax; ++il) for(int ic=1; ic <= grid.ncpu; ++ic) {
        int nca = grid.numbl(ic, il);
        if (nca > 0) {
            std::vector<int32_t> igs(nca), nxt(nca), prv(nca), ftr(nca);
            int curr = grid.headl(ic, il);
            for(int i=0; i<nca; ++i) { igs[i]=curr; nxt[i]=grid.next[curr-1]; prv[i]=grid.prev[curr-1]; ftr[i]=grid.father[curr-1]; curr=grid.next[curr-1]; }
            write_record(igs.data(), nca); write_record(nxt.data(), nca); write_record(prv.data(), nca);
            for(int d=1; d<=NDIM; ++d) { std::vector<double> xgd(nca); curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ xgd[i]=grid.xg[(d-1)*grid.ngridmax+(curr-1)]; curr=grid.next[curr-1]; } write_record(xgd.data(), nca); }
            write_record(ftr.data(), nca);
            for(int n=1; n<=2*NDIM; ++n) { std::vector<int32_t> nbd(nca); curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ nbd[i]=grid.nbor[(n-1)*grid.ngridmax+(curr-1)]; curr=grid.next[curr-1]; } write_record(nbd.data(), nca); }
            for(int s=1; s<=constants::twotondim; ++s) { std::vector<int32_t> sd(nca); curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ sd[i]=grid.son[grid.ncoarse+(s-1)*grid.ngridmax+curr]; curr=grid.next[curr-1]; } write_record(sd.data(), nca); }
            for(int s=1; s<=constants::twotondim; ++s) { std::vector<int32_t> f1d(nca); curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ f1d[i]=grid.flag1[grid.ncoarse+(s-1)*grid.ngridmax+curr]; curr=grid.next[curr-1]; } write_record(f1d.data(), nca); }
            for(int s=1; s<=constants::twotondim; ++s) { std::vector<int32_t> cmd(nca); curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ cmd[i]=grid.cpu_map[grid.ncoarse+(s-1)*grid.ngridmax+curr]; curr=grid.next[curr-1]; } write_record(cmd.data(), nca); }
        }
    }
    file_.close();
}

void RamsesWriter::write_hydro(const AmrGrid& grid, const SnapshotInfo& info) {
    int32_t val;
    val = grid.ncpu; write_record(&val, 1);
    
    bool mhd = false;
#ifdef MHD
    mhd = true;
#endif

    int nener = info.nener;
    int npscal = (grid.nvar > (mhd ? 8 : 5) + nener) ? grid.nvar - ((mhd ? 8 : 5) + nener) : 0;
    int nvar_file = 1 + NDIM + nener + 1 + npscal;
    if (mhd) nvar_file += 6; 
    write_record(&nvar_file, 1);
    
    val = NDIM; write_record(&val, 1);
    val = grid.nlevelmax; write_record(&val, 1);
    val = 0; write_record(&val, 1); // nboundary
    double gamma = info.gamma; write_record(&gamma, 1);

    real_t smallr = 1e-10;

    for(int il=1; il <= grid.nlevelmax; ++il) for(int ic=1; ic <= grid.ncpu; ++ic) {
        int nca = grid.numbl(ic, il);
        val = il; write_record(&val, 1);
        val = nca; write_record(&val, 1);
        if (nca > 0) {
            for(int is=1; is<=constants::twotondim; ++is) {
                std::vector<double> data(nca);
                int curr = grid.headl(ic, il);
                for(int i=0; i<nca; ++i) { data[i] = grid.uold(grid.ncoarse + (is-1)*grid.ngridmax + curr, 1); curr = grid.next[curr-1]; }
                write_record(data.data(), nca);

                for(int d=1; d<=NDIM; ++d) {
                    curr = grid.headl(ic, il);
                    for(int i=0; i<nca; ++i) {
                        int idc = grid.ncoarse + (is-1)*grid.ngridmax + curr;
                        data[i] = grid.uold(idc, 1+d) / std::max(grid.uold(idc, 1), smallr);
                        curr = grid.next[curr-1];
                    }
                    write_record(data.data(), nca);
                }

                if (mhd) {
                    for(int d=1; d<=3; ++d) {
                        curr = grid.headl(ic, il);
                        for(int i=0; i<nca; ++i) { data[i] = grid.uold(grid.ncoarse + (is-1)*grid.ngridmax + curr, 5+d); curr = grid.next[curr-1]; }
                        write_record(data.data(), nca);
                    }
                    for(int d=1; d<=3; ++d) {
                        curr = grid.headl(ic, il);
                        for(int i=0; i<nca; ++i) { data[i] = grid.uold(grid.ncoarse + (is-1)*grid.ngridmax + curr, grid.nvar - 3 + d); curr = grid.next[curr-1]; }
                        write_record(data.data(), nca);
                    }
                }

                for(int ie=1; ie<=nener; ++ie) {
                    curr = grid.headl(ic, il);
                    for(int i=0; i<nca; ++i) {
                        int idc = grid.ncoarse + (is-1)*grid.ngridmax + curr;
                        int iv_in = (mhd ? 8 : 5) + ie;
                        data[i] = (gamma - 1.0) * grid.uold(idc, iv_in);
                        curr = grid.next[curr-1];
                    }
                    write_record(data.data(), nca);
                }

                curr = grid.headl(ic, il);
                for(int i=0; i<nca; ++i) {
                    int idc = grid.ncoarse + (is-1)*grid.ngridmax + curr;
                    real_t d = std::max(grid.uold(idc, 1), smallr);
                    real_t ek = 0.5*(grid.uold(idc,2)*grid.uold(idc,2)+grid.uold(idc,3)*grid.uold(idc,3)+grid.uold(idc,4)*grid.uold(idc,4))/d;
                    real_t em = 0, er = 0;
                    if (mhd) {
                        real_t A = 0.5*(grid.uold(idc,6)+grid.uold(idc,grid.nvar-2)), B = 0.5*(grid.uold(idc,7)+grid.uold(idc,grid.nvar-1)), C = 0.5*(grid.uold(idc,8)+grid.uold(idc,grid.nvar));
                        em = 0.5*(A*A+B*B+C*C);
                    }
                    for(int ie=0; ie<nener; ++ie) er += grid.uold(idc, (mhd ? 8 : 5) + 1 + ie);
                    data[i] = (gamma - 1.0) * (grid.uold(idc, 5) - ek - em - er);
                    curr = grid.next[curr-1];
                }
                write_record(data.data(), nca);

                for(int iscal=1; iscal<=npscal; ++iscal) {
                    curr = grid.headl(ic, il);
                    for(int i=0; i<nca; ++i) {
                        int idc = grid.ncoarse + (is-1)*grid.ngridmax + curr;
                        int iv_in = (mhd ? 8 : 5) + nener + iscal;
                        data[i] = grid.uold(idc, iv_in) / std::max(grid.uold(idc, 1), smallr);
                        curr = grid.next[curr-1];
                    }
                    write_record(data.data(), nca);
                }
            }
        }
    }
    file_.close();
}

void RamsesWriter::write_grav(const AmrGrid& grid, const SnapshotInfo& info) {
    int32_t val;
    val = grid.ncpu; write_record(&val, 1);
    val = NDIM + 1; write_record(&val, 1);
    val = grid.nlevelmax; write_record(&val, 1);
    val = 0; write_record(&val, 1);
    for(int il=1; il <= grid.nlevelmax; ++il) for(int ic=1; ic <= grid.ncpu; ++ic) {
        int nca = grid.numbl(ic, il);
        val = il; write_record(&val, 1);
        val = nca; write_record(&val, 1);
        if (nca > 0) for(int is=1; is<=constants::twotondim; ++is) {
            std::vector<double> pd(nca); int curr = grid.headl(ic, il);
            for(int i=0; i<nca; ++i){ pd[i]=grid.phi[grid.ncoarse+(is-1)*grid.ngridmax+curr]; curr=grid.next[curr-1]; } write_record(pd.data(), nca);
            for(int d=1; d<=NDIM; ++d) { std::vector<double> fd(nca); curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ fd[i]=grid.f(grid.ncoarse+(is-1)*grid.ngridmax+curr, d); curr=grid.next[curr-1]; } write_record(fd.data(), nca); }
        }
    }
    file_.close();
}

} // namespace ramses
