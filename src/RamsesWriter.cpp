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

template <typename T>
void RamsesWriter::write_record_internal(std::ofstream& file, const T* data, size_t count) {
    int32_t size_bytes = static_cast<int32_t>(count * sizeof(T));
    file.write(reinterpret_cast<const char*>(&size_bytes), sizeof(int32_t));
    if (count > 0) file.write(reinterpret_cast<const char*>(data), size_bytes);
    file.write(reinterpret_cast<const char*>(&size_bytes), sizeof(int32_t));
}

template <typename T>
void RamsesWriter::write_record(const T* data, size_t count) {
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

    int32_t val;
    val = grid.ncpu; write_record_internal(file, &val, 1);
    val = NDIM; write_record_internal(file, &val, 1);
    int32_t nxnyz[3] = {params::nx, params::ny, params::nz}; write_record_internal(file, nxnyz, 3);
    val = grid.nlevelmax; write_record_internal(file, &val, 1);
    val = grid.ngridmax; write_record_internal(file, &val, 1);
    val = 0; write_record_internal(file, &val, 1); // nboundary
    val = grid.ngridmax - grid.numbf; write_record_internal(file, &val, 1); // ngrid_current
    double boxlen = params::boxlen; write_record_internal(file, &boxlen, 1);
    
    int32_t nout_iout_ifout[3] = {info.noutput, info.iout, 1}; // ifout starts at 1
    write_record_internal(file, nout_iout_ifout, 3);
    
    std::vector<double> tout(info.noutput, 0.0);
    for(int i=0; i<std::min((int)info.tout.size(), info.noutput); ++i) tout[i] = info.tout[i];
    write_record_internal(file, tout.data(), info.noutput);
    std::vector<double> aout(info.noutput, 0.0); write_record_internal(file, aout.data(), info.noutput);
    double t = info.t; write_record_internal(file, &t, 1);
    std::vector<double> dtold(grid.nlevelmax, 0.0); write_record_internal(file, dtold.data(), grid.nlevelmax);
    std::vector<double> dtnew(grid.nlevelmax, 0.0); write_record_internal(file, dtnew.data(), grid.nlevelmax);
    int32_t nstep_vals[2] = {info.nstep, info.nstep_coarse}; write_record_internal(file, nstep_vals, 2);
    
    double einit_data[3] = {info.einit, info.mass_tot_0, info.rho_tot}; write_record_internal(file, einit_data, 3);
    double cosmo[7] = {info.omega_m, info.omega_l, info.omega_k, info.omega_b, info.h0, info.aexp_ini, info.boxlen_ini}; write_record_internal(file, cosmo, 7);
    double aexp_hexp[5] = {info.aexp, info.hexp, info.aexp_old, info.epot_tot_int, info.epot_tot_old}; write_record_internal(file, aexp_hexp, 5);
    double msph = info.mass_sph; write_record_internal(file, &msph, 1);
    
    std::vector<int32_t> hl(grid.ncpu * grid.nlevelmax), tl(grid.ncpu * grid.nlevelmax), nl(grid.ncpu * grid.nlevelmax);
    for(int il=0; il<grid.nlevelmax; ++il) for(int ic=0; ic<grid.ncpu; ++ic) {
        hl[il * grid.ncpu + ic] = grid.headl(ic+1, il+1);
        tl[il * grid.ncpu + ic] = grid.taill(ic+1, il+1);
        nl[il * grid.ncpu + ic] = grid.numbl(ic+1, il+1);
    }
    write_record_internal(file, hl.data(), hl.size()); write_record_internal(file, tl.data(), tl.size()); write_record_internal(file, nl.data(), nl.size());
    
    std::vector<int32_t> nt(10 * grid.nlevelmax, 0); 
    for(int il=0; il<grid.nlevelmax; ++il) nt[il * 10] = grid.count_grids_at_level(il+1);
    write_record_internal(file, nt.data(), nt.size());
    
    // Write free memory
    int32_t mem_info[5] = {grid.headf, grid.tailf, grid.numbf, grid.ngridmax - grid.numbf, grid.ngridmax - grid.numbf}; 
    write_record_internal(file, mem_info, 5);
    
    char ord[128] = {0}; std::strncpy(ord, "hilbert", 127); write_record_internal(file, ord, 128);
    
    std::vector<double> bound_key(grid.ncpu + 1, 0.0);
    bound_key[grid.ncpu] = 1.0;
    write_record_internal(file, bound_key.data(), bound_key.size());

    // Write coarse level
    write_record_internal(file, &grid.son[1], grid.ncoarse);
    write_record_internal(file, &grid.flag1[1], grid.ncoarse);
    write_record_internal(file, &grid.cpu_map[1], grid.ncoarse);

    for(int il=1; il <= grid.nlevelmax; ++il) for(int ic=1; ic <= grid.ncpu; ++ic) {
        int nca = grid.numbl(ic, il);
        if (nca > 0) {
            std::vector<int32_t> igs(nca), nxt(nca), prv(nca), ftr(nca);
            std::vector<double> xgd(nca);
            std::vector<int32_t> nbd(nca), sd(nca), f1d(nca), cmd(nca);
            int curr = grid.headl(ic, il);
            for(int i=0; i<nca; ++i) { igs[i]=curr; nxt[i]=grid.next[curr-1]; prv[i]=grid.prev[curr-1]; ftr[i]=grid.father[curr-1]; curr=grid.next[curr-1]; }
            write_record_internal(file, igs.data(), nca);
            write_record_internal(file, nxt.data(), nca);
            write_record_internal(file, prv.data(), nca);
            for(int d=1; d<=NDIM; ++d) { curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ xgd[i]=grid.xg[(d-1)*grid.ngridmax+(curr-1)]; curr=grid.next[curr-1]; } write_record_internal(file, xgd.data(), nca); }
            write_record_internal(file, ftr.data(), nca);
            for(int n=1; n<=2*NDIM; ++n) { curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ nbd[i]=grid.nbor[(n-1)*grid.ngridmax+(curr-1)]; curr=grid.next[curr-1]; } write_record_internal(file, nbd.data(), nca); }
            for(int s=1; s<=constants::twotondim; ++s) { curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ sd[i]=grid.son[grid.ncoarse+(s-1)*grid.ngridmax+curr]; curr=grid.next[curr-1]; } write_record_internal(file, sd.data(), nca); }
            for(int s=1; s<=constants::twotondim; ++s) { curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ cmd[i]=grid.cpu_map[grid.ncoarse+(s-1)*grid.ngridmax+curr]; curr=grid.next[curr-1]; } write_record_internal(file, cmd.data(), nca); }
            for(int s=1; s<=constants::twotondim; ++s) { curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ f1d[i]=grid.flag1[grid.ncoarse+(s-1)*grid.ngridmax+curr]; curr=grid.next[curr-1]; } write_record_internal(file, f1d.data(), nca); }
        }
    }
    file.close();
}

void RamsesWriter::write_hydro(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary | std::ios::trunc);
    if (!file.is_open()) return;

    int32_t val;
    val = grid.ncpu; write_record_internal(file, &val, 1);
    
    bool mhd = false;
#ifdef MHD
    mhd = true;
#endif

    int nener = info.nener;
    int npscal = (grid.nvar > (mhd ? 8 : 5) + nener) ? grid.nvar - ((mhd ? 8 : 5) + nener) : 0;
    int nvar_file = 1 + NDIM + nener + 1 + npscal;
    if (mhd) nvar_file += 6; 
    write_record_internal(file, &nvar_file, 1);
    
    val = NDIM; write_record_internal(file, &val, 1);
    val = grid.nlevelmax; write_record_internal(file, &val, 1);
    val = 0; write_record_internal(file, &val, 1); // nboundary
    double gamma = info.gamma; write_record_internal(file, &gamma, 1);

    real_t smallr = 1e-10;

    for(int il=1; il <= grid.nlevelmax; ++il) for(int ic=1; ic <= grid.ncpu; ++ic) {
        int nca = grid.numbl(ic, il);
        val = il; write_record_internal(file, &val, 1);
        val = nca; write_record_internal(file, &val, 1);
        if (nca > 0) {
            for(int is=1; is<=constants::twotondim; ++is) {
                std::vector<double> data(nca);
                int curr = grid.headl(ic, il);
                for(int i=0; i<nca; ++i) { data[i] = grid.uold(grid.ncoarse + (is-1)*grid.ngridmax + curr, 1); curr = grid.next[curr-1]; }
                write_record_internal(file, data.data(), nca);

                for(int d=1; d<=NDIM; ++d) {
                    curr = grid.headl(ic, il);
                    for(int i=0; i<nca; ++i) {
                        int idc = grid.ncoarse + (is-1)*grid.ngridmax + curr;
                        data[i] = grid.uold(idc, 1+d) / std::max(grid.uold(idc, 1), smallr);
                        curr = grid.next[curr-1];
                    }
                    write_record_internal(file, data.data(), nca);
                }

                if (mhd) {
                    for(int d=1; d<=3; ++d) {
                        curr = grid.headl(ic, il);
                        for(int i=0; i<nca; ++i) { data[i] = grid.uold(grid.ncoarse + (is-1)*grid.ngridmax + curr, 5+d); curr = grid.next[curr-1]; }
                        write_record_internal(file, data.data(), nca);
                    }
                    for(int d=1; d<=3; ++d) {
                        curr = grid.headl(ic, il);
                        for(int i=0; i<nca; ++i) { data[i] = grid.uold(grid.ncoarse + (is-1)*grid.ngridmax + curr, grid.nvar - 3 + d); curr = grid.next[curr-1]; }
                        write_record_internal(file, data.data(), nca);
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
                    write_record_internal(file, data.data(), nca);
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
                write_record_internal(file, data.data(), nca);

                for(int iscal=1; iscal<=npscal; ++iscal) {
                    curr = grid.headl(ic, il);
                    for(int i=0; i<nca; ++i) {
                        int idc = grid.ncoarse + (is-1)*grid.ngridmax + curr;
                        int iv_in = (mhd ? 8 : 5) + nener + iscal;
                        data[i] = grid.uold(idc, iv_in) / std::max(grid.uold(idc, 1), smallr);
                        curr = grid.next[curr-1];
                    }
                    write_record_internal(file, data.data(), nca);
                }
            }
        }
    }
    file.close();
}

void RamsesWriter::write_grav(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary | std::ios::trunc);
    if (!file.is_open()) return;

    int32_t val;
    val = grid.ncpu; write_record_internal(file, &val, 1);
    val = NDIM + 1; write_record_internal(file, &val, 1);
    val = grid.nlevelmax; write_record_internal(file, &val, 1);
    val = 0; write_record_internal(file, &val, 1);
    for(int il=1; il <= grid.nlevelmax; ++il) for(int ic=1; ic <= grid.ncpu; ++ic) {
        int nca = grid.numbl(ic, il);
        val = il; write_record_internal(file, &val, 1);
        val = nca; write_record_internal(file, &val, 1);
        if (nca > 0) for(int is=1; is<=constants::twotondim; ++is) {
            std::vector<double> pd(nca); int curr = grid.headl(ic, il);
            for(int i=0; i<nca; ++i){ pd[i]=grid.phi[grid.ncoarse+(is-1)*grid.ngridmax+curr]; curr=grid.next[curr-1]; } write_record_internal(file, pd.data(), nca);
            for(int d=1; d<=NDIM; ++d) { std::vector<double> fd(nca); curr=grid.headl(ic, il); for(int i=0; i<nca; ++i){ fd[i]=grid.f(grid.ncoarse+(is-1)*grid.ngridmax+curr, d); curr=grid.next[curr-1]; } write_record_internal(file, fd.data(), nca); }
        }
    }
    file.close();
}

void RamsesWriter::write_header(const SnapshotInfo& info) {
    std::ofstream file(filename_);
    if (!file.is_open()) return;
    file << "ncpu        = " << std::setw(15) << 1 << std::endl;
    file << "ndim        = " << std::setw(15) << NDIM << std::endl;
    file << "levelmin    = " << std::setw(15) << params::levelmin << std::endl;
    file << "levelmax    = " << std::setw(15) << params::nlevelmax << std::endl;
    file << "ngridmax    = " << std::setw(15) << params::ngridmax << std::endl;
    file << "nstep       = " << std::setw(15) << info.nstep << std::endl;
    file << "nstep_coarse=" << std::setw(15) << info.nstep_coarse << std::endl;
    file << "noutput     =" << std::setw(15) << info.noutput << std::endl;
    file << "time        = " << std::setw(15) << std::scientific << info.t << std::endl;
    file << "boxlen      = " << std::setw(15) << std::scientific << params::boxlen << std::endl;
    file << "omega_m     = " << std::setw(15) << std::scientific << info.omega_m << std::endl;
    file << "omega_l     = " << std::setw(15) << std::scientific << info.omega_l << std::endl;
    file << "omega_k     = " << std::setw(15) << std::scientific << info.omega_k << std::endl;
    file << "omega_b     = " << std::setw(15) << std::scientific << info.omega_b << std::endl;
    file << "h0          = " << std::setw(15) << std::scientific << info.h0 << std::endl;
    file << "aexp        = " << std::setw(15) << std::scientific << info.aexp << std::endl;
    file << "unit_l      = " << std::setw(15) << std::scientific << params::units_length << std::endl;
    file << "unit_d      = " << std::setw(15) << std::scientific << params::units_density << std::endl;
    file << "unit_t      = " << std::setw(15) << std::scientific << params::units_time << std::endl;
    file << std::endl;
    file << "ordering type=hilbert" << std::endl;
    file << std::endl;
}

void RamsesWriter::write_hydro_descriptor(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_);
    if (!file.is_open()) return;
    
    bool mhd = false;
#ifdef MHD
    mhd = true;
#endif

    int ivar = 1;
    file << ivar++ << ", density, double" << std::endl;
    for(int d=1; d<=NDIM; ++d) {
        std::string vname = "velocity_";
        if(d==1) vname += "x"; else if(d==2) vname += "y"; else vname += "z";
        file << ivar++ << ", " << vname << ", double" << std::endl;
    }
    if(mhd) {
        file << ivar++ << ", B_x_left, double" << std::endl;
        file << ivar++ << ", B_y_left, double" << std::endl;
        file << ivar++ << ", B_z_left, double" << std::endl;
        file << ivar++ << ", B_x_right, double" << std::endl;
        file << ivar++ << ", B_y_right, double" << std::endl;
        file << ivar++ << ", B_z_right, double" << std::endl;
    }
    for(int ie=1; ie<=info.nener; ++ie) {
        std::stringstream ss;
        ss << "non_thermal_pressure_" << std::setfill('0') << std::setw(2) << ie;
        file << ivar++ << ", " << ss.str() << ", double" << std::endl;
    }
    file << ivar++ << ", pressure, double" << std::endl;
    int npscal = (grid.nvar > (mhd ? 8 : 5) + info.nener) ? grid.nvar - ((mhd ? 8 : 5) + info.nener) : 0;
    for(int iscal=1; iscal<=npscal; ++iscal) {
        file << ivar++ << ", scalar_" << iscal << ", double" << std::endl;
    }
}

void RamsesWriter::write_extra_headers(const SnapshotInfo& info) {
    std::ofstream file(filename_);
    if (!file.is_open()) return;
    file << "# Generated by RAMSES-CPP" << std::endl;
    file << "total       0" << std::endl;
    file << "0" << std::endl;
}

} // namespace ramses
