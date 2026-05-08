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
#include <algorithm>

namespace ramses {

RamsesWriter::RamsesWriter(const std::string& filename) : filename_(filename) {}

void RamsesWriter::write_amr(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary);
    auto write_rec = [&](const void* data, size_t size) { 
        int32_t s = (int32_t)size; 
        file.write((char*)&s, 4); 
        file.write((char*)data, size); 
        file.write((char*)&s, 4); 
    };
    
    // Header - Integers (Records 1-7)
    int32_t ncpu = grid.ncpu; write_rec(&ncpu, 4);
    int32_t ndim = NDIM; write_rec(&ndim, 4);
    int32_t nxnyz[3] = {(int32_t)params::nx, (int32_t)params::ny, (int32_t)params::nz}; write_rec(nxnyz, 12);
    int32_t nlevelmax = (int32_t)params::nlevelmax; write_rec(&nlevelmax, 4);
    int32_t ngridmax = (int32_t)grid.ngridmax; write_rec(&ngridmax, 4);
    int32_t nboundary = (int32_t)grid.nboundary; write_rec(&nboundary, 4);
    int32_t ngrid_tot = 0; for(int l=1; l<=grid.nlevelmax; ++l) { for(int c=1; c<=grid.ncpu; ++c) ngrid_tot += grid.numbl(c, l); }
    write_rec(&ngrid_tot, 4);
    
    // Header - Reals (Record 8)
    double boxlen = params::boxlen; write_rec(&boxlen, 8);
    
    // Time variables (Records 9-14)
    int32_t nout_vars[3] = {(int32_t)info.noutput, (int32_t)info.iout, (int32_t)info.iout}; 
    write_rec(nout_vars, 12);
    
    std::vector<double> tout(info.noutput, info.t); write_rec(tout.data(), tout.size() * 8);
    std::vector<double> aout(info.noutput, 1.0); write_rec(aout.data(), aout.size() * 8);
    double t_sim = info.t; write_rec(&t_sim, 8);
    
    std::vector<double> dtold(params::nlevelmax, 1e-3); write_rec(dtold.data(), dtold.size() * 8);
    std::vector<double> dtnew(params::nlevelmax, 1e-3); write_rec(dtnew.data(), dtnew.size() * 8);
    
    // Step and Cosmo (Records 15-19)
    int32_t step_vars[2] = {(int32_t)info.nstep, (int32_t)info.nstep_coarse}; write_rec(step_vars, 8);
    double cosmo1[3] = {info.einit, info.mass_tot_0, info.rho_tot}; write_rec(cosmo1, 24);
    double cosmo2[7] = {info.omega_m, info.omega_l, info.omega_k, info.omega_b, info.h0, info.aexp_ini, info.boxlen_ini}; write_rec(cosmo2, 56);
    double cosmo3[5] = {info.aexp, info.hexp, info.aexp_old, info.epot_tot_int, info.epot_tot_old}; write_rec(cosmo3, 40);
    double mass_sph_val = info.mass_sph; write_rec(&mass_sph_val, 8);
    
    // Level lists (Records 20-23)
    std::vector<int32_t> headl_list(grid.ncpu * grid.nlevelmax);
    for (int l = 1; l <= grid.nlevelmax; ++l) {
        for(int c=1; c<=grid.ncpu; ++c) headl_list[(l-1)*grid.ncpu + (c-1)] = grid.get_headl(c, l);
    }
    write_rec(headl_list.data(), headl_list.size() * 4);
    
    std::vector<int32_t> taill_list(grid.ncpu * grid.nlevelmax, 0); 
    for (int l = 1; l <= grid.nlevelmax; ++l) {
        for(int c=1; c<=grid.ncpu; ++c) taill_list[(l-1)*grid.ncpu + (c-1)] = grid.taill(c, l);
    }
    write_rec(taill_list.data(), taill_list.size() * 4);
    
    std::vector<int32_t> numbl_list(grid.ncpu * grid.nlevelmax);
    for (int l = 1; l <= grid.nlevelmax; ++l) {
        for(int c=1; c<=grid.ncpu; ++c) numbl_list[(l-1)*grid.ncpu + (c-1)] = grid.numbl(c, l);
    }
    write_rec(numbl_list.data(), numbl_list.size() * 4);
    
    std::vector<int32_t> numbtot_list(10 * grid.nlevelmax, 0);
    for (int l = 1; l <= grid.nlevelmax; ++l) {
        numbtot_list[(l-1)*10] = grid.numbl(1, l);
    }
    write_rec(numbtot_list.data(), numbtot_list.size() * 4);

    // Boundary lists (Records 24-26 if nboundary > 0)
    if (grid.nboundary > 0) {
        std::vector<int32_t> headb(grid.nboundary * grid.nlevelmax, 0); write_rec(headb.data(), headb.size() * 4);
        std::vector<int32_t> tailb(grid.nboundary * grid.nlevelmax, 0); write_rec(tailb.data(), tailb.size() * 4);
        std::vector<int32_t> numbb(grid.nboundary * grid.nlevelmax, 0); write_rec(numbb.data(), numbb.size() * 4);
    }

    // Memory (Record 27)
    int32_t mem_rec[5] = { (int32_t)grid.headf, (int32_t)grid.tailf, (int32_t)grid.numbf, 0, 0 };
    write_rec(mem_rec, 20);

    // Ordering (Record 28)
    char ordering_str[128] = {0}; 
    std::strncpy(ordering_str, "hilbert", 128); 
    write_rec(ordering_str, 128);
    
    // Domain key (Record 29)
    std::vector<double> bound_key(grid.ncpu + 1, 0.0); 
    bound_key[grid.ncpu] = 1.0; 
    write_rec(bound_key.data(), (grid.ncpu + 1) * 8);

    // Coarse level (Records 30-32)
    std::vector<int32_t> coarse_son(grid.ncoarse); for(int i=0; i<grid.ncoarse; ++i) coarse_son[i] = grid.son[i];
    write_rec(coarse_son.data(), grid.ncoarse * 4);
    std::vector<int32_t> coarse_flag(grid.ncoarse); for(int i=0; i<grid.ncoarse; ++i) coarse_flag[i] = grid.flag1[i];
    write_rec(coarse_flag.data(), grid.ncoarse * 4);
    std::vector<int32_t> coarse_cpu(grid.ncoarse); for(int i=0; i<grid.ncoarse; ++i) coarse_cpu[i] = grid.cpu_map[i];
    write_rec(coarse_cpu.data(), grid.ncoarse * 4);

    // Fine levels
    for (int il = 1; il <= grid.nlevelmax; ++il) {
        for (int ib = 1; ib <= grid.nboundary + grid.ncpu; ++ib) {
            int ncache = (ib <= grid.ncpu) ? grid.numbl(ib, il) : 0;
            if (ncache > 0) {
                std::vector<int> g_list; 
                int ig = (ib <= grid.ncpu) ? grid.get_headl(ib, il) : 0; 
                while(ig > 0) { g_list.push_back(ig); ig = grid.next[ig-1]; }
                
                int n = g_list.size();
                std::vector<int32_t> tmp_i(n);
                
                // 1. Grid index
                for(int i=0; i<n; ++i) tmp_i[i] = g_list[i]; write_rec(tmp_i.data(), n * 4);
                // 2. Next
                for(int i=0; i<n; ++i) tmp_i[i] = grid.next[g_list[i]-1]; write_rec(tmp_i.data(), n * 4);
                // 3. Prev
                for(int i=0; i<n; ++i) tmp_i[i] = grid.prev[g_list[i]-1]; write_rec(tmp_i.data(), n * 4);
                
                // 4. Grid center (xg)
                std::vector<double> tmp_d(n);
                for(int d=0; d<ndim; ++d) { 
                    for(int i=0; i<n; ++i) tmp_d[i] = grid.xg[d*grid.ngridmax + g_list[i]-1]; 
                    write_rec(tmp_d.data(), n * 8); 
                }
                
                // 5. Father index
                for(int i=0; i<n; ++i) tmp_i[i] = grid.father[g_list[i]-1]; write_rec(tmp_i.data(), n * 4);
                
                // 6. Neighbor indices (twondim records)
                for(int d=0; d<ndim*2; ++d) { 
                    for(int i=0; i<n; ++i) tmp_i[i] = grid.nbor[d*grid.ngridmax + g_list[i]-1]; 
                    write_rec(tmp_i.data(), n * 4); 
                }
                
                // 7. Son indices (twotondim records)
                for(int d=0; d<constants::twotondim; ++d) { 
                    for(int i=0; i<n; ++i) tmp_i[i] = grid.son[grid.ncoarse + d*grid.ngridmax + g_list[i]-1]; 
                    write_rec(tmp_i.data(), n * 4); 
                }
                
                // 8. CPU map (twotondim records)
                for(int d=0; d<constants::twotondim; ++d) { 
                    for(int i=0; i<n; ++i) tmp_i[i] = grid.cpu_map[grid.ncoarse + d*grid.ngridmax + g_list[i]-1]; 
                    write_rec(tmp_i.data(), n * 4); 
                }
                
                // 9. Flag map (twotondim records)
                for(int d=0; d<constants::twotondim; ++d) { 
                    for(int i=0; i<n; ++i) tmp_i[i] = grid.flag1[grid.ncoarse + d*grid.ngridmax + g_list[i]-1]; 
                    write_rec(tmp_i.data(), n * 4); 
                }
            }
        }
    }
}

void RamsesWriter::write_hydro(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary);
    auto write_rec = [&](const void* data, size_t size) { 
        int32_t s = (int32_t)size; 
        file.write((char*)&s, 4); 
        file.write((char*)data, size); 
        file.write((char*)&s, 4); 
    };
    
    int32_t ncpu = (int32_t)grid.ncpu; write_rec(&ncpu, 4); 
    int32_t nvar = (int32_t)grid.nvar; write_rec(&nvar, 4);
    int32_t ngrid_tot = 0; for(int l=1; l<=grid.nlevelmax; ++l) { for(int c=1; c<=grid.ncpu; ++c) ngrid_tot += grid.numbl(c, l); }
    write_rec(&ngrid_tot, 4);
    int32_t nlevelmax = (int32_t)grid.nlevelmax; write_rec(&nlevelmax, 4);
    int32_t nboundary = (int32_t)grid.nboundary; write_rec(&nboundary, 4);
    double gamma_val = grid.gamma; write_rec(&gamma_val, 8);

    // Coarse level hydro
    for (int iv = 1; iv <= grid.nvar; ++iv) {
        std::vector<double> tmp(grid.ncoarse);
        for (int i = 1; i <= grid.ncoarse; ++i) tmp[i-1] = grid.uold(i, iv);
        write_rec(tmp.data(), grid.ncoarse * 8);
    }

    for (int il = 1; il <= grid.nlevelmax; ++il) {
        for (int icpu = 1; icpu <= grid.ncpu + grid.nboundary; ++icpu) {
            int ncache = (icpu <= grid.ncpu) ? grid.numbl(icpu, il) : 0;
            int32_t il_rec = il, nc_rec = ncache;
            write_rec(&il_rec, 4);
            write_rec(&nc_rec, 4);
            if (ncache > 0) {
                std::vector<int> g_list; 
                int ig = (icpu <= grid.ncpu) ? grid.get_headl(icpu, il) : 0; 
                while(ig > 0) { g_list.push_back(ig); ig = grid.next[ig-1]; }
                
                for (int ic = 0; ic < constants::twotondim; ++ic) {
                    for(int iv=1; iv<=grid.nvar; ++iv) {
                        std::vector<double> tmp(ncache);
                        for(int i=0; i<ncache; ++i) tmp[i] = grid.uold(grid.ncoarse + ic*grid.ngridmax + g_list[i], iv);
                        write_rec(tmp.data(), ncache * 8);
                    }
                }
            }
        }
    }
}

void RamsesWriter::write_rt(const AmrGrid& grid, const SnapshotInfo& info, int nGroups, real_t rt_c) {
    std::ofstream file(filename_, std::ios::binary);
    auto write_rec = [&](const void* data, size_t size) { 
        int32_t s = (int32_t)size; 
        file.write((char*)&s, 4); 
        file.write((char*)data, size); 
        file.write((char*)&s, 4); 
    };
    
    int nrtvar_fields = nGroups * (1 + NDIM);
    int32_t ncpu = grid.ncpu; write_rec(&ncpu, 4);
    int32_t nrv_rec = (int32_t)nrtvar_fields; write_rec(&nrv_rec, 4);
    int32_t ndim = (int32_t)NDIM; write_rec(&ndim, 4);
    int32_t nlevelmax = (int32_t)params::nlevelmax; write_rec(&nlevelmax, 4);
    int32_t nboundary = (int32_t)grid.nboundary; write_rec(&nboundary, 4);
    double gamma_val = grid.gamma; write_rec(&gamma_val, 8);

    int nvar_hydro_base = 5 + info.nener;
#ifdef MHD
    nvar_hydro_base = 11 + info.nener;
#endif
    int nrtvar_total = grid.nvar - nvar_hydro_base;
    int nIons = nrtvar_total - nrtvar_fields;
    int rt_start_var = nvar_hydro_base + nIons;

    // Coarse level RT
    for (int igr = 0; igr < nGroups; ++igr) {
        int iv_n = rt_start_var + (igr * (1 + NDIM)) + 1;
        std::vector<double> tmp_n(grid.ncoarse);
        for (int i = 1; i <= grid.ncoarse; ++i) tmp_n[i-1] = grid.uold(i, iv_n) * rt_c;
        write_rec(tmp_n.data(), grid.ncoarse * 8);
        for(int d=0; d<NDIM; ++d) {
            std::vector<double> tmp_f(grid.ncoarse);
            int iv_f = iv_n + 1 + d;
            for (int i = 1; i <= grid.ncoarse; ++i) tmp_f[i-1] = grid.uold(i, iv_f);
            write_rec(tmp_f.data(), grid.ncoarse * 8);
        }
    }

    for (int il = 1; il <= grid.nlevelmax; ++il) {
        for (int icpu = 1; icpu <= grid.ncpu + grid.nboundary; ++icpu) {
            int ncache = (icpu <= grid.ncpu) ? grid.numbl(icpu, il) : 0;
            int32_t il_rec = il, nc_rec = ncache;
            write_rec(&il_rec, 4);
            write_rec(&nc_rec, 4);
            if (ncache > 0) {
                std::vector<int> g_list; 
                int ig = (icpu <= grid.ncpu) ? grid.get_headl(icpu, il) : 0; 
                while(ig > 0) { g_list.push_back(ig); ig = grid.next[ig-1]; }
                
                for (int ic = 0; ic < constants::twotondim; ++ic) {
                    for (int igr = 0; igr < nGroups; ++igr) {
                        // Np
                        std::vector<double> tmp_n(ncache);
                        int iv_n = rt_start_var + (igr * (1 + NDIM)) + 1;
                        for (int i = 0; i < ncache; ++i) tmp_n[i] = grid.uold(grid.ncoarse + ic * grid.ngridmax + g_list[i], iv_n) * rt_c;
                        write_rec(tmp_n.data(), ncache * 8);
                        
                        // Fp
                        for(int d=0; d<NDIM; ++d) {
                            std::vector<double> tmp_f(ncache);
                            int iv_f = iv_n + 1 + d;
                            for (int i = 0; i < ncache; ++i) tmp_f[i] = grid.uold(grid.ncoarse + ic * grid.ngridmax + g_list[i], iv_f);
                            write_rec(tmp_f.data(), ncache * 8);
                        }
                    }
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
    file << "nstep_coarse= " << info.nstep_coarse << std::endl;
    file << "noutput     = " << info.noutput << std::endl;
    file << "time        = " << std::scientific << std::setprecision(16) << info.t << std::endl;
    file << "boxlen      = " << params::boxlen << std::endl;
    file << "omega_m     = " << info.omega_m << std::endl;
    file << "omega_l     = " << info.omega_l << std::endl;
    file << "omega_k     = " << info.omega_k << std::endl;
    file << "omega_b     = " << info.omega_b << std::endl;
    file << "h0          = " << info.h0 << std::endl;
    file << "aexp        = " << info.aexp << std::endl;
    file << "unit_l      = 1.0" << std::endl;
    file << "unit_d      = 1.0" << std::endl;
    file << "unit_t      = 1.0" << std::endl;
    file << "ordering type= hilbert" << std::endl;
}

void RamsesWriter::write_header_file(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_);
    file << "# Family Count" << std::endl;
    file << "undefined 0" << std::endl;
    file << "Particle fields" << std::endl;
}

void RamsesWriter::write_hydro_descriptor(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_);
    file << "# version: 1" << std::endl;
    int ivar = 1;
    file << ivar++ << ", density, double" << std::endl;
    file << ivar++ << ", momentum_x, double" << std::endl;
    if (NDIM > 1) file << ivar++ << ", momentum_y, double" << std::endl;
    if (NDIM > 2) file << ivar++ << ", momentum_z, double" << std::endl;
    file << ivar++ << ", total_energy, double" << std::endl;
    for (int ie = 1; ie <= info.nener; ++ie) { 
        std::stringstream ss; ss << "non_thermal_energy_" << std::setfill('0') << std::setw(2) << (ie - 1); 
        file << ivar++ << ", " << ss.str() << ", double" << std::endl; 
    }
    
    int nvar_hydro_base = 5 + info.nener;
#ifdef MHD
    nvar_hydro_base = 11 + info.nener;
#endif
    int nrtvar_all = grid.nvar - nvar_hydro_base;
    for (int ip = 1; ip <= nrtvar_all; ++ip) {
        std::stringstream ss; 
        if (ip == 1) ss << "xHII";
        else if (ip == 2) ss << "xHeII";
        else if (ip == 3) ss << "xHeIII";
        else ss << "scalar_" << std::setfill('0') << std::setw(2) << (ip - 1);
        file << ivar++ << ", " << ss.str() << ", double" << std::endl;
    }
}

void RamsesWriter::write_rt_descriptor(const AmrGrid& grid, const SnapshotInfo& info, int nGroups) {
    std::ofstream file(filename_);
    file << "# version: 1" << std::endl;
    int ivar = 1;
    for (int ig = 1; ig <= nGroups; ++ig) {
        file << ivar++ << ", photon_density_" << std::setfill('0') << std::setw(2) << ig << ", double" << std::endl;
        file << ivar++ << ", photon_flux_" << std::setfill('0') << std::setw(2) << ig << "_x, double" << std::endl;
        if (NDIM > 1) file << ivar++ << ", photon_flux_" << std::setfill('0') << std::setw(2) << ig << "_y, double" << std::endl;
        if (NDIM > 2) file << ivar++ << ", photon_flux_" << std::setfill('0') << std::setw(2) << ig << "_z, double" << std::endl;
    }
}

void RamsesWriter::write_particles(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary);
    if (!file.is_open()) return;
    int32_t ncpu = grid.ncpu, ndim = NDIM, npart = grid.npart;
    write_record_internal(file, &ncpu, 1);
    write_record_internal(file, &ndim, 1);
    write_record_internal(file, &npart, 1);
    
    int32_t seed = 12345; write_record_internal(file, &seed, 1);
    int32_t nstar = 0; real_t mstar = 0.0;
    write_record_internal(file, &nstar, 1);
    write_record_internal(file, &mstar, 1);
    write_record_internal(file, &mstar, 1);
    int32_t nsink = 0; write_record_internal(file, &nsink, 1);

    if (npart > 0) {
        for (int d = 0; d < NDIM; ++d) {
            std::vector<double> buf(npart);
            for (int i = 0; i < npart; ++i) buf[i] = grid.xp[d * grid.npartmax + i];
            write_record_internal(file, buf.data(), npart);
        }
        for (int d = 0; d < NDIM; ++d) {
            std::vector<double> buf(npart);
            for (int i = 0; i < npart; ++i) buf[i] = grid.vp[d * grid.npartmax + i];
            write_record_internal(file, buf.data(), npart);
        }
        {
            std::vector<double> buf(npart);
            for (int i = 0; i < npart; ++i) buf[i] = grid.mp[i];
            write_record_internal(file, buf.data(), npart);
        }
        {
            std::vector<int32_t> buf(npart);
            for (int i = 0; i < npart; ++i) buf[i] = grid.idp[i];
            write_record_internal(file, buf.data(), npart);
        }
        {
            std::vector<int32_t> buf(npart);
            for (int i = 0; i < npart; ++i) buf[i] = grid.levelp[i];
            write_record_internal(file, buf.data(), npart);
        }
    }
}

void RamsesWriter::write_particles_descriptor(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_);
    if (!file.is_open()) return;
    file << "# version: 1" << std::endl;
    int ivar = 1;
    for(int d=1; d<=NDIM; ++d) file << ivar++ << ", position_" << (d==1?'x':(d==2?'y':'z')) << ", double" << std::endl;
    for(int d=1; d<=NDIM; ++d) file << ivar++ << ", velocity_" << (d==1?'x':(d==2?'y':'z')) << ", double" << std::endl;
    file << ivar++ << ", mass, double" << std::endl;
    file << ivar++ << ", identity, int" << std::endl;
    file << ivar++ << ", level, int" << std::endl;
}

void RamsesWriter::write_grav(const AmrGrid&, const SnapshotInfo&) {}

} // namespace ramses
