#include "ramses/RamsesWriter.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Parameters.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <sstream>
#include <cstring>
#include <sys/stat.h>
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

    // 1-8. Basic Header
    int32_t ncpu = (int32_t)grid.ncpu; write_rec(&ncpu, 4);
    int32_t ndim = (int32_t)NDIM; write_rec(&ndim, 4);
    int32_t nx_y_z[3] = {(int32_t)grid.nx, (int32_t)grid.ny, (int32_t)grid.nz}; write_rec(nx_y_z, 12);
    int32_t nlevelmax = (int32_t)grid.nlevelmax; write_rec(&nlevelmax, 4);
    int32_t ngridmax = (int32_t)grid.ngridmax; write_rec(&ngridmax, 4);
    int32_t nboundary = (int32_t)grid.nboundary; write_rec(&nboundary, 4);
    int32_t ngrid_tot = 0; for(int l=1; l<=grid.nlevelmax; ++l) { for(int c=1; c<=grid.ncpu; ++c) ngrid_tot += grid.numbl(c, l); }
    write_rec(&ngrid_tot, 4);
    double boxlen_val = params::boxlen / (double)grid.nx; write_rec(&boxlen_val, 8);
    
    // 9-14. Time and Units
    int32_t nout_rec[3] = {(int32_t)info.noutput, (int32_t)info.iout, (int32_t)info.iout};
    write_rec(nout_rec, 12);
    std::vector<double> tout_v(info.noutput, 0.0); write_rec(tout_v.data(), tout_v.size() * 8);
    std::vector<double> aout_v(info.noutput, 1.0); write_rec(aout_v.data(), aout_v.size() * 8);
    double time_val = info.t; write_rec(&time_val, 8);
    std::vector<double> dtold(grid.nlevelmax, 0.0); write_rec(dtold.data(), dtold.size() * 8);
    std::vector<double> dtnew(grid.nlevelmax, 0.0); write_rec(dtnew.data(), dtnew.size() * 8);
    
    // 15-19. Cosmo/Mass
    int32_t step_vars[2] = {(int32_t)info.nstep, (int32_t)info.nstep_coarse}; write_rec(step_vars, 8);
    double cosmo1[3] = {info.einit, info.mass_tot_0, info.rho_tot}; write_rec(cosmo1, 24);
    double cosmo2[7] = {info.omega_m, info.omega_l, info.omega_k, info.omega_b, info.h0, info.aexp_ini, info.boxlen_ini}; write_rec(cosmo2, 56);
    double cosmo3[5] = {info.aexp, info.hexp, info.aexp_old, info.epot_tot_int, info.epot_tot_old}; write_rec(cosmo3, 40);
    double mass_sph_val = info.mass_sph; write_rec(&mass_sph_val, 8);
    
    // 20-23. Level lists
    std::vector<int32_t> headl_list(grid.ncpu * grid.nlevelmax);
    for (int l = 1; l <= grid.nlevelmax; ++l) {
        for(int c=1; c<=grid.ncpu; ++c) headl_list[(l-1)*grid.ncpu + (c-1)] = grid.get_headl(c, l);
    }
    write_rec(headl_list.data(), headl_list.size() * 4);
    
    std::vector<int32_t> taill_list(grid.ncpu * grid.nlevelmax); 
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

    // 24-26. Boundary (optional)
    if (grid.nboundary > 0) {
        std::vector<int32_t> headb(grid.nboundary * grid.nlevelmax, 0); write_rec(headb.data(), headb.size() * 4);
        std::vector<int32_t> tailb(grid.nboundary * grid.nlevelmax, 0); write_rec(tailb.data(), tailb.size() * 4);
        std::vector<int32_t> numbb(grid.nboundary * grid.nlevelmax, 0); write_rec(numbb.data(), numbb.size() * 4);
    }

    // 27-29. Memory and Domain
    int32_t mem_rec[5] = { (int32_t)grid.headf, (int32_t)grid.tailf, (int32_t)grid.numbf, 0, 0 };
    write_rec(mem_rec, 20);
    char ordering_str[128] = {0}; std::strncpy(ordering_str, "hilbert", 128); write_rec(ordering_str, 128);
    std::vector<double> bound_key(grid.ncpu + 1, 0.0); bound_key[grid.ncpu] = 1.0; write_rec(bound_key.data(), (grid.ncpu + 1) * 8);

    // 30-32. Coarse Level
    std::vector<int32_t> coarse_son(grid.ncoarse); for(int i=0; i<grid.ncoarse; ++i) coarse_son[i] = grid.son[i];
    write_rec(coarse_son.data(), grid.ncoarse * 4);
    std::vector<int32_t> coarse_flag(grid.ncoarse); for(int i=0; i<grid.ncoarse; ++i) coarse_flag[i] = grid.flag1[i];
    write_rec(coarse_flag.data(), grid.ncoarse * 4);
    std::vector<int32_t> coarse_cpu(grid.ncoarse); for(int i=0; i<grid.ncoarse; ++i) coarse_cpu[i] = grid.cpu_map[i];
    write_rec(coarse_cpu.data(), grid.ncoarse * 4);

    // Grid data
    for (int il = 1; il <= grid.nlevelmax; ++il) {
        for (int ib = 1; ib <= grid.nboundary + grid.ncpu; ++ib) {
            int ncache = (ib <= grid.ncpu) ? grid.numbl(ib, il) : 0;
            int32_t il_rec = il, nc_rec = ncache;
            write_rec(&il_rec, 4); write_rec(&nc_rec, 4);

            if (ncache > 0) {
                std::vector<int> g_list; 
                int ig = (ib <= grid.ncpu) ? grid.get_headl(ib, il) : 0; 
                while(ig > 0) { g_list.push_back(ig); ig = grid.next[ig-1]; }
                int n = g_list.size();
                std::vector<int32_t> tmp_i(n);
                
                for(int i=0; i<n; ++i) tmp_i[i] = g_list[i]; write_rec(tmp_i.data(), n * 4);
                for(int i=0; i<n; ++i) tmp_i[i] = grid.next[g_list[i]-1]; write_rec(tmp_i.data(), n * 4);
                for(int i=0; i<n; ++i) tmp_i[i] = grid.prev[g_list[i]-1]; write_rec(tmp_i.data(), n * 4);
                
                std::vector<double> tmp_d(n);
                for(int d=0; d<NDIM; ++d) { 
                    for(int i=0; i<n; ++i) tmp_d[i] = grid.xg[d*grid.ngridmax + g_list[i]-1]; 
                    write_rec(tmp_d.data(), n * 8); 
                }
                for(int i=0; i<n; ++i) tmp_i[i] = grid.father[g_list[i]-1]; write_rec(tmp_i.data(), n * 4);
                for(int d=0; d<NDIM*2; ++d) { 
                    for(int i=0; i<n; ++i) tmp_i[i] = grid.nbor[d*grid.ngridmax + g_list[i]-1]; 
                    write_rec(tmp_i.data(), n * 4); 
                }
                for(int d=0; d<constants::twotondim; ++d) { 
                    int iskip = grid.ncoarse + d*grid.ngridmax;
                    for(int i=0; i<n; ++i) tmp_i[i] = grid.son[iskip + g_list[i]-1]; 
                    write_rec(tmp_i.data(), n * 4); 
                }
                for(int d=0; d<constants::twotondim; ++d) { 
                    int iskip = grid.ncoarse + d*grid.ngridmax;
                    for(int i=0; i<n; ++i) tmp_i[i] = grid.cpu_map[iskip + g_list[i]-1]; 
                    write_rec(tmp_i.data(), n * 4); 
                }
                for(int d=0; d<constants::twotondim; ++d) { 
                    int iskip = grid.ncoarse + d*grid.ngridmax;
                    for(int i=0; i<n; ++i) tmp_i[i] = grid.flag1[iskip + g_list[i]-1]; 
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
    int32_t ndim = (int32_t)NDIM; write_rec(&ndim, 4);
    int32_t nlevelmax = (int32_t)grid.nlevelmax; write_rec(&nlevelmax, 4);
    int32_t nboundary = (int32_t)grid.nboundary; write_rec(&nboundary, 4);
    double gamma_val = grid.gamma; write_rec(&gamma_val, 8);

    real_t smallr = 1e-10;
    
    // Fine levels
    for (int il = 1; il <= grid.nlevelmax; ++il) {
        for (int icpu = 1; icpu <= grid.ncpu + grid.nboundary; ++icpu) {
            int ncache = (icpu <= grid.ncpu) ? grid.numbl(icpu, il) : 0;
            int32_t il_rec = il, nc_rec = ncache;
            write_rec(&il_rec, 4); write_rec(&nc_rec, 4);
            
            if (ncache > 0) {
                std::vector<int> g_list; 
                int ig = (icpu <= grid.ncpu) ? grid.get_headl(icpu, il) : 0; 
                while(ig > 0) { g_list.push_back(ig); ig = grid.next[ig-1]; }
                
                for (int ic = 0; ic < constants::twotondim; ++ic) {
                    for(int iv=1; iv<=grid.nvar; ++iv) {
                        std::vector<double> tmp(ncache);
                        for(int i=0; i<ncache; ++i) {
                            int idc = grid.ncoarse + ic*grid.ngridmax + g_list[i];
                            real_t d = std::max(grid.uold(idc, 1), smallr);
                            if (iv == 1) tmp[i] = d;
                            else if (iv >= 2 && iv <= 1 + NDIM) tmp[i] = grid.uold(idc, iv) / d; 
                            else if (iv == NDIM + 2) { 
                                real_t v2 = 0; for (int j = 1; j <= NDIM; ++j) { real_t v = grid.uold(idc, 1+j)/d; v2 += v*v; }
                                tmp[i] = std::max((grid.uold(idc, iv) - 0.5 * d * v2) * (grid.gamma - 1.0), 0.0);
                            }
                            else tmp[i] = grid.uold(idc, iv);
                        }
                        write_rec(tmp.data(), ncache * 8);
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
    file << "nx          = " << grid.nx << std::endl;
    file << "ny          = " << grid.ny << std::endl;
    file << "nz          = " << grid.nz << std::endl;
    file << "ncoarse     = " << grid.ncoarse << std::endl;
    file << "nboundary   = " << grid.nboundary << std::endl;
    file << "levelmin    = " << params::levelmin << std::endl;
    file << "levelmax    = " << params::nlevelmax << std::endl;
    file << "ngridmax    = " << grid.ngridmax << std::endl;
    file << "nstep       = " << info.nstep << std::endl;
    file << "nstep_coarse= " << info.nstep_coarse << std::endl;
    file << "noutput     = " << info.noutput << std::endl;
    file << "time        = " << std::scientific << std::setprecision(16) << info.t << std::endl;
    file << "boxlen      = " << std::scientific << std::setprecision(16) << params::boxlen / (double)grid.nx << std::endl;
    file << "omega_m     = " << info.omega_m << std::endl;
    file << "omega_l     = " << info.omega_l << std::endl;
    file << "omega_k     = " << info.omega_k << std::endl;
    file << "omega_b     = " << info.omega_b << std::endl;
    file << "h0          = " << info.h0 << std::endl;
    file << "aexp        = " << info.aexp << std::endl;
    file << "unit_l      = " << params::units_length << std::endl;
    file << "unit_d      = " << params::units_density << std::endl;
    file << "unit_t      = " << params::units_time << std::endl;
    file << "ordering type= hilbert" << std::endl;
}

void RamsesWriter::write_header_file(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_);
    file << "# Family Count" << std::endl;
    file << "undefined " << grid.npart << std::endl;
    file << "Particle fields" << std::endl;
}

void RamsesWriter::write_hydro_descriptor(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_);
    file << "# version: 1" << std::endl;
    int ivar = 1;
    file << ivar++ << ", density, double" << std::endl;
    file << ivar++ << ", velocity_x, double" << std::endl;
    if (NDIM > 1) file << ivar++ << ", velocity_y, double" << std::endl;
    if (NDIM > 2) file << ivar++ << ", velocity_z, double" << std::endl;
    file << ivar++ << ", pressure, double" << std::endl;
    for (int ie = 1; ie <= info.nener; ++ie) { 
        std::stringstream ss; ss << "non_thermal_pressure_" << std::setfill('0') << std::setw(2) << ie; 
        file << ivar++ << ", " << ss.str() << ", double" << std::endl; 
    }
    int nvar_hydro_base = NDIM + 2 + info.nener;
#ifdef MHD
    nvar_hydro_base = 11 + info.nener;
#endif
    int npassive_all = grid.nvar - nvar_hydro_base;
    for (int ip = 1; ip <= npassive_all; ++ip) {
        std::stringstream ss; ss << "scalar_" << std::setfill('0') << std::setw(2) << ip;
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

void RamsesWriter::write_rt(const AmrGrid& grid, const SnapshotInfo& info, int nGroups, int nIons, real_t rt_c) {
    std::ofstream file(filename_, std::ios::binary);
    auto write_rec = [&](const void* data, size_t size) {
        int32_t s = (int32_t)size;
        file.write((char*)&s, 4);
        file.write((char*)data, size);
        file.write((char*)&s, 4);
    };

    int32_t ncpu = (int32_t)grid.ncpu; write_rec(&ncpu, 4);
    int32_t nrtvar = (int32_t)nGroups * (1 + NDIM); write_rec(&nrtvar, 4);
    int32_t ndim = (int32_t)NDIM; write_rec(&ndim, 4);
    int32_t nlevelmax = (int32_t)grid.nlevelmax; write_rec(&nlevelmax, 4);
    int32_t nboundary = (int32_t)grid.nboundary; write_rec(&nboundary, 4);
    double gamma_val = grid.gamma; write_rec(&gamma_val, 8);

    int nvar_hydro = grid.nvar - nIons - nGroups * (1 + NDIM);
    int ioffset = nvar_hydro + nIons;

    for (int il = 1; il <= grid.nlevelmax; ++il) {
        for (int icpu = 1; icpu <= grid.ncpu + grid.nboundary; ++icpu) {
            int ncache = (icpu <= grid.ncpu) ? grid.numbl(icpu, il) : 0;
            int32_t il_rec = il, nc_rec = ncache;
            write_rec(&il_rec, 4); write_rec(&nc_rec, 4);

            if (ncache > 0) {
                std::vector<int> g_list;
                int ig_head = (icpu <= grid.ncpu) ? grid.get_headl(icpu, il) : 0;
                int ig_curr = ig_head;
                while(ig_curr > 0) { g_list.push_back(ig_curr); ig_curr = grid.next[ig_curr-1]; }

                for (int ic = 0; ic < constants::twotondim; ++ic) {
                    for (int igr = 1; igr <= nGroups; ++igr) {
                        // Photon density (convert to flux units: density * c)
                        std::vector<double> tmp_n(ncache);
                        int iv_n = ioffset + (igr-1)*(1+NDIM) + 1;
                        for(int i=0; i<ncache; ++i) {
                            int idc = grid.ncoarse + ic*grid.ngridmax + g_list[i];
                            tmp_n[i] = grid.uold(idc, iv_n) * rt_c;
                        }
                        write_rec(tmp_n.data(), ncache * 8);

                        // Photon flux
                        for (int d = 1; d <= NDIM; ++d) {
                            std::vector<double> tmp_f(ncache);
                            int iv_f = ioffset + (igr-1)*(1+NDIM) + 1 + d;
                            for(int i=0; i<ncache; ++i) {
                                int idc = grid.ncoarse + ic*grid.ngridmax + g_list[i];
                                tmp_f[i] = grid.uold(idc, iv_f);
                            }
                            write_rec(tmp_f.data(), ncache * 8);
                        }
                    }
                }
            }
        }
    }
}

void RamsesWriter::write_particles(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary);
    auto write_rec = [&](const void* data, size_t size) { 
        int32_t s = (int32_t)size; 
        file.write((char*)&s, 4); 
        file.write((char*)data, size); 
        file.write((char*)&s, 4); 
    };
    int32_t ncpu = grid.ncpu, ndim = NDIM, npart = grid.npart;
    write_rec(&ncpu, 4); write_rec(&ndim, 4); write_rec(&npart, 4);
    int32_t seed = 12345; write_rec(&seed, 4);
    int32_t nstar = 0; double mstar = 0.0;
    write_rec(&nstar, 4); write_rec(&mstar, 8); write_rec(&mstar, 8);
    int32_t nsink = 0; write_rec(&nsink, 4);
    if (npart > 0) {
        std::vector<int> active_indices;
        for (int i = 0; i < grid.npartmax; ++i) {
            if (grid.idp[i] > 0) active_indices.push_back(i);
        }
        // Safety check: count might not match if something is inconsistent
        int actual_npart = active_indices.size();

        for (int d = 0; d < NDIM; ++d) {
            std::vector<double> buf(actual_npart);
            for (int i = 0; i < actual_npart; ++i) buf[i] = grid.xp[d * grid.npartmax + active_indices[i]];
            write_rec(buf.data(), actual_npart * 8);
        }
        for (int d = 0; d < NDIM; ++d) {
            std::vector<double> buf(actual_npart);
            for (int i = 0; i < actual_npart; ++i) buf[i] = grid.vp[d * grid.npartmax + active_indices[i]];
            write_rec(buf.data(), actual_npart * 8);
        }
        std::vector<double> mp_buf(actual_npart); for (int i = 0; i < actual_npart; ++i) mp_buf[i] = grid.mp[active_indices[i]]; write_rec(mp_buf.data(), actual_npart * 8);
        std::vector<int32_t> idp_buf(actual_npart); for (int i = 0; i < actual_npart; ++i) idp_buf[i] = (int32_t)grid.idp[active_indices[i]]; write_rec(idp_buf.data(), actual_npart * 4);
        std::vector<int32_t> levp_buf(actual_npart); for (int i = 0; i < actual_npart; ++i) levp_buf[i] = (int32_t)grid.levelp[active_indices[i]]; write_rec(levp_buf.data(), actual_npart * 4);
        std::vector<int32_t> fam_buf(actual_npart); for (int i = 0; i < actual_npart; ++i) fam_buf[i] = (int32_t)grid.family[active_indices[i]]; write_rec(fam_buf.data(), actual_npart * 4);
        std::vector<int32_t> tag_buf(actual_npart); for (int i = 0; i < actual_npart; ++i) tag_buf[i] = (int32_t)grid.tag[active_indices[i]]; write_rec(tag_buf.data(), actual_npart * 4);
        std::vector<double> tp_buf(actual_npart); for (int i = 0; i < actual_npart; ++i) tp_buf[i] = grid.tp[active_indices[i]]; write_rec(tp_buf.data(), actual_npart * 8);
        std::vector<double> zp_buf(actual_npart); for (int i = 0; i < actual_npart; ++i) zp_buf[i] = grid.zp[active_indices[i]]; write_rec(zp_buf.data(), actual_npart * 8);
    }
}

void RamsesWriter::write_particles_descriptor(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_);
    file << "# version: 1" << std::endl;
    int ivar = 1;
    for(int d=1; d<=NDIM; ++d) file << ivar++ << ", position_" << (d==1?'x':(d==2?'y':'z')) << ", double" << std::endl;
    for(int d=1; d<=NDIM; ++d) file << ivar++ << ", velocity_" << (d==1?'x':(d==2?'y':'z')) << ", double" << std::endl;
    file << ivar++ << ", mass, double" << std::endl;
    file << ivar++ << ", identity, int" << std::endl;
    file << ivar++ << ", level, int" << std::endl;
    file << ivar++ << ", family, int" << std::endl;
    file << ivar++ << ", tag, int" << std::endl;
    file << ivar++ << ", birth_time, double" << std::endl;
    file << ivar++ << ", metallicity, double" << std::endl;
}

void RamsesWriter::write_grav(const AmrGrid& grid, const SnapshotInfo& info) {
    std::ofstream file(filename_, std::ios::binary);
    auto write_rec = [&](const void* data, size_t size) { 
        int32_t s = (int32_t)size; 
        file.write((char*)&s, 4); 
        file.write((char*)data, size); 
        file.write((char*)&s, 4); 
    };
    int32_t ncpu = (int32_t)grid.ncpu, nvar_grav = NDIM + 1;
    write_rec(&ncpu, 4); write_rec(&nvar_grav, 4);
    int32_t ngrid_tot = 0; for(int l=1; l<=grid.nlevelmax; ++l) { for(int c=1; c<=grid.ncpu; ++c) ngrid_tot += grid.numbl(c, l); }
    write_rec(&ngrid_tot, 4);
    int32_t nlevelmax = (int32_t)grid.nlevelmax, nboundary = (int32_t)grid.nboundary;
    write_rec(&nlevelmax, 4); write_rec(&nboundary, 4);
    { std::vector<double> tmp(grid.ncoarse); for (int i = 1; i <= grid.ncoarse; ++i) tmp[i-1] = grid.phi[i-1]; write_rec(tmp.data(), grid.ncoarse * 8); }
    for (int d = 1; d <= NDIM; ++d) { std::vector<double> tmp(grid.ncoarse); for (int i = 1; i <= grid.ncoarse; ++i) tmp[i-1] = grid.f(i, d); write_rec(tmp.data(), grid.ncoarse * 8); }
    for (int il = 1; il <= grid.nlevelmax; ++il) {
        for (int icpu = 1; icpu <= grid.ncpu + grid.nboundary; ++icpu) {
            int ncache = (icpu <= grid.ncpu) ? grid.numbl(icpu, il) : 0;
            if (ncache > 0) {
                int32_t il_rec = il, nc_rec = ncache;
                write_rec(&il_rec, 4); write_rec(&nc_rec, 4);
                std::vector<int> g_list; int ig = (icpu <= grid.ncpu) ? grid.get_headl(icpu, il) : 0; while(ig > 0) { g_list.push_back(ig); ig = grid.next[ig-1]; }
                for (int ic = 0; ic < constants::twotondim; ++ic) {
                    std::vector<double> tmp_p(ncache); for(int i=0; i<ncache; ++i) tmp_p[i] = grid.phi[grid.ncoarse + ic*grid.ngridmax + g_list[i] - 1]; write_rec(tmp_p.data(), ncache * 8);
                    for(int d=1; d<=NDIM; ++d) { std::vector<double> tmp_f(ncache); for(int i=0; i<ncache; ++i) tmp_f[i] = grid.f(grid.ncoarse + ic*grid.ngridmax + g_list[i], d); write_rec(tmp_f.data(), ncache * 8); }
                }
            }
        }
    }
}

} // namespace ramses
