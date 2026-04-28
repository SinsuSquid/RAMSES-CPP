#include "ramses/RamsesWriter.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <vector>
#include <cstring>

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
    // R1-R7
    val = grid.ncpu; write_record(&val, 1);
    val = NDIM; write_record(&val, 1);
    int32_t nxnyz[3] = {params::nx, params::ny, params::nz}; write_record(nxnyz, 3);
    val = grid.nlevelmax; write_record(&val, 1);
    val = grid.ngridmax; write_record(&val, 1);
    val = 0; write_record(&val, 1); // nboundary
    val = 0; write_record(&val, 1); // ngrid_current
    
    // R8
    double boxlen = params::boxlen; write_record(&boxlen, 1);
    
    // R9
    int32_t nout_iout_ifout[3] = {info.noutput, info.iout, 0}; write_record(nout_iout_ifout, 3);
    
    // R10
    std::vector<double> tout(info.noutput, 0.0);
    for(int i=0; i<std::min((int)info.tout.size(), info.noutput); ++i) tout[i] = info.tout[i];
    write_record(tout.data(), info.noutput);
    
    // R11
    std::vector<double> aout(info.noutput, 0.0); write_record(aout.data(), info.noutput);
    
    // R12
    double t = info.t; write_record(&t, 1);
    
    // R13, R14
    std::vector<double> dtold(grid.nlevelmax, 0.0); write_record(dtold.data(), grid.nlevelmax);
    std::vector<double> dtnew(grid.nlevelmax, 0.0); write_record(dtnew.data(), grid.nlevelmax);
    
    // R15
    int32_t nstep_coarse[2] = {info.nstep, info.nstep}; write_record(nstep_coarse, 2);
    
    // R16
    double einit[3] = {0.0, 0.0, 0.0}; write_record(einit, 3);
    
    // R17
    double cosmo[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; write_record(cosmo, 7);
    
    // R18
    double aexp_hexp[5] = {1.0, 0.0, 1.0, 0.0, 0.0}; write_record(aexp_hexp, 5);
    
    // R19
    double msph = 0.0; write_record(&msph, 1);
    
    // R20, R21, R22
    std::vector<int32_t> headl(grid.ncpu * grid.nlevelmax);
    std::vector<int32_t> taill(grid.ncpu * grid.nlevelmax);
    std::vector<int32_t> numbl(grid.ncpu * grid.nlevelmax);
    for(int ilevel=0; ilevel<grid.nlevelmax; ++ilevel) {
        for(int icpu=0; icpu<grid.ncpu; ++icpu) {
            headl[ilevel * grid.ncpu + icpu] = grid.headl(icpu+1, ilevel+1);
            taill[ilevel * grid.ncpu + icpu] = grid.taill(icpu+1, ilevel+1);
            numbl[ilevel * grid.ncpu + icpu] = grid.numbl(icpu+1, ilevel+1);
        }
    }
    write_record(headl.data(), headl.size());
    write_record(taill.data(), taill.size());
    write_record(numbl.data(), numbl.size());

    // R23
    std::vector<int32_t> numbtot(10 * grid.nlevelmax, 0); 
    for(int ilevel=0; ilevel<grid.nlevelmax; ++ilevel) {
        numbtot[ilevel] = grid.count_grids_at_level(ilevel+1);
    }
    write_record(numbtot.data(), numbtot.size());
    
    // R24-R26 skipped since nboundary=0
    
    // R27
    int32_t free_mem[5] = {grid.headf, grid.tailf, grid.numbf, 0, 0}; write_record(free_mem, 5);
    
    // R28
    char ordering[128] = {0}; std::strncpy(ordering, "hilbert", 127);
    write_record(ordering, 128);
    
    // R29
    int32_t key_size = 0; write_record(&key_size, 1);
    
    // R30, R31, R32 (Coarse level)
    std::vector<int32_t> son_coarse(grid.ncoarse);
    std::vector<int32_t> flag1_coarse(grid.ncoarse);
    std::vector<int32_t> cpu_map_coarse(grid.ncoarse);
    for(int i=0; i<grid.ncoarse; ++i) {
        son_coarse[i] = grid.son[i+1];
        flag1_coarse[i] = grid.flag1[i+1];
        cpu_map_coarse[i] = grid.cpu_map[i+1];
    }
    write_record(son_coarse.data(), grid.ncoarse);
    write_record(flag1_coarse.data(), grid.ncoarse);
    write_record(cpu_map_coarse.data(), grid.ncoarse);

    // Fine levels
    for(int ilevel=1; ilevel <= grid.nlevelmax; ++ilevel) {
        for(int icpu=1; icpu <= grid.ncpu; ++icpu) {
            int ncache = grid.numbl(icpu, ilevel);
            if (ncache > 0) {
                std::vector<int32_t> ind_grid(ncache);
                int curr = grid.headl(icpu, ilevel);
                for(int i=0; i<ncache; ++i) {
                    if (curr <= 0 || curr > grid.ngridmax) {
                        break;
                    }
                    ind_grid[i] = curr;
                    curr = grid.next[curr-1];
                }
                write_record(ind_grid.data(), ncache); // ind_grid
                
                std::vector<int32_t> next_grid(ncache);
                curr = grid.headl(icpu, ilevel);
                for(int i=0; i<ncache; ++i) {
                    next_grid[i] = grid.next[curr-1];
                    curr = grid.next[curr-1];
                }
                write_record(next_grid.data(), ncache); // next
                
                std::vector<int32_t> prev_grid(ncache);
                curr = grid.headl(icpu, ilevel);
                for(int i=0; i<ncache; ++i) {
                    prev_grid[i] = grid.prev[curr-1];
                    curr = grid.next[curr-1];
                }
                write_record(prev_grid.data(), ncache); // prev
                
                for(int idim=1; idim<=NDIM; ++idim) {
                    std::vector<double> xg_dim(ncache);
                    curr = grid.headl(icpu, ilevel);
                    for(int i=0; i<ncache; ++i) {
                        xg_dim[i] = grid.xg[(idim - 1) * grid.ngridmax + (curr - 1)];
                        curr = grid.next[curr-1];
                    }
                    write_record(xg_dim.data(), ncache);
                }
                
                std::vector<int32_t> father(ncache);
                curr = grid.headl(icpu, ilevel);
                for(int i=0; i<ncache; ++i) {
                    father[i] = grid.father[curr-1];
                    curr = grid.next[curr-1];
                }
                write_record(father.data(), ncache);
                
                for(int inbor=1; inbor<=2*NDIM; ++inbor) {
                    std::vector<int32_t> nbor_dim(ncache);
                    curr = grid.headl(icpu, ilevel);
                    for(int i=0; i<ncache; ++i) {
                        nbor_dim[i] = grid.nbor[(inbor - 1) * grid.ngridmax + (curr - 1)];
                        curr = grid.next[curr-1];
                    }
                    write_record(nbor_dim.data(), ncache);
                }
                
                for(int ison=1; ison<=constants::twotondim; ++ison) {
                    std::vector<int32_t> son_dim(ncache);
                    curr = grid.headl(icpu, ilevel);
                    for(int i=0; i<ncache; ++i) {
                        son_dim[i] = grid.son[grid.ncoarse + (ison-1)*grid.ngridmax + curr];
                        curr = grid.next[curr-1];
                    }
                    write_record(son_dim.data(), ncache);
                }
                
                for(int ison=1; ison<=constants::twotondim; ++ison) {
                    std::vector<int32_t> flag1_dim(ncache);
                    curr = grid.headl(icpu, ilevel);
                    for(int i=0; i<ncache; ++i) {
                        flag1_dim[i] = grid.flag1[grid.ncoarse + (ison-1)*grid.ngridmax + curr];
                        curr = grid.next[curr-1];
                    }
                    write_record(flag1_dim.data(), ncache);
                }

                for(int ison=1; ison<=constants::twotondim; ++ison) {
                    std::vector<int32_t> cpu_map_dim(ncache);
                    curr = grid.headl(icpu, ilevel);
                    for(int i=0; i<ncache; ++i) {
                        cpu_map_dim[i] = grid.cpu_map[grid.ncoarse + (ison-1)*grid.ngridmax + curr];
                        curr = grid.next[curr-1];
                    }
                    write_record(cpu_map_dim.data(), ncache);
                }
            }
        }
    }
    file_.close();
}

void RamsesWriter::write_hydro(const AmrGrid& grid, const SnapshotInfo& info) {
    int32_t val;
    val = grid.ncpu; write_record(&val, 1);
    val = grid.nvar; write_record(&val, 1);
    val = NDIM; write_record(&val, 1);
    val = grid.nlevelmax; write_record(&val, 1);
    val = 0; write_record(&val, 1); // nboundary
    double gamma = info.gamma; write_record(&gamma, 1);

    real_t smallr = 1e-10;

    for(int ilevel=1; ilevel <= grid.nlevelmax; ++ilevel) {
        for(int icpu=1; icpu <= grid.ncpu; ++icpu) {
            int ncache = grid.numbl(icpu, ilevel);
            val = ilevel; write_record(&val, 1);
            val = ncache; write_record(&val, 1);
            if (ncache > 0) {
                for(int ison=1; ison<=constants::twotondim; ++ison) {
                    // We need to write variables in the correct order:
                    // 1. density
                    // 2-4. velocity
                    // 5-7. B_left (if MHD)
                    // 8-10. B_right (if MHD)
                    // 11. pressure
                    // ... others ...
                    
                    int nvar_out = grid.nvar;
                    for (int ivar_out = 1; ivar_out <= nvar_out; ++ivar_out) {
                        std::vector<double> data(ncache);
                        int curr = grid.headl(icpu, ilevel);
                        for(int i=0; i<ncache; ++i) {
                            int ind_cell = grid.ncoarse + (ison-1)*grid.ngridmax + curr;
                            
                            // Map ivar_out to physical variable and transform
#ifndef MHD
                            if (ivar_out == 1) data[i] = grid.unew(ind_cell, 1); // density
                            else if (ivar_out <= 1 + NDIM) {
                                real_t d = std::max(grid.unew(ind_cell, 1), smallr);
                                data[i] = grid.unew(ind_cell, ivar_out) / d; // velocity
                            } else if (ivar_out == 2 + NDIM) { // pressure
                                real_t d = std::max(grid.unew(ind_cell, 1), smallr);
                                real_t ekin = 0;
                                for(int id=1; id<=NDIM; ++id) ekin += 0.5 * grid.unew(ind_cell, 1+id)*grid.unew(ind_cell, 1+id)/d;
                                data[i] = (gamma - 1.0) * (grid.unew(ind_cell, 2+NDIM) - ekin);
                            } else {
                                data[i] = grid.unew(ind_cell, ivar_out);
                            }
#else
                            // MHD Mapping (grid.nvar = nhydro + nener + 3)
                            // nhydro = 8: rho, mx, my, mz, E, BxL, ByL, BzL
                            int nener = grid.nvar - 11;
                            if (ivar_out == 1) data[i] = grid.unew(ind_cell, 1); // density
                            else if (ivar_out <= 4) { // velocity
                                real_t d = std::max(grid.unew(ind_cell, 1), smallr);
                                data[i] = grid.unew(ind_cell, ivar_out) / d;
                            } else if (ivar_out >= 5 && ivar_out <= 7) { // B_left
                                data[i] = grid.unew(ind_cell, ivar_out + 1); // 6,7,8
                            } else if (ivar_out >= 8 && ivar_out <= 10) { // B_right
                                data[i] = grid.unew(ind_cell, grid.nvar - 3 + (ivar_out - 7));
                            } else if (ivar_out == 11) { // pressure
                                real_t d = std::max(grid.unew(ind_cell, 1), smallr);
                                real_t ekin = 0.5*(grid.unew(ind_cell, 2)*grid.unew(ind_cell, 2) + 
                                                   grid.unew(ind_cell, 3)*grid.unew(ind_cell, 3) + 
                                                   grid.unew(ind_cell, 4)*grid.unew(ind_cell, 4))/d;
                                real_t Bx = 0.5*(grid.unew(ind_cell, 6) + grid.unew(ind_cell, grid.nvar-2));
                                real_t By = 0.5*(grid.unew(ind_cell, 7) + grid.unew(ind_cell, grid.nvar-1));
                                real_t Bz = 0.5*(grid.unew(ind_cell, 8) + grid.unew(ind_cell, grid.nvar));
                                real_t emag = 0.5*(Bx*Bx + By*By + Bz*Bz);
                                data[i] = (gamma - 1.0) * (grid.unew(ind_cell, 5) - ekin - emag);
                            } else { // nener or scalars
                                data[i] = grid.unew(ind_cell, ivar_out > 11 ? ivar_out - 3 : ivar_out); // this needs careful mapping if nener > 0
                                // For now, assume ivar_out > 11 are nener/scalars which were at 9...nvar
                                if (ivar_out > 11) data[i] = grid.unew(ind_cell, 9 + (ivar_out - 12));
                            }
#endif
                            curr = grid.next[curr-1];
                        }
                        write_record(data.data(), ncache);
                    }
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
    val = 0; write_record(&val, 1); // nboundary

    for(int ilevel=1; ilevel <= grid.nlevelmax; ++ilevel) {
        for(int icpu=1; icpu <= grid.ncpu; ++icpu) {
            int ncache = grid.numbl(icpu, ilevel);
            val = ilevel; write_record(&val, 1);
            val = ncache; write_record(&val, 1);
            if (ncache > 0) {
                for(int ison=1; ison<=constants::twotondim; ++ison) {
                    std::vector<double> phi_data(ncache);
                    int curr = grid.headl(icpu, ilevel);
                    for(int i=0; i<ncache; ++i) {
                        int ind_cell = grid.ncoarse + (ison-1)*grid.ngridmax + curr;
                        phi_data[i] = grid.phi[ind_cell];
                        curr = grid.next[curr-1];
                    }
                    write_record(phi_data.data(), ncache);

                    for(int idim=1; idim<=NDIM; ++idim) {
                        std::vector<double> f_data(ncache);
                        curr = grid.headl(icpu, ilevel);
                        for(int i=0; i<ncache; ++i) {
                            int ind_cell = grid.ncoarse + (ison-1)*grid.ngridmax + curr;
                            f_data[i] = grid.f(ind_cell, idim);
                            curr = grid.next[curr-1];
                        }
                        write_record(f_data.data(), ncache);
                    }
                }
            }
        }
    }
    file_.close();
}

} // namespace ramses
