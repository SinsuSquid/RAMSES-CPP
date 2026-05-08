#include "ramses/RamsesReader.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <vector>
#include <cstring>
#include <algorithm>

namespace ramses {
RamsesReader::RamsesReader(const std::string& filename) : filename_(filename) {}

bool RamsesReader::load_amr(AmrGrid& grid) {
    std::ifstream file(filename_, std::ios::binary);
    if (!file.is_open()) return false;
    int32_t ncpu, ndim, nx, ny, nz, nlevelmax, ngridmax;
    read_record(file, ncpu);
    read_record(file, ndim);
    std::vector<int> nxnyz; read_record(file, nxnyz); 
    if(nxnyz.size()>=3) { nx=nxnyz[0]; ny=nxnyz[1]; nz=nxnyz[2]; } else { nx=ny=nz=1; }
    read_record(file, nlevelmax);
    read_record(file, ngridmax);
    int32_t nboundary, ngrid_current; double boxlen;
    read_record(file, nboundary);
    read_record(file, ngrid_current);
    read_record(file, boxlen);

    grid.allocate(nx, ny, nz, ngridmax, grid.nvar, ncpu, nlevelmax);
    grid.nboundary = nboundary;

    int32_t noutput, iout, ifout;
    std::vector<int> rec9; read_record(file, rec9); 
    if(rec9.size()>=3) { noutput=rec9[0]; iout=rec9[1]; ifout=rec9[2]; } else { noutput=1; }
    
    std::vector<double> tout; read_record(file, tout);
    std::vector<double> aout; read_record(file, aout);
    double t; read_record(file, t);
    std::vector<double> dtold; read_record(file, dtold);
    std::vector<double> dtnew; read_record(file, dtnew);
    std::vector<int> rec15; read_record(file, rec15);
    
    // EINIT, MASS_TOT_0, RHO_TOT
    std::vector<double> rec16; read_record(file, rec16);
    // COSMO
    std::vector<double> rec17; read_record(file, rec17);
    // AEXP, HEXP...
    std::vector<double> rec18; read_record(file, rec18);
    // MASS_SPH
    double msph; read_record(file, msph);
    
    std::vector<int> hl_buf, tl_buf, nl_buf;
    read_record(file, hl_buf); 
    read_record(file, tl_buf);
    read_record(file, nl_buf);
    for(int i=0; i<std::min((int)hl_buf.size(), (int)grid.headl_vec.size()); ++i) grid.headl_vec[i] = hl_buf[i];
    for(int i=0; i<std::min((int)tl_buf.size(), (int)grid.taill_vec.size()); ++i) grid.taill_vec[i] = tl_buf[i];
    for(int i=0; i<std::min((int)nl_buf.size(), (int)grid.numbl_vec.size()); ++i) grid.numbl_vec[i] = nl_buf[i];
    
    std::vector<int> nt; read_record(file, nt); 
    if (nboundary > 0) { std::vector<int> hb, tb, nb; read_record(file, hb); read_record(file, tb); read_record(file, nb); }
    
    std::vector<int> mem; read_record(file, mem);
    char ord[128]; read_record_raw(file, ord, 128);
    std::vector<double> bk; read_record(file, bk);
    
    read_record_raw(file, grid.son.data(), grid.ncoarse * sizeof(int));
    read_record_raw(file, grid.flag1.data(), grid.ncoarse * sizeof(int));
    read_record_raw(file, grid.cpu_map.data(), grid.ncoarse * sizeof(int));
    
    for(int il=1; il<=nlevelmax; ++il) for(int ic=1; ic<=ncpu + nboundary; ++ic) {
        int nca = (ic <= ncpu) ? grid.numbl(ic, il) : 0;
        if (nca > 0) {
            std::vector<int> igs; read_record(file, igs);
            std::vector<int> nxt; read_record(file, nxt);
            std::vector<int> prv; read_record(file, prv);
            for(int i=0; i<nca; ++i) { grid.next[igs[i]-1] = nxt[i]; grid.prev[igs[i]-1] = prv[i]; }
            for(int d=1; d<=ndim; ++d) { std::vector<double> xgd; read_record(file, xgd); for(int i=0; i<nca; ++i) grid.xg[(d-1)*ngridmax + (igs[i]-1)] = xgd[i]; }
            std::vector<int> ftr; read_record(file, ftr); for(int i=0; i<nca; ++i) grid.father[igs[i]-1] = ftr[i];
            for(int n=1; n<=2*ndim; ++n) { std::vector<int> nbd; read_record(file, nbd); }
            for(int s=1; s<=constants::twotondim; ++s) { std::vector<int> sd; read_record(file, sd); for(int i=0; i<nca; ++i) grid.son[grid.ncoarse+(s-1)*ngridmax + igs[i]-1] = sd[i]; }
            for(int s=1; s<=constants::twotondim; ++s) { std::vector<int> cmd; read_record(file, cmd); for(int i=0; i<nca; ++i) grid.cpu_map[grid.ncoarse+(s-1)*ngridmax + igs[i]-1] = cmd[i]; }
            for(int s=1; s<=constants::twotondim; ++s) { std::vector<int> f1d; read_record(file, f1d); for(int i=0; i<nca; ++i) grid.flag1[grid.ncoarse+(s-1)*ngridmax + igs[i]-1] = f1d[i]; }
        }
    }
    return true;
}

bool RamsesReader::load_hydro(AmrGrid& grid) {
    std::ifstream file(filename_, std::ios::binary);
    if (!file.is_open()) return false;
    int32_t ncpu, nvar, ndim, nlevelmax, nboundary; double gamma;
    read_record(file, ncpu);
    read_record(file, nvar);
    read_record(file, ndim);
    read_record(file, nlevelmax);
    read_record(file, nboundary);
    read_record(file, gamma);

    // Coarse level hydro
    for (int iv = 1; iv <= nvar; ++iv) {
        std::vector<double> buf; read_record(file, buf);
        if (iv <= grid.nvar) for (int i = 1; i <= grid.ncoarse; ++i) grid.uold(i, iv) = buf[i-1];
    }

    for(int il=1; il<=nlevelmax; ++il) for(int ic=1; ic<=ncpu + nboundary; ++ic) {
        int ilevel2, nca; read_record(file, ilevel2); read_record(file, nca);
        if (nca > 0) {
            std::vector<int> igs(nca); int head = grid.headl(ic, il); 
            for(int i=0; i<nca; ++i) { 
                if (head <= 0) break;
                igs[i] = head; head = grid.next[head-1]; 
            }
            for(int ind=1; ind<=constants::twotondim; ++ind) {
                int iskip = grid.ncoarse + (ind-1)*grid.ngridmax;
                for(int iv=1; iv<=nvar; ++iv) { 
                    std::vector<double> buf; read_record(file, buf); 
                    if(iv<=grid.nvar) {
                        for(int i=0; i<(int)std::min((int)nca, (int)buf.size()); ++i) {
                            if (igs[i] > 0) grid.uold(igs[i]+iskip, iv) = buf[i];
                        }
                    }
                }
            }
        }
    }
    return true;
}

template <typename T> void RamsesReader::read_record(std::ifstream& f, T& d) { int32_t s; f.read((char*)&s, 4); f.read((char*)&d, sizeof(T)); f.read((char*)&s, 4); }
template <typename T> void RamsesReader::read_record(std::ifstream& f, std::vector<T>& d) { int32_t s; f.read((char*)&s, 4); d.resize(s/sizeof(T)); if(s>0) f.read((char*)d.data(), s); f.read((char*)&s, 4); }
bool RamsesReader::load_particles(AmrGrid& grid) {
    std::ifstream file(filename_, std::ios::binary);
    if (!file.is_open()) return false;
    int32_t ncpu, ndim, npart;
    read_record(file, ncpu);
    read_record(file, ndim);
    read_record(file, npart);
    
    if (npart > 0) {
        grid.resize_particles(npart);
        grid.npart = npart;
        
        // Skip seeds and other header info
        int32_t s; file.read((char*)&s, 4); file.seekg(s, std::ios::cur); file.read((char*)&s, 4); // localseed
        file.read((char*)&s, 4); file.seekg(s, std::ios::cur); file.read((char*)&s, 4); // nstar etc.
        
        // Positions
        for (int d = 0; d < ndim; ++d) {
            std::vector<double> buf; read_record(file, buf);
            for (int i = 0; i < npart; ++i) grid.xp[d * grid.npartmax + i] = buf[i];
        }
        // Velocities
        for (int d = 0; d < ndim; ++d) {
            std::vector<double> buf; read_record(file, buf);
            for (int i = 0; i < npart; ++i) grid.vp[d * grid.npartmax + i] = buf[i];
        }
        // Mass
        {
            std::vector<double> buf; read_record(file, buf);
            for (int i = 0; i < npart; ++i) grid.mp[i] = buf[i];
        }
        // ID
        {
            std::vector<int32_t> buf; read_record(file, buf);
            for (int i = 0; i < npart; ++i) grid.idp[i] = buf[i];
        }
        // Level
        {
            std::vector<int32_t> buf; read_record(file, buf);
            for (int i = 0; i < npart; ++i) grid.levelp[i] = buf[i];
        }
    }
    return true;
}

void RamsesReader::read_record_raw(std::ifstream& f, void* d, size_t s_max) { int32_t s; f.read((char*)&s, 4); size_t r = std::min((size_t)s, s_max); if(r>0) f.read((char*)d, r); if((size_t)s > s_max) f.seekg(s - s_max, std::ios::cur); f.read((char*)&s, 4); }
void RamsesReader::read_record(std::ifstream& f, void* d, size_t s) { read_record_raw(f, d, s); }
}
