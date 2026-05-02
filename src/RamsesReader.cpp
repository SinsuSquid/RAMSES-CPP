#include "ramses/RamsesReader.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <vector>
#include <cstring>

namespace ramses {
RamsesReader::RamsesReader(const std::string& filename) : filename_(filename) {}
bool RamsesReader::load_amr(AmrGrid& grid) {
    std::ifstream file(filename_, std::ios::binary);
    if (!file.is_open()) return false;
    int32_t ncpu, ndim, nx, ny, nz, nlevelmax, ngridmax;
    read_record(file, ncpu); read_record(file, ndim);
    std::vector<int> nxnyz; read_record(file, nxnyz); nx=nxnyz[0]; ny=nxnyz[1]; nz=nxnyz[2];
    read_record(file, nlevelmax); read_record(file, ngridmax);
    int32_t nboundary, ngrid_current; double boxlen;
    read_record(file, nboundary); read_record(file, ngrid_current); read_record(file, boxlen);
    grid.allocate(nx, ny, nz, ngridmax, grid.nvar, ncpu, nlevelmax);
    grid.nboundary = nboundary;
    int32_t noutput, iout, ifout; read_record(file, noutput); // wait, grouped? no.
    // wait, I need to match the hexdump exactly.
    // Hex dump showed: ncpu(4), ndim(4), nxnyz(12), nlevelmax(4).
    // All separate except nxnyz.
    // My code above does that.
    
    // Resume after nlevelmax:
    // ngridmax(Rec 5), nboundary(Rec 6), ngrid_current(Rec 7), boxlen(Rec 8)
    // noutput, iout, ifout (Rec 9) - grouped in hexdump?
    // Let's check hexdump at 0x2c: 04 00 00 00 | 0a 00 00 00 | 04 00 00 00
    // Record 4 is nlevelmax=10.
    // Next: Record 5 starts at 0x38.
    // 0x38: 04 00 00 00 | e8 03 00 00 | 04 00 00 00
    // Record 5 is ngridmax=1000.
    
    // Okay, so they ARE separate.
    
    // Skip to Record 9:
    for(int i=0; i<3; ++i) { int dummy_i; read_record(file, dummy_i); } // Skip 6, 7, 8 (nboundary, ngrid_current, boxlen)
    std::vector<int> rec9; read_record(file, rec9); // Rec 9: noutput, iout, ifout
    std::vector<double> rec10, rec11; read_record(file, rec10); read_record(file, rec11); // Rec 10, 11: tout, aout
    double t; read_record(file, t); // Rec 12: t
    std::vector<double> rec13, rec14; read_record(file, rec13); read_record(file, rec14); // Rec 13, 14: dtold, dtnew
    std::vector<int> rec15; read_record(file, rec15); // Rec 15: nstep, nstep_coarse
    for(int i=0; i<4; ++i) { std::vector<double> dummy_f; read_record(file, dummy_f); } // Rec 16-19: einit, cosmo, aexp, msph
    
    std::vector<int> hl_buf, tl_buf, nl_buf;
    read_record(file, hl_buf); for(int i=0; i<ncpu*nlevelmax; ++i) grid.headl_vec[i] = hl_buf[i];
    read_record(file, tl_buf); for(int i=0; i<ncpu*nlevelmax; ++i) grid.taill_vec[i] = tl_buf[i];
    read_record(file, nl_buf); for(int i=0; i<ncpu*nlevelmax; ++i) grid.numbl_vec[i] = nl_buf[i];
    
    std::vector<int> nt; read_record(file, nt); // Rec 23: numbtot
    if (nboundary > 0) { std::vector<int> hb, tb, nb; read_record(file, hb); read_record(file, tb); read_record(file, nb); }
    
    std::vector<int> mem; read_record(file, mem); // Rec 27/24: headf...
    char ord[128]; read_record_raw(file, ord, 128); // Rec 28/25: ordering
    std::vector<double> bk; read_record(file, bk); // Rec 29/26: bound_key
    
    read_record_raw(file, grid.son.data(), grid.ncoarse * sizeof(int));
    read_record_raw(file, grid.flag1.data(), grid.ncoarse * sizeof(int));
    read_record_raw(file, grid.cpu_map.data(), grid.ncoarse * sizeof(int));
    
    for(int il=1; il<=nlevelmax; ++il) for(int ic=1; ic<=ncpu + nboundary; ++ic) {
        int nca = (ic <= ncpu) ? grid.numbl(ic, il) : 0;
        if (nca > 0) {
            std::vector<int> igs; read_record(file, igs);
            std::vector<int> nxt; read_record(file, nxt); for(int i=0; i<nca; ++i) grid.next[igs[i]-1] = nxt[i];
            std::vector<int> prv; read_record(file, prv); for(int i=0; i<nca; ++i) grid.prev[igs[i]-1] = prv[i];
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
    int32_t s; file.read((char*)&s, 4); file.read((char*)&ncpu, 4); file.read((char*)&nvar, 4); file.read((char*)&ndim, 4); file.read((char*)&nlevelmax, 4); file.read((char*)&nboundary, 4); file.read((char*)&gamma, 8); file.read((char*)&s, 4);
    for(int il=1; il<=nlevelmax; ++il) for(int ic=1; ic<=ncpu + nboundary; ++ic) {
        int ilevel2, nca; read_record(file, ilevel2); read_record(file, nca);
        if (nca > 0) {
            std::vector<int> igs(nca); int head = grid.headl(ic, il); for(int i=0; i<nca; ++i) { igs[i] = head; head = grid.next[head-1]; }
            for(int ind=1; ind<=constants::twotondim; ++ind) {
                int iskip = grid.ncoarse + (ind-1)*grid.ngridmax;
                for(int iv=1; iv<=nvar; ++iv) { std::vector<double> buf; read_record(file, buf); if(iv<=grid.nvar) for(int i=0; i<nca; ++i) grid.uold(igs[i]+iskip, iv) = buf[i]; }
            }
        }
    }
    return true;
}
template <typename T> void RamsesReader::read_record(std::ifstream& f, T& d) { int32_t s; f.read((char*)&s, 4); f.read((char*)&d, sizeof(T)); f.read((char*)&s, 4); }
template <typename T> void RamsesReader::read_record(std::ifstream& f, std::vector<T>& d) { int32_t s; f.read((char*)&s, 4); d.resize(s/sizeof(T)); f.read((char*)d.data(), s); f.read((char*)&s, 4); }
void RamsesReader::read_record_raw(std::ifstream& f, void* d, size_t s_max) { int32_t s; f.read((char*)&s, 4); size_t r = std::min((size_t)s, s_max); f.read((char*)d, r); if((size_t)s > s_max) f.seekg(s - s_max, std::ios::cur); f.read((char*)&s, 4); }
}
