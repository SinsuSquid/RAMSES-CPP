#include "ramses/RamsesWriter.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <vector>

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

template void RamsesWriter::write_record<int>(const int*, size_t);
template void RamsesWriter::write_record<real_t>(const real_t*, size_t);
template void RamsesWriter::write_record<char>(const char*, size_t);

void RamsesWriter::write_amr(const AmrGrid& grid) {
    int32_t val;
    val = grid.ncpu; write_record(&val, 1);
    val = NDIM; write_record(&val, 1);
    int32_t nxnyz[3] = {params::nx, params::ny, params::nz}; write_record(nxnyz, 3);
    val = grid.nlevelmax; write_record(&val, 1);
    val = grid.ngridmax; write_record(&val, 1);
    val = 0; write_record(&val, 1); 
    val = 0; write_record(&val, 1); 
    double boxlen = params::boxlen; write_record(&boxlen, 1);
    int32_t nout_vars[3] = {1, 1, 1}; write_record(nout_vars, 3);
    double tout = 0.0; write_record(&tout, 1);
    double aout = 0.0; write_record(&aout, 1);
    double t = 0.0; write_record(&t, 1);
    std::vector<double> dtold(grid.nlevelmax, 0.0); write_record(dtold.data(), grid.nlevelmax);
    std::vector<double> dtnew(grid.nlevelmax, 0.0); write_record(dtnew.data(), grid.nlevelmax);
    int32_t steps[2] = {0, 0}; write_record(steps, 2);
    double einit[3] = {0.0, 0.0, 0.0}; write_record(einit, 3);
    double cosmo[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; write_record(cosmo, 7);
    double aexp_hexp[5] = {1.0, 0.0, 1.0, 0.0, 0.0}; write_record(aexp_hexp, 5);
    double msph = 0.0; write_record(&msph, 1);
    std::vector<int32_t> headl(grid.ncpu * grid.nlevelmax, 0);
    std::vector<int32_t> taill(grid.ncpu * grid.nlevelmax, 0);
    std::vector<int32_t> numbl(grid.ncpu * grid.nlevelmax, 0);
    for(int ilevel=0; ilevel<grid.nlevelmax; ++ilevel) {
        for(int icpu=0; icpu<grid.ncpu; ++icpu) {
            numbl[ilevel * grid.ncpu + icpu] = grid.numbl(icpu+1, ilevel+1);
            headl[ilevel * grid.ncpu + icpu] = grid.headl(icpu+1, ilevel+1);
            taill[ilevel * grid.ncpu + icpu] = grid.taill(icpu+1, ilevel+1);
        }
    }
    write_record(headl.data(), headl.size());
    write_record(taill.data(), taill.size());
    write_record(numbl.data(), numbl.size());
    std::vector<int32_t> numbtot(10 * grid.nlevelmax, 0); write_record(numbtot.data(), numbtot.size());
    int32_t free_mem[5] = {grid.headf, grid.tailf, grid.numbf, 0, 0}; write_record(free_mem, 5);
    char ordering[128] = "hilbert"; write_record(ordering, 128);
    int32_t key_size = 8; write_record(&key_size, 1);
    file_.close();
}

void RamsesWriter::write_hydro(const AmrGrid& grid, int nlevelmax) {
    int32_t val;
    val = grid.ncpu; write_record(&val, 1);
    val = grid.nvar; write_record(&val, 1);
    val = NDIM; write_record(&val, 1);
    val = nlevelmax; write_record(&val, 1);
    val = 0; write_record(&val, 1); 
    file_.close();
}

} // namespace ramses
