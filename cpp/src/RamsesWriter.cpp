#include "ramses/RamsesWriter.hpp"
#include "ramses/Parameters.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace ramses {

void RamsesWriter::write_amr(const AmrGrid& grid) {
    if (!file_.is_open()) return;
    write_single<int>(grid.ncpu);
    uint32_t header_tag = 9 * 4 + sizeof(real_t);
    write_tag(header_tag);
    file_.write(reinterpret_cast<const char*>(&grid.ncpu), 4);
    int ndim = NDIM; file_.write(reinterpret_cast<const char*>(&ndim), 4);
    file_.write(reinterpret_cast<const char*>(&params::nx), 4);
    file_.write(reinterpret_cast<const char*>(&params::ny), 4);
    file_.write(reinterpret_cast<const char*>(&params::nz), 4);
    file_.write(reinterpret_cast<const char*>(&grid.nlevelmax), 4);
    file_.write(reinterpret_cast<const char*>(&grid.ngridmax), 4);
    int nboundary = 0; file_.write(reinterpret_cast<const char*>(&nboundary), 4);
    int ngrid_current = 0; file_.write(reinterpret_cast<const char*>(&ngrid_current), 4);
    real_t boxlen = 1.0; file_.write(reinterpret_cast<const char*>(&boxlen), sizeof(real_t));
    write_tag(header_tag);
    for(int i=0; i<11; ++i) write_single<int>(0);
    std::vector<int> buf(grid.ncpu * grid.nlevelmax);
    for(size_t i=0; i<buf.size(); ++i) buf[i] = grid.headl.data()[i];
    write_record(buf);
    for(size_t i=0; i<buf.size(); ++i) buf[i] = grid.taill.data()[i];
    write_record(buf);
    for(size_t i=0; i<buf.size(); ++i) buf[i] = grid.numbl.data()[i];
    write_record(buf);
    write_single<int>(0);
    write_single<int>(0); write_single<int>(0); write_single<int>(0);
    std::vector<int> c_buf(grid.ncoarse);
    for(int i=0; i<grid.ncoarse; ++i) c_buf[i] = grid.son[i+1];
    write_record(c_buf);
    for(int i=0; i<grid.ncoarse; ++i) c_buf[i] = grid.flag1[i+1];
    write_record(c_buf);
    for(int i=0; i<grid.ncoarse; ++i) c_buf[i] = grid.cpu_map[i+1];
    write_record(c_buf);
}

void RamsesWriter::write_hydro(const AmrGrid& grid, int ilevel_max) {
    if (!file_.is_open()) return;

    // 1. Header
    uint32_t header_tag = 5 * 4;
    write_tag(header_tag);
    file_.write(reinterpret_cast<const char*>(&grid.ncpu), 4);
    file_.write(reinterpret_cast<const char*>(&grid.nvar), 4);
    int ndim = NDIM; file_.write(reinterpret_cast<const char*>(&ndim), 4);
    file_.write(reinterpret_cast<const char*>(&grid.nlevelmax), 4);
    int nboundary = 0; file_.write(reinterpret_cast<const char*>(&nboundary), 4);
    write_tag(header_tag);

    write_single<real_t>(1.4); // gamma dummy

    for (int ilevel = 1; ilevel <= grid.nlevelmax; ++ilevel) {
        for (int icpu = 1; icpu <= grid.ncpu; ++icpu) {
            int ncache = grid.numbl(icpu, ilevel);
            write_single<int>(ilevel);
            write_single<int>(ncache);

            if (ncache > 0) {
                std::vector<int> ind_grid;
                int curr = grid.headl(icpu, ilevel);
                while(curr > 0) {
                    ind_grid.push_back(curr);
                    curr = grid.next[curr-1];
                }

                write_record(ind_grid);

                for (int ind = 1; ind <= constants::twotondim; ++ind) {
                    int iskip = grid.ncoarse + (ind - 1) * grid.ngridmax;
                    
                    // RAMSES expects PRIMITIVE: density, velocity, pressure
                    std::vector<real_t> dens(ncache);
                    for(int i=0; i<ncache; ++i) dens[i] = grid.uold(ind_grid[i] + iskip, 1);
                    write_record(dens);

                    for (int idim = 1; idim <= NDIM; ++idim) {
                        std::vector<real_t> vel(ncache);
                        for(int i=0; i<ncache; ++i) vel[i] = grid.uold(ind_grid[i] + iskip, idim + 1) / std::max(dens[i], static_cast<real_t>(1e-10));
                        write_record(vel);
                    }

                    std::vector<real_t> pres(ncache);
                    for(int i=0; i<ncache; ++i) {
                        real_t d = dens[i];
                        real_t e_kin = 0;
                        for(int idim=1; idim<=NDIM; ++idim) {
                            real_t mom = grid.uold(ind_grid[i] + iskip, idim + 1);
                            e_kin += 0.5 * mom * mom / std::max(d, static_cast<real_t>(1e-10));
                        }
                        real_t e_int = grid.uold(ind_grid[i] + iskip, NDIM + 2) - e_kin;
                        pres[i] = e_int * (1.4 - 1.0); // gamma=1.4 dummy
                    }
                    write_record(pres);
                }
            }
        }
    }
}

} // namespace ramses
