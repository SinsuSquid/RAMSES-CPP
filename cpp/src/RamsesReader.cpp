#include "ramses/RamsesReader.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace ramses {

bool RamsesReader::load_amr(AmrGrid& grid) {
    if (!file_.is_open()) return false;

    // 1. Grid Variables
    uint32_t t1 = read_tag();
    int ncpu, ndim, nx, ny, nz, nlevelmax, ngridmax, nboundary, ngrid_current;
    real_t boxlen;

    file_.read(reinterpret_cast<char*>(&ncpu), 4);
    file_.read(reinterpret_cast<char*>(&ndim), 4);
    file_.read(reinterpret_cast<char*>(&nx), 4);
    file_.read(reinterpret_cast<char*>(&ny), 4);
    file_.read(reinterpret_cast<char*>(&nz), 4);
    file_.read(reinterpret_cast<char*>(&nlevelmax), 4);
    file_.read(reinterpret_cast<char*>(&ngridmax), 4);
    file_.read(reinterpret_cast<char*>(&nboundary), 4);
    file_.read(reinterpret_cast<char*>(&ngrid_current), 4);
    file_.read(reinterpret_cast<char*>(&boxlen), sizeof(real_t));
    uint32_t t2 = read_tag();

    std::cout << "Loading AMR Snapshot: ncpu=" << ncpu << " ngrid_current=" << ngrid_current << std::endl;

    // Allocate grid (nvar=5 default)
    grid.allocate(nx, ny, nz, ngridmax, 5, ncpu, nlevelmax);

    // 2. Time Variables
    for(int i=0; i<11; ++i) skip_record(); 

    // 3. Level/List Variables
    std::vector<int> buf;
    read_record(buf); // headl
    for(int i=0; i<ncpu*nlevelmax; ++i) grid.headl.data()[i] = buf[i];
    
    read_record(buf); // taill
    for(int i=0; i<ncpu*nlevelmax; ++i) grid.taill.data()[i] = buf[i];

    read_record(buf); // numbl
    for(int i=0; i<ncpu*nlevelmax; ++i) grid.numbl.data()[i] = buf[i];

    skip_record(); // numbtot
    
    // 4. Free Memory / Ordering
    skip_record(); // headf, tailf, ...
    skip_record(); // ordering
    skip_record(); // bound_key or bisec info

    // 5. Coarse Level
    read_record(buf); // son(1:ncoarse)
    for(int i=0; i<grid.ncoarse; ++i) grid.son[i+1] = buf[i];

    read_record(buf); // flag1(1:ncoarse)
    for(int i=0; i<grid.ncoarse; ++i) grid.flag1[i+1] = buf[i];

    read_record(buf); // cpu_map(1:ncoarse)
    for(int i=0; i<grid.ncoarse; ++i) grid.cpu_map[i+1] = buf[i];

    // 6. Fine Levels
    for (int ilevel = 1; ilevel <= nlevelmax; ++ilevel) {
        for (int ibound = 1; ibound <= nboundary + ncpu; ++ibound) {
            int ncache;
            if (ibound <= ncpu) {
                ncache = grid.numbl(ibound, ilevel);
            } else {
                ncache = 0; // Boundary reading skipped
            }

            if (ncache > 0) {
                std::vector<int> ind_grid;
                read_record(ind_grid);
                
                std::vector<int> iig;
                read_record(iig); // next
                for(int i=0; i<ncache; ++i) grid.next[ind_grid[i]-1] = iig[i];

                read_record(iig); // prev
                for(int i=0; i<ncache; ++i) grid.prev[ind_grid[i]-1] = iig[i];

                for (int idim = 1; idim <= NDIM; ++idim) {
                    std::vector<real_t> xdp;
                    read_record(xdp); // xg
                    for(int i=0; i<ncache; ++i) grid.get_xg(ind_grid[i], idim) = xdp[i];
                }

                read_record(iig); // father
                for(int i=0; i<ncache; ++i) grid.father[ind_grid[i]-1] = iig[i];

                for (int n = 1; n <= constants::twondim; ++n) {
                    read_record(iig); // nbor
                    for(int i=0; i<ncache; ++i) grid.nbor[(n-1)*grid.ngridmax + (ind_grid[i]-1)] = iig[i];
                }

                for (int n = 1; n <= constants::twotondim; ++n) {
                    read_record(iig); // son
                    int iskip = grid.ncoarse + (n-1)*grid.ngridmax;
                    for(int i=0; i<ncache; ++i) grid.son[ind_grid[i] + iskip] = iig[i];
                }

                for (int n = 1; n <= constants::twotondim; ++n) {
                    read_record(iig); // cpu_map
                    int iskip = grid.ncoarse + (n-1)*grid.ngridmax;
                    for(int i=0; i<ncache; ++i) grid.cpu_map[ind_grid[i] + iskip] = iig[i];
                }

                for (int n = 1; n <= constants::twotondim; ++n) {
                    read_record(iig); // flag1
                    int iskip = grid.ncoarse + (n-1)*grid.ngridmax;
                    for(int i=0; i<ncache; ++i) grid.flag1[ind_grid[i] + iskip] = iig[i];
                }
            }
        }
    }
    return true;
}

bool RamsesReader::load_hydro(AmrGrid& grid) {
    if (!file_.is_open()) return false;

    // 1. Header
    uint32_t t1 = read_tag();
    int ncpu, nvar, ndim, nlevelmax, nboundary;
    file_.read(reinterpret_cast<char*>(&ncpu), 4);
    file_.read(reinterpret_cast<char*>(&nvar), 4);
    file_.read(reinterpret_cast<char*>(&ndim), 4);
    file_.read(reinterpret_cast<char*>(&nlevelmax), 4);
    file_.read(reinterpret_cast<char*>(&nboundary), 4);
    uint32_t t2 = read_tag();

    real_t gamma = read_single<real_t>();
    std::cout << "Loading Hydro Snapshot: nvar=" << nvar << " gamma=" << gamma << std::endl;

    const real_t smallr = 1e-10;

    for (int ilevel = 1; ilevel <= nlevelmax; ++ilevel) {
        for (int ibound = 1; ibound <= nboundary + ncpu; ++ibound) {
            int lv = read_single<int>();
            int ncache = read_single<int>();

            if (ncache > 0) {
                std::vector<int> ind_grid;
                read_record(ind_grid);

                for (int ind = 1; ind <= constants::twotondim; ++ind) {
                    int iskip = grid.ncoarse + (ind - 1) * grid.ngridmax;
                    
                    std::vector<real_t> dens;
                    read_record(dens);
                    for(int i=0; i<ncache; ++i) grid.uold(ind_grid[i] + iskip, 1) = dens[i];

                    for (int idim = 1; idim <= NDIM; ++idim) {
                        std::vector<real_t> vel;
                        read_record(vel);
                        for(int i=0; i<ncache; ++i) {
                            grid.uold(ind_grid[i] + iskip, idim + 1) = vel[i] * grid.uold(ind_grid[i] + iskip, 1);
                        }
                    }

                    std::vector<real_t> pres;
                    read_record(pres);
                    for(int i=0; i<ncache; ++i) {
                        real_t d = grid.uold(ind_grid[i] + iskip, 1);
                        real_t e_kin = 0;
                        for(int idim=1; idim<=NDIM; ++idim) {
                            real_t mom = grid.uold(ind_grid[i] + iskip, idim + 1);
                            e_kin += 0.5 * mom * mom / std::max(d, smallr);
                        }
                        real_t e_int = pres[i] / (gamma - 1.0);
                        grid.uold(ind_grid[i] + iskip, NDIM + 2) = e_kin + e_int;
                    }

                    for (int iv = NDIM + 3; iv <= nvar; ++iv) {
                        std::vector<real_t> scalar;
                        read_record(scalar);
                        for(int i=0; i<ncache; ++i) {
                            grid.uold(ind_grid[i] + iskip, iv) = scalar[i] * grid.uold(ind_grid[i] + iskip, 1);
                        }
                    }
                }
            }
        }
    }
    return true;
}

} // namespace ramses
