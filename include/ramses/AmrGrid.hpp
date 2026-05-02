#ifndef RAMSES_AMR_GRID_HPP
#define RAMSES_AMR_GRID_HPP

#include "Types.hpp"
#include <vector>
#include <string>

namespace ramses {

class AmrGrid {
public:
    AmrGrid() : ncpu(1), ngridmax(1000), nlevelmax(10), nvar(5), boxlen(1.0) {}

    void allocate(int nx, int ny, int nz, int ngridmax, int nvar, int ncpu, int nlevelmax);

    int ncpu, ngridmax, nlevelmax, nvar, ncoarse, ncell;
    real_t boxlen, gamma;

    std::vector<int> headl_vec, taill_vec, numbl_vec; 
    std::vector<int> next, prev, father, son, flag1, cpu_map;
    std::vector<real_t> xg, uold_vec, unew_vec;
    
    std::vector<real_t> phi, rho;
    std::vector<real_t> f_vec; 

    int nboundary = 0;
    std::vector<int> ibound_min, ibound_max, bound_type;

    int get_headl(int icpu, int ilevel) const { return headl_vec[(ilevel-1)*ncpu + (icpu-1)]; }
    int headl(int icpu, int ilevel) const { return get_headl(icpu, ilevel); }
    int taill(int icpu, int ilevel) const { return taill_vec[(ilevel-1)*ncpu + (icpu-1)]; }
    int numbl(int icpu, int ilevel) const { return numbl_vec[(ilevel-1)*ncpu + (icpu-1)]; }

    int headf, tailf, numbf;
    int get_free_grid();
    void add_to_level_list(int igrid, int ilevel);
    int count_grids_at_level(int ilevel) const;

    real_t& uold(int icell, int ivar) { return uold_vec[(ivar-1)*ncell + (icell-1)]; }
    real_t& unew(int icell, int ivar) { return unew_vec[(ivar-1)*ncell + (icell-1)]; }
    const real_t& uold(int icell, int ivar) const { return uold_vec[(ivar-1)*ncell + (icell-1)]; }
    
    real_t& f(int icell, int idim) { return f_vec[(idim-1)*ncell + (icell-1)]; }
    const real_t& f(int icell, int idim) const { return f_vec[(idim-1)*ncell + (icell-1)]; }

    void get_nbor_grids(int igrid, int ign[7]) const;
    void get_nbor_cells(const int ign[7], int ic, int icn[6], int igrid) const;
};

} // namespace ramses

#endif // RAMSES_AMR_GRID_HPP
