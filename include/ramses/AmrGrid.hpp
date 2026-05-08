#ifndef RAMSES_AMR_GRID_HPP
#define RAMSES_AMR_GRID_HPP

#include "Types.hpp"
#include <vector>
#include <string>
#include <functional>

namespace ramses {

class AmrGrid {
public:
    AmrGrid() : ncpu(1), ngridmax(1000), nlevelmax(10), nvar(5), ncoarse(1), ncell(2001) {}

    int ncpu;
    int ngridmax;
    int nlevelmax;
    int nvar;
    int ncoarse;
    int ncell;
    real_t boxlen = 1.0;
    real_t gamma = 1.4;

    std::vector<int> headl_vec, taill_vec, numbl_vec;
    std::vector<int> next, prev, father, son, nbor, flag1, cpu_map;
    std::vector<real_t> xg, uold_vec, unew_vec;
    int headf, tailf, numbf;

    std::vector<real_t> rho, phi, f_vec;
    real_t f(int icell, int idim) const { return f_vec[(idim-1)*ncell + (icell-1)]; }
    real_t& f(int icell, int idim) { return f_vec[(idim-1)*ncell + (icell-1)]; }

    // Particles
    int npart = 0, npartmax = 0;
    std::vector<real_t> xp, vp, mp;
    std::vector<int> idp, levelp, headp, nextp, prevp, tailp, numbp;
    int headp_free, tailp_free, numbp_free;

    std::vector<int> ibound_min, ibound_max, jbound_min, jbound_max, kbound_min, kbound_max, bound_type;
    int nboundary;

    void allocate(int nx, int ny, int nz, int ngridmax, int nvar, int ncpu, int nlevelmax);
    void resize_particles(int new_npartmax);
    int get_free_grid();
    void free_grid(int igrid);
    void add_to_level_list(int igrid, int ilevel);
    int get_headl(int icpu, int ilevel) const { return headl_vec[(ilevel-1)*ncpu + (icpu-1)]; }
    int taill(int icpu, int ilevel) const { return taill_vec[(ilevel-1)*ncpu + (icpu-1)]; }
    int numbl(int icpu, int ilevel) const { return numbl_vec[(ilevel-1)*ncpu + (icpu-1)]; }
    int headl(int icpu, int ilevel) const { return get_headl(icpu, ilevel); }
    int count_grids_at_level(int ilevel) const;

    real_t& uold(int icell, int ivar) { return uold_vec[(ivar-1)*ncell + (icell-1)]; }
    real_t& unew(int icell, int ivar) { return unew_vec[(ivar-1)*ncell + (icell-1)]; }
    real_t uold(int icell, int ivar) const { return uold_vec[(ivar-1)*ncell + (icell-1)]; }
    real_t unew(int icell, int ivar) const { return unew_vec[(ivar-1)*ncell + (icell-1)]; }

    void get_nbor_grids(int igrid, int ign[7]) const;
    void get_nbor_cells(const int ign[7], int ic, int icn[6], int igrid) const;
    void get_27_cell_neighbors(int icell, int nbors[27]) const;
    
    // Interpolation hook
    using InterpolHook = std::function<void(const real_t[7][64], real_t[8][64])>;
    void set_interpol_hook(InterpolHook hook) { interpol_hook_ = hook; }
    void interpol(const real_t u1[7][64], real_t u2[8][64]) const { if(interpol_hook_) interpol_hook_(u1, u2); }

    int nx, ny, nz;

private:
    InterpolHook interpol_hook_;
};

} // namespace ramses

#endif
