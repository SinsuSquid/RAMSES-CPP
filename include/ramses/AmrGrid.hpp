#ifndef RAMSES_AMR_GRID_HPP
#define RAMSES_AMR_GRID_HPP

#include "Types.hpp"
#include <vector>
#include <string>
#include <functional>

namespace ramses {

class AmrGrid {
public:
    AmrGrid() : ncpu(1), ngridmax(1000), nlevelmax(10), nvar(5), ncoarse(1), ncell(2001), nboundary(0) {}

    int ncpu;
    int ngridmax;
    int nlevelmax;
    int nvar;
    int ncoarse;
    int ncell;
    int nboundary;
    real_t boxlen = 1.0;
    real_t gamma = 1.4;

    std::vector<int> headl_vec, taill_vec, numbl_vec;
    std::vector<int> next, prev, father, son, nbor, flag1, flag2, cpu_map;
    std::vector<real_t> xg, uold_vec, unew_vec;
    int headf, tailf, numbf;

    std::vector<real_t> rho, phi, f_vec;
    std::vector<qdp_t> hilbert_keys;
    real_t f(int icell, int idim) const { return f_vec.at((idim-1)*ncell + (icell-1)); }
    real_t& f(int icell, int idim) { return f_vec.at((idim-1)*ncell + (icell-1)); }

    // Coordinate utilities
    void get_cell_center(int icell, real_t xc[3]) const;
    int find_cell_by_coords(const real_t x[3], int level) const;

    // Particles
    int npart = 0, npartmax = 0;
    std::vector<real_t> xp, vp, mp, tp, zp;
    std::vector<uint8_t> family, tag;
    std::vector<int> idp, levelp, headp, nextp, prevp, tailp, numbp;
    int headp_free = 0, tailp_free = 0, numbp_free = 0;

    // Field indices
    int imetal = 0;
    int idelay = 0;

    std::vector<int> ibound_min, ibound_max, jbound_min, jbound_max, kbound_min, kbound_max, bound_type;

    void set_nboundary(int nb) { nboundary = nb; }
    void allocate(int nx, int ny, int nz, int ngridmax, int nvar, int ncpu, int nlevelmax);
    void resize_particles(int new_npartmax);
    int get_free_grid();
    void free_grid(int igrid);
    int get_free_particle();
    void free_particle(int ip);

    void add_to_level_list(int igrid, int ilevel);

    int get_headl(int icpu, int ilevel) const { 
        if (ncpu == 0) return 0;
        if (ilevel < 0 || ilevel > nlevelmax || icpu < 1 || icpu > ncpu + nboundary) return 0;
        return headl_vec.at((size_t)ilevel * ncpu + (icpu - 1)); 
    }
    int get_cell_level(int icell) const {
        if (icell <= ncoarse) return 0;
        int ig = ((icell - ncoarse - 1) % ngridmax) + 1;
        int level = 1; int curr_ig = ig;
        while(curr_ig > 0) {
            int father_cell = father[curr_ig-1];
            if (father_cell <= ncoarse) break;
            curr_ig = ((father_cell - ncoarse - 1) % ngridmax) + 1;
            level++;
        }
        return level;
    }
    int& headl(int icpu, int ilevel) { 
        return headl_vec.at((size_t)ilevel * ncpu + (icpu - 1)); 
    }
    int taill(int icpu, int ilevel) const { 
        if (ncpu == 0) return 0;
        if (ilevel < 0 || ilevel > nlevelmax || icpu < 1 || icpu > ncpu + nboundary) return 0;
        return taill_vec.at((size_t)ilevel * ncpu + (icpu - 1)); 
    }
    int& taill(int icpu, int ilevel) { 
        return taill_vec.at((size_t)ilevel * ncpu + (icpu - 1)); 
    }
    int numbl(int icpu, int ilevel) const { 
        if (ncpu == 0) return 0;
        if (ilevel < 0 || ilevel > nlevelmax || icpu < 1 || icpu > ncpu + nboundary) return 0;
        return numbl_vec.at((size_t)ilevel * ncpu + (icpu - 1)); 
    }
    int& numbl(int icpu, int ilevel) { 
        return numbl_vec.at((size_t)ilevel * ncpu + (icpu - 1)); 
    }
    int headl(int icpu, int ilevel) const { return get_headl(icpu, ilevel); }
    int count_grids_at_level(int ilevel) const;
    void synchronize_level_counts();

    real_t& uold(int icell, int ivar) { return uold_vec.at((size_t)(ivar-1)*ncell + (icell-1)); }
    real_t& unew(int icell, int ivar) { return unew_vec.at((size_t)(ivar-1)*ncell + (icell-1)); }
    real_t uold(int icell, int ivar) const { return uold_vec.at((size_t)(ivar-1)*ncell + (icell-1)); }
    real_t unew(int icell, int ivar) const { return unew_vec.at((size_t)(ivar-1)*ncell + (icell-1)); }

    void get_nbor_grids(int igrid, int ign[7]) const;
    void get_nbor_cells(const int ign[7], int ic, int icn[6], int igrid) const;
    void get_nbor_cells_coarse(int icell, int icn[6]) const;
    void get_27_cell_neighbors(int icell, int nbors[27]) const;
    
    // Interpolation hook
    using InterpolHook = std::function<void(const real_t[7][64], real_t[8][64])>;
    void set_interpol_hook(InterpolHook hook) { interpol_hook_ = hook; }
    void interpol(const real_t u1[7][64], real_t u2[8][64]) const { 
        if(interpol_hook_) interpol_hook_(u1, u2); 
        else {
            for(int i=0; i<8; ++i) for(int iv=0; iv<nvar; ++iv) u2[i][iv] = u1[0][iv];
        }
    }

    int nx = 1, ny = 1, nz = 1;

private:
    InterpolHook interpol_hook_;
};

} // namespace ramses

#endif
