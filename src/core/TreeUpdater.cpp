#include "ramses/core/TreeUpdater.hpp"
#include "ramses/core/AmrGrid.hpp"
#include "ramses/core/Constants.hpp"
#include "ramses/core/MpiManager.hpp"
#include "ramses/utils/Logger.hpp"
#include <algorithm>
#include <cmath>
#include <vector>
#include "ramses/core/Parameters.hpp"
#include <iostream>

namespace ramses {

TreeUpdater::TreeUpdater(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

int get_nbor_of_coarse(const AmrGrid& grid, int ic, int idim, int side) {
    int icn[6] = {0};
    grid.get_nbor_cells_coarse(ic, icn);
    return icn[idim * 2 + side];
}

void TreeUpdater::make_grid_fine(int ilevel) {
    if (ilevel >= grid_.nlevelmax) return;
    authorize_fine(ilevel);
    int myid = MpiManager::instance().rank() + 1, n2d = (1 << NDIM);
    std::vector<int> cells_to_refine;
    
    if (ilevel == 0) {
        for (int ic = 1; ic <= grid_.ncoarse; ++ic) {
            if (grid_.flag1[ic - 1] == 1 && grid_.flag2[ic - 1] == 1 && grid_.son[ic - 1] == 0) cells_to_refine.push_back(ic);
        }
    } else {
        int ig = grid_.get_headl(myid, ilevel);
        while (ig > 0) {
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                if (grid_.flag1[idc] == 1 && grid_.flag2[idc] == 1 && grid_.son[idc] == 0) cells_to_refine.push_back(idc + 1);
            }
            ig = grid_.next[ig - 1];
        }
    }

    int ncreate = cells_to_refine.size();
    if (ramses::params::verbose) {
        RAMSES_INFO("    Entering refine_fine for level  {}", ilevel);
        RAMSES_INFO("   ==> Make      {} sub-grids", ncreate);
    }

    for (size_t idx = 0; idx < cells_to_refine.size(); ++idx) {
        int ic_coarse = cells_to_refine[idx];
        int idc = ic_coarse - 1;
        int old_ngrid = grid_.ngridmax;
        int new_ig = grid_.get_free_grid(); if (new_ig == 0) return;
        
        if (grid_.ngridmax != old_ngrid) {
            // A grid resize happened! Shift all cell indices in cells_to_refine.
            for (size_t j = 0; j < cells_to_refine.size(); ++j) {
                if (cells_to_refine[j] > grid_.ncoarse) {
                    int cell_ig = (cells_to_refine[j] - 1 - grid_.ncoarse) % old_ngrid + 1;
                    int ic = (cells_to_refine[j] - 1 - grid_.ncoarse) / old_ngrid + 1;
                    cells_to_refine[j] = grid_.ncoarse + (ic - 1) * grid_.ngridmax + cell_ig;
                }
            }
            ic_coarse = cells_to_refine[idx];
            idc = ic_coarse - 1;
        }

        grid_.son[idc] = new_ig; grid_.father[new_ig - 1] = idc + 1;
        
        if (ilevel == 0) {
            for (int idim = 0; idim < 3; ++idim) {
                for (int side = 0; side < 2; ++side) {
                    grid_.nbor[(idim * 2 + side) * grid_.ngridmax + new_ig - 1] = get_nbor_of_coarse(grid_, ic_coarse, idim, side);
                }
            }
            int ixyz[3], idx = ic_coarse - 1; 
            ixyz[2] = idx / (grid_.nx * grid_.ny); 
            idx %= (grid_.nx * grid_.ny); 
            ixyz[1] = idx / grid_.nx; 
            ixyz[0] = idx % grid_.nx;
            for (int d = 1; d <= NDIM; ++d) grid_.xg[(d - 1) * grid_.ngridmax + new_ig - 1] = (real_t)ixyz[d - 1] + 0.5;
        } else {
            int ig = ((idc - grid_.ncoarse) % grid_.ngridmax) + 1;
            int ic = ((idc - grid_.ncoarse) / grid_.ngridmax) + 1;
            int ign[7]; grid_.get_nbor_grids(ig, ign);
            int icn_nb[6]; grid_.get_nbor_cells(ign, ic, icn_nb, ig);
            for (int i = 0; i < 6; ++i) grid_.nbor[i * grid_.ngridmax + new_ig - 1] = icn_nb[i];
            for (int d = 1; d <= NDIM; ++d) {
                real_t dx_level = 1.0 / (real_t)(1 << ilevel);
                int ix = (ic - 1) & 1, iy = ((ic - 1) & 2) >> 1, iz = ((ic - 1) & 4) >> 2; int ixyz[3] = {ix, iy, iz};
                grid_.xg[(d - 1) * grid_.ngridmax + (new_ig - 1)] = grid_.xg[(d - 1) * grid_.ngridmax + (ig - 1)] + (real_t)(ixyz[d - 1] - 0.5) * dx_level;
            }
        }
        
        int icn_coarse_nb[6] = {0};
        if (ilevel == 0) {
            for(int idim=0; idim<NDIM; ++idim) {
                icn_coarse_nb[idim*2] = get_nbor_of_coarse(grid_, ic_coarse, idim, 0);
                icn_coarse_nb[idim*2+1] = get_nbor_of_coarse(grid_, ic_coarse, idim, 1);
            }
        } else {
            int ig_coarse = ((ic_coarse - 1 - grid_.ncoarse) % grid_.ngridmax) + 1;
            int ic_pos = ((ic_coarse - 1 - grid_.ncoarse) / grid_.ngridmax) + 1;
            int ign_coarse[7]; grid_.get_nbor_grids(ig_coarse, ign_coarse);
            grid_.get_nbor_cells(ign_coarse, ic_pos, icn_coarse_nb, ig_coarse);
        }

        real_t u1[7][64] = {0}, u2[8][64] = {0};
        for(int iv=1; iv<=grid_.nvar; ++iv) u1[0][iv-1] = grid_.uold(ic_coarse, iv);
        
        for (int i = 1; i <= 6; ++i) {
            int id_n = icn_coarse_nb[i - 1];
            if (id_n > 0) {
                for (int iv = 1; iv <= grid_.nvar; ++iv) {
                    u1[i][iv - 1] = grid_.uold(id_n, iv);
                }
            } else {
                for (int iv = 1; iv <= grid_.nvar; ++iv) {
                    u1[i][iv - 1] = u1[0][iv - 1];
                }
                int ib = -id_n;
                if (ib > 0 && ib <= (int)grid_.bound_type.size() && grid_.bound_type.at(ib-1) == 1) {
                    int idim = (i - 1) / 2;
                    u1[i][1 + idim] *= -1.0;
                }
            }
        }
        
        grid_.interpol(u1, u2);
        
        for (int isc = 1; isc <= n2d; ++isc) {
            int id_child = grid_.ncoarse + (isc - 1) * grid_.ngridmax + new_ig - 1;
            for (int iv = 1; iv <= grid_.nvar; ++iv) {
                grid_.uold(id_child + 1, iv) = u2[isc-1][iv-1];
                grid_.unew(id_child + 1, iv) = u2[isc-1][iv-1];
            }
            grid_.cpu_map[id_child] = myid;
        }
        grid_.add_to_level_list(new_ig, ilevel + 1);
    }
}

void TreeUpdater::remove_grid_fine(int ilevel) {
    if (ilevel >= grid_.nlevelmax) return;
    int myid = MpiManager::instance().rank() + 1;
    int ig = grid_.get_headl(myid, ilevel + 1);
    int nkill = 0;
    while (ig > 0) {
        int next_ig = grid_.next[ig - 1];
        int n2d = (1 << NDIM);
        int father_cell = grid_.father[ig - 1];
        
        bool should_remove = false;
        if (father_cell > 0 && grid_.flag1[father_cell - 1] == 0) {
            should_remove = true;
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                if (grid_.son[idc] > 0) {
                    should_remove = false;
                    break;
                }
            }
        }
        
        if (should_remove) {
            nkill++;
            if (father_cell > 0) {
                grid_.son[father_cell - 1] = 0;
            }
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                grid_.flag1[idc] = 0;
            }
            int p = grid_.prev[ig - 1];
            int n = grid_.next[ig - 1];
            if (p > 0) grid_.next[p - 1] = n; else grid_.headl(myid, ilevel + 1) = n;
            if (n > 0) grid_.prev[n - 1] = p; else grid_.taill(myid, ilevel + 1) = p;
            grid_.numbl(myid, ilevel + 1)--;
            grid_.free_grid(ig);
        }
        ig = next_ig;
    }
    if (ramses::params::verbose) {
        RAMSES_INFO("    Entering refine_fine for level  {}", ilevel);
        RAMSES_INFO("   ==> Kill      {} sub-grids", nkill);
    }
}

void TreeUpdater::restrict_fine(int ilevel) {
    int myid = MpiManager::instance().rank() + 1;
    int ig = grid_.get_headl(myid, ilevel + 1);
    while (ig > 0) {
        int next_ig = grid_.next[ig - 1];
        int n2d = (1 << NDIM);
        int father_cell = grid_.father[ig - 1];
        
        if (father_cell > 0) {
            for (int iv = 1; iv <= grid_.nvar; ++iv) {
                real_t sum = 0;
                for (int ic = 1; ic <= n2d; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                    sum += grid_.uold(idc + 1, iv);
                }
                grid_.uold(father_cell, iv) = sum / n2d;
            }
        }
        ig = next_ig;
    }
}

void TreeUpdater::smooth_fine(int ilevel) {
    int myid = MpiManager::instance().rank() + 1, n2d = (1 << NDIM);
    int ngridmax = grid_.ngridmax;
    int n_nbor[3] = {1, 2, 2};

    for (int ismooth = 0; ismooth < NDIM; ++ismooth) {
        int thresh = n_nbor[ismooth];
        std::vector<int> ok_flag(grid_.ncell, 0);

        if (ilevel == 0) {
            for (int ic = 1; ic <= grid_.ncoarse; ++ic) {
                if (grid_.flag1[ic - 1] == 0) {
                    int icn_nb[6]; grid_.get_nbor_cells_coarse(ic, icn_nb);
                    int num_flagged_nbors = 0;
                    for (int inbor = 0; inbor < 2 * NDIM; ++inbor) {
                        int neighbor_cell = icn_nb[inbor];
                        if (neighbor_cell > 0 && neighbor_cell <= grid_.ncoarse) {
                            if (grid_.flag1[neighbor_cell - 1] == 1) num_flagged_nbors++;
                        }
                    }
                    if (num_flagged_nbors >= thresh) ok_flag[ic - 1] = 1;
                }
            }
        } else {
            int ig = grid_.get_headl(myid, ilevel);
            while (ig > 0) {
                int ign[7]; grid_.get_nbor_grids(ig, ign);
                for (int ic = 1; ic <= n2d; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * ngridmax + ig - 1;
                    if (grid_.flag1[idc] == 0) {
                        int icn_nb[6]; grid_.get_nbor_cells_exact(ign, ic, icn_nb);
                        int num_flagged_nbors = 0;
                        for (int inbor = 0; inbor < 2 * NDIM; ++inbor) {
                            int neighbor_cell = icn_nb[inbor];
                            if (neighbor_cell > 0 && neighbor_cell <= grid_.ncell) {
                                if (grid_.flag1[neighbor_cell - 1] == 1) num_flagged_nbors++;
                            }
                        }
                        if (num_flagged_nbors >= thresh) ok_flag[idc] = 1;
                    }
                }
                ig = grid_.next[ig - 1];
            }
        }

        for (int i = 0; i < (int)grid_.ncell; ++i) {
            if (ok_flag[i] == 1) grid_.flag1[i] = 1;
        }
    }
}


void TreeUpdater::init_flag(int ilevel) {}
void TreeUpdater::test_flag(int ilevel) {}

void TreeUpdater::flag_all(real_t ed, real_t ep, real_t ev, real_t eb2, const std::vector<real_t>& evar, const std::vector<int>& nexpand, int icount, int nsubcycle_val) {
    int lmax = grid_.nlevelmax;
    for (int il = lmax - 1; il >= 0; --il) {
        int nexp = nexpand[std::min(il + 1, 32) - 1];
        flag_fine(il, ed, ep, ev, eb2, evar, nexp, icount, nsubcycle_val);
    }
}

void TreeUpdater::refine_all() {
    int lmax = grid_.nlevelmax;
    for (int il = 0; il < lmax; ++il) {
        make_grid_fine(il);
    }
}

void TreeUpdater::flag_fine(int ilevel, real_t ed, real_t ep, real_t ev, real_t eb2, const std::vector<real_t>& evar, int nexp, int icount, int nsubcycle_val) {
    int lmin = config_.get_int("amr_params", "levelmin", 1), lmax = grid_.nlevelmax;
    if (ilevel >= lmax) return;
    int myid = MpiManager::instance().rank() + 1, n2d = (1 << NDIM);
    
    struct CellState {
        real_t d = 0;
        real_t v[3] = {0};
        real_t p = 0;
        real_t b2 = 0;
    };

    auto get_cell_state = [&](int cell_idx) -> CellState {
        CellState state;
        if (cell_idx <= 0 || cell_idx > grid_.ncell) return state;
        real_t gamma = grid_.gamma;
        bool is_mhd = (grid_.nvar >= 11);

        state.d = std::max(grid_.uold(cell_idx, 1), (real_t)1e-20);
        
        real_t v2 = 0.0;
        int ndim_vel = is_mhd ? 3 : NDIM;
        for (int idim = 0; idim < ndim_vel; ++idim) {
            state.v[idim] = grid_.uold(cell_idx, 2 + idim) / state.d;
            v2 += state.v[idim] * state.v[idim];
        }
        
        real_t ek = 0.5 * state.d * v2;
        real_t em = 0.0;
        real_t e_tot = 0.0;
        
        int nener = 0;
#ifdef RAMSES_NENER
        if (RAMSES_NENER > 0) nener = RAMSES_NENER;
#endif
        if (nener == 0) nener = config_.get_int("hydro_params", "nener", 0);
        real_t e_nonthermal = 0.0;
        int iener = 0;
        if (nener > 0) {
            iener = is_mhd ? (ndim_vel + 2 + 3) : (ndim_vel + 2);
            for (int ie = 0; ie < nener; ++ie) {
                real_t gamma_rad_val = (ie < (int)grid_.gamma_rad.size()) ? grid_.gamma_rad[ie] : 1.3333333333333333;
                real_t erad = grid_.uold(cell_idx, iener + 1 + ie);
                e_nonthermal += erad / (gamma_rad_val - 1.0);
            }
        }
        
        if (is_mhd) {
            real_t bx = 0.5 * (grid_.uold(cell_idx, 6) + grid_.uold(cell_idx, grid_.nvar - 2));
            real_t by = 0.5 * (grid_.uold(cell_idx, 7) + grid_.uold(cell_idx, grid_.nvar - 1));
            real_t bz = 0.5 * (grid_.uold(cell_idx, 8) + grid_.uold(cell_idx, grid_.nvar));
            state.b2 = bx*bx + by*by + bz*bz;
            em = 0.5 * state.b2;
            e_tot = grid_.uold(cell_idx, 5);
        } else {
            e_tot = grid_.uold(cell_idx, NDIM + 2);
        }
        
        real_t e_int = e_tot - ek - em - e_nonthermal;
        state.p = std::max((gamma - 1.0) * e_int, state.d * (real_t)1e-20);
        return state;
    };

    auto get_err_grad = [](real_t vl, real_t vc, real_t vr, real_t floor_val) -> real_t {
        real_t d1 = std::abs(vr - vc) / (vr + vc + floor_val);
        real_t d2 = std::abs(vc - vl) / (vc + vl + floor_val);
        return 2.0 * std::max(d1, d2);
    };

    if (ilevel == 0) {
        for (int i = 0; i < grid_.ncoarse; ++i) grid_.flag1[i] = 0;
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            grid_.cpu_map[i-1] = myid;
            if (ilevel < lmin) {
                grid_.flag1[i-1] = 1;
            } else {
                bool ok = false;
                int child_ig = grid_.son[i - 1];
                if (child_ig > 0) {
                    for (int ic_son = 1; ic_son <= n2d; ++ic_son) {
                        int idc_son = grid_.ncoarse + (ic_son - 1) * grid_.ngridmax + child_ig - 1;
                        if (grid_.son[idc_son] > 0 || grid_.flag1[idc_son] == 1) {
                            ok = true;
                            break;
                        }
                    }
                }
                if (ok) grid_.flag1[i-1] = 1;
            }
        }
    } else {
        int ig = grid_.get_headl(myid, ilevel);
        while (ig > 0) {
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                grid_.flag1[idc] = 0;
                if (ilevel < lmin) {
                    grid_.flag1[idc] = 1;
                } else {
                    bool ok = false;
                    int child_ig = grid_.son[idc];
                    if (child_ig > 0) {
                        for (int ic_son = 1; ic_son <= n2d; ++ic_son) {
                            int idc_son = grid_.ncoarse + (ic_son - 1) * grid_.ngridmax + child_ig - 1;
                            if (grid_.son[idc_son] > 0 || grid_.flag1[idc_son] == 1) {
                                ok = true;
                                break;
                            }
                        }
                    }
                    if (ok) grid_.flag1[idc] = 1;
                }
            }
            ig = grid_.next[ig - 1];
        }
    }

    if (ilevel < lmin) return;

    smooth_fine(ilevel);

    std::vector<double> r_refine = config_.get_double_array("refine_params", "r_refine");
    std::vector<double> x_refine = config_.get_double_array("refine_params", "x_refine");
    std::vector<double> y_refine = config_.get_double_array("refine_params", "y_refine");
    std::vector<double> z_refine = config_.get_double_array("refine_params", "z_refine");
    std::vector<double> exp_refine = config_.get_double_array("refine_params", "exp_refine");

    int geom_idx = ilevel - lmin;
    bool has_geom = (geom_idx >= 0 && geom_idx < (int)r_refine.size() && r_refine[geom_idx] > 0.0);
    real_t rr = has_geom ? r_refine[geom_idx] : -1.0;
    real_t xr = (has_geom && geom_idx < (int)x_refine.size()) ? x_refine[geom_idx] : 0.0;
    real_t yr = (has_geom && geom_idx < (int)y_refine.size()) ? y_refine[geom_idx] : 0.0;
    real_t zr = (has_geom && geom_idx < (int)z_refine.size()) ? z_refine[geom_idx] : 0.0;
    real_t er = (has_geom && geom_idx < (int)exp_refine.size()) ? exp_refine[geom_idx] : 2.0;

    auto check_geom = [&](real_t x, real_t y, real_t z) -> bool {
        if (!has_geom) return true;
        real_t xn = 2.0 * std::abs(x - xr) / rr;
        real_t yn = (NDIM > 1) ? (2.0 * std::abs(y - yr) / rr) : 0.0;
        real_t zn = (NDIM > 2) ? (2.0 * std::abs(z - zr) / rr) : 0.0;
        real_t r = 0.0;
        if (er < 10.0) {
            real_t sum = std::pow(xn, er);
            if (NDIM > 1) sum += std::pow(yn, er);
            if (NDIM > 2) sum += std::pow(zn, er);
            r = std::pow(sum, 1.0 / er);
        } else {
            r = std::max({xn, yn, zn});
        }
        return (r < 1.0);
    };

    auto check_gradient = [&](const CellState& sl, const CellState& sc, const CellState& sr, int idim) -> bool {
        if (ed > 0.0) {
            real_t error = get_err_grad(sl.d, sc.d, sr.d, 1e-10);
            if (error > ed) return true;
        }
        if (ep > 0.0) {
            real_t error = get_err_grad(sl.p, sc.p, sr.p, 1e-10);
            if (error > ep) return true;
        }
        if (ev > 0.0) {
            real_t vl = sl.v[idim], vc = sc.v[idim], vr = sr.v[idim];
            real_t cl = std::sqrt(std::max(grid_.gamma * sl.p / sl.d, (real_t)1e-20));
            real_t cc = std::sqrt(std::max(grid_.gamma * sc.p / sc.d, (real_t)1e-20));
            real_t cr = std::sqrt(std::max(grid_.gamma * sr.p / sr.d, (real_t)1e-20));
            real_t d1 = std::abs(vr - vc) / (cr + cc + std::abs(vr) + std::abs(vc) + 1e-10);
            real_t d2 = std::abs(vc - vl) / (cc + cl + std::abs(vc) + std::abs(vl) + 1e-10);
            real_t error = 2.0 * std::max(d1, d2);
            if (error > ev) return true;
        }
        return false;
    };

    auto get_boundary_or_cell_state = [&](int id_n, const CellState& sc, int idim) -> CellState {
        if (id_n > 0 && id_n <= grid_.ncell) {
            return get_cell_state(id_n);
        }
        CellState s = sc;
        int ib = -id_n;
        if (ib > 0 && ib <= (int)grid_.bound_type.size() && grid_.bound_type[ib - 1] == 1) {
            s.v[idim] = -sc.v[idim];
        }
        return s;
    };

    if (ed > 0.0 || ep > 0.0 || ev > 0.0 || eb2 > 0.0) {
        if (ilevel == 0) {
            for (int i = 1; i <= grid_.ncoarse; ++i) {
                if (grid_.flag1[i-1] == 0) {
                    bool ok = false;
                    CellState state_c = get_cell_state(i);
                    int icn_coarse[6] = {0};
                    grid_.get_nbor_cells_coarse(i, icn_coarse);
                    for (int idim = 0; idim < NDIM; ++idim) {
                        int id_l = icn_coarse[idim * 2];
                        int id_r = icn_coarse[idim * 2 + 1];
                        CellState state_l = get_boundary_or_cell_state(id_l, state_c, idim);
                        CellState state_r = get_boundary_or_cell_state(id_r, state_c, idim);
                        
                        if (check_gradient(state_l, state_c, state_r, idim)) {
                            ok = true; break;
                        }
                    }
                    if (ok) {
                        real_t xc[3]; grid_.get_cell_center(i, xc);
                        if (check_geom(xc[0], xc[1], xc[2])) grid_.flag1[i-1] = 1;
                    }
                }
            }
        } else {
            int ig = grid_.get_headl(myid, ilevel);
            while (ig > 0) {
                int ign[7]; grid_.get_nbor_grids(ig, ign);
                for (int ic = 1; ic <= n2d; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                    if (grid_.flag1[idc] == 0) {
                        bool ok = false;
                        CellState state_c = get_cell_state(idc + 1);
                        int icn_nb[6]; grid_.get_nbor_cells_exact(ign, ic, icn_nb);
                        for (int idim = 0; idim < NDIM; ++idim) {
                            int id_l = icn_nb[idim*2];
                            if (id_l == 0) id_l = grid_.nbor[(idim*2) * grid_.ngridmax + ig - 1];
                            int id_r = icn_nb[idim*2+1];
                            if (id_r == 0) id_r = grid_.nbor[(idim*2+1) * grid_.ngridmax + ig - 1];
                            
                            CellState state_l = get_boundary_or_cell_state(id_l, state_c, idim);
                            CellState state_r = get_boundary_or_cell_state(id_r, state_c, idim);
                            
                            if (check_gradient(state_l, state_c, state_r, idim)) {
                                ok = true; break;
                            }
                        }
                        if (ok) {
                            real_t xc[3]; grid_.get_cell_center(idc + 1, xc);
                            if (check_geom(xc[0], xc[1], xc[2])) grid_.flag1[idc] = 1;
                        }
                    }
                }
                ig = grid_.next[ig - 1];
            }
        }
    }

    for (int i = 0; i < nexp; ++i) {
        smooth_fine(ilevel);
    }

    if (ilevel > lmin && icount < nsubcycle_val) {
        ensure_ref_rules(ilevel);
    }
}


void TreeUpdater::authorize_fine(int ilevel) {
    if (ilevel > grid_.nlevelmax) return;
    int myid = MpiManager::instance().rank() + 1;
    for (int ic = 0; ic < grid_.ncell; ++ic) {
        grid_.flag2[ic] = (grid_.cpu_map[ic] == myid) ? 1 : 0;
    }
}

void TreeUpdater::ensure_ref_rules(int ilevel) {
    if (ilevel == 0) return;
    int myid = MpiManager::instance().rank() + 1;
    int ngridmax = grid_.ngridmax;
    int n2d = (1 << NDIM);
    int ig = grid_.get_headl(myid, ilevel);
    while (ig > 0) {
        int father_cell = grid_.father[ig - 1];
        bool ok = true;
        if (father_cell > 0) {
            int nbors[27];
            grid_.get_27_cell_neighbors(father_cell, nbors);
            for (int dz = -1; dz <= 1; ++dz) {
                if (NDIM < 3 && dz != 0) continue;
                for (int dy = -1; dy <= 1; ++dy) {
                    if (NDIM < 2 && dy != 0) continue;
                    for (int dx = -1; dx <= 1; ++dx) {
                        int i = (dx + 1) * 9 + (dy + 1) * 3 + (dz + 1);
                        if (nbors[i] == 0) {
                            ok = false;
                            break;
                        } else if (nbors[i] > 0) {
                            if (grid_.son[nbors[i] - 1] == 0) {
                                ok = false;
                                break;
                            }
                        }
                    }
                    if (!ok) break;
                }
                if (!ok) break;
            }
        } else {
            ok = false;
        }

        if (!ok) {
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * ngridmax + ig - 1;
                grid_.flag1[idc] = 0;
            }
        }
        ig = grid_.next[ig - 1];
    }
}

} // namespace ramses
