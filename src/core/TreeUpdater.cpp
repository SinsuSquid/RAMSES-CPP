#include "ramses/core/TreeUpdater.hpp"
#include "ramses/core/AmrGrid.hpp"
#include "ramses/core/Constants.hpp"
#include "ramses/core/MpiManager.hpp"
#include "ramses/utils/Logger.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace ramses {

TreeUpdater::TreeUpdater(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

int get_nbor_of_coarse(const AmrGrid& grid, int ic, int idim, int side) {
    int nx = grid.nx, ny = grid.ny, nz = grid.nz;
    int idx = ic - 1;
    int iz = idx / (nx * ny); int rem = idx % (nx * ny);
    int iy = rem / nx; int ix = rem % nx;
    int ixyz[3] = {ix, iy, iz};
    int n[3] = {nx, ny, nz};
    if (side == 0) ixyz[idim] = (ixyz[idim] - 1 + n[idim]) % n[idim];
    else ixyz[idim] = (ixyz[idim] + 1) % n[idim];
    return ixyz[2] * nx * ny + ixyz[1] * nx + ixyz[0] + 1;
}

void TreeUpdater::make_grid_fine(int ilevel) {
    if (ilevel > grid_.nlevelmax) return;
    int myid = MpiManager::instance().rank() + 1, n2d = (1 << NDIM);
    std::vector<int> cells_to_refine;
    
    if (ilevel == 1) {
        for (int ic = 1; ic <= grid_.ncoarse; ++ic) {
            if (grid_.flag1[ic - 1] == 1 && grid_.flag2[ic - 1] == 1 && grid_.son[ic - 1] == 0) cells_to_refine.push_back(ic);
        }
    } else {
        int ig = grid_.get_headl(myid, ilevel - 1);
        while (ig > 0) {
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                if (grid_.flag1[idc] == 1 && grid_.flag2[idc] == 1 && grid_.son[idc] == 0) cells_to_refine.push_back(idc + 1);
            }
            ig = grid_.next[ig - 1];
        }
    }

    int ncreate = cells_to_refine.size();
    if (config_.get_bool("run_params", "verbose", false)) {
        RAMSES_INFO("    Entering refine_fine for level  {}", ilevel - 1);
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
        
        if (ilevel == 1) {
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
                real_t dx_level = 1.0 / (real_t)(1 << (ilevel - 1));
                int ix = (ic - 1) & 1, iy = ((ic - 1) & 2) >> 1, iz = ((ic - 1) & 4) >> 2; int ixyz[3] = {ix, iy, iz};
                grid_.xg[(d - 1) * grid_.ngridmax + (new_ig - 1)] = grid_.xg[(d - 1) * grid_.ngridmax + (ig - 1)] + (real_t)(ixyz[d - 1] - 0.5) * dx_level;
            }
        }
        
        int icn_coarse_nb[6] = {0};
        if (ilevel == 1) {
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
        
        int icn_ref[6] = {0};
        grid_.interpol(u1, u2);
        
        for (int isc = 1; isc <= n2d; ++isc) {
            int id_child = grid_.ncoarse + (isc - 1) * grid_.ngridmax + new_ig - 1;
            for (int iv = 1; iv <= grid_.nvar; ++iv) {
                grid_.uold(id_child + 1, iv) = u2[isc-1][iv-1];
                grid_.unew(id_child + 1, iv) = u2[isc-1][iv-1];
            }
            grid_.cpu_map[id_child] = myid;
        }
        grid_.add_to_level_list(new_ig, ilevel);
    }
}

void TreeUpdater::remove_grid_fine(int ilevel) {
    if (ilevel > grid_.nlevelmax) return;
    int myid = MpiManager::instance().rank() + 1;
    int ig = grid_.get_headl(myid, ilevel);
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
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                grid_.flag1[idc] = 0;
            }
            int p = grid_.prev[ig - 1];
            int n = grid_.next[ig - 1];
            if (p > 0) grid_.next[p - 1] = n; else grid_.headl(myid, ilevel) = n;
            if (n > 0) grid_.prev[n - 1] = p; else grid_.taill(myid, ilevel) = p;
            grid_.numbl(myid, ilevel)--;
            grid_.free_grid(ig);
        }
        ig = next_ig;
    }
    if (config_.get_bool("run_params", "verbose", false)) {
        RAMSES_INFO("    Entering refine_fine for level  {}", ilevel - 1);
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

    // flag2 is used as workspace
    for (int i = 0; i < (int)grid_.ncell; ++i) grid_.flag2[i] = 0;

    if (ilevel == 1) {
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
                if (num_flagged_nbors > 0) grid_.flag2[ic - 1] = 1;
            }
        }
    } else {
        int ig = grid_.get_headl(myid, ilevel - 1);
        while (ig > 0) {
            int ign[7]; grid_.get_nbor_grids(ig, ign);
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * ngridmax + ig - 1;
                if (grid_.flag1[idc] == 0) {
                    // Use exact same-level lookup (mirrors Fortran getnborcells):
                    // returns 0 for non-existent neighbor grids so stale flag1
                    // values at coarser levels are never counted.
                    int icn_nb[6]; grid_.get_nbor_cells_exact(ign, ic, icn_nb);
                    int num_flagged_nbors = 0;
                    for (int inbor = 0; inbor < 2 * NDIM; ++inbor) {
                        int neighbor_cell = icn_nb[inbor];
                        if (neighbor_cell > 0 && neighbor_cell <= grid_.ncell) {
                            if (grid_.flag1[neighbor_cell - 1] == 1) num_flagged_nbors++;
                        }
                    }
                    if (num_flagged_nbors > 0) grid_.flag2[idc] = 1;
                }
            }
            ig = grid_.next[ig - 1];
        }
    }

    // Apply flag2 to flag1
    for (int i = 0; i < (int)grid_.ncell; ++i) {
        if (grid_.flag2[i] == 1) grid_.flag1[i] = 1;
    }
}

void TreeUpdater::flag_fine(int ilevel, real_t ed, real_t ep, real_t ev, real_t eb2, const std::vector<real_t>& evar, int nexp, int icount, int nsubcycle_val) {
    if (ilevel > grid_.nlevelmax) return;
    int myid = MpiManager::instance().rank() + 1, n2d = (1 << NDIM), lmin = config_.get_int("amr_params", "levelmin", 1);
    
    bool verbose = config_.get_bool("run_params", "verbose", false);
    if (verbose) {
        if (ilevel == 1) RAMSES_INFO("  Entering flag_coarse");
        else RAMSES_INFO("Entering flag");
    }

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

        state.d = std::max(grid_.uold(cell_idx, 1), (real_t)1e-10);
        
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
        if (is_mhd) {
            e_tot = grid_.uold(cell_idx, 5);
            real_t bx = 0.5 * (grid_.uold(cell_idx, 6) + grid_.uold(cell_idx, grid_.nvar - 2));
            real_t by = 0.5 * (grid_.uold(cell_idx, 7) + grid_.uold(cell_idx, grid_.nvar - 1));
            real_t bz = 0.5 * (grid_.uold(cell_idx, 8) + grid_.uold(cell_idx, grid_.nvar));
            state.b2 = 0.5 * (bx*bx + by*by + bz*bz);
            em = state.b2;
            iener = 5;
        } else {
            e_tot = grid_.uold(cell_idx, NDIM + 2);
            iener = NDIM + 2;
        }
        for (int ie = 0; ie < nener; ++ie) {
            e_nonthermal += grid_.uold(cell_idx, iener + 1 + ie);
        }
        
        state.p = std::max((e_tot - ek - em - e_nonthermal) * (gamma - 1.0), state.d * (real_t)1e-10);
        return state;
    };

    auto get_err_grad = [&](real_t val_l, real_t val_c, real_t val_r, real_t floor_val) {
        return 2.0 * std::max(
            std::abs(val_r - val_c) / (val_r + val_c + floor_val),
            std::abs(val_c - val_l) / (val_c + val_l + floor_val)
        );
    };

    auto get_err_grad_u = [&](real_t vl, real_t vc, real_t vr, real_t pl, real_t pc, real_t pr, real_t dl, real_t dc, real_t dr) {
        real_t gamma = grid_.gamma;
        real_t floor_val = 1e-10;
        real_t cg = std::sqrt(std::max(gamma * pl / dl, floor_val * floor_val));
        real_t cm = std::sqrt(std::max(gamma * pc / dc, floor_val * floor_val));
        real_t cd = std::sqrt(std::max(gamma * pr / dr, floor_val * floor_val));
        return 2.0 * std::max(
            std::abs(vr - vc) / (cd + cm + std::abs(vr) + std::abs(vc) + floor_val),
            std::abs(vc - vl) / (cm + cg + std::abs(vc) + std::abs(vl) + floor_val)
        );
    };

    if (ilevel == 1) {
        for (int i = 0; i < grid_.ncoarse; ++i) grid_.flag1[i] = 0;
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            grid_.cpu_map[i-1] = myid;
            if (ilevel <= lmin) {
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
        int ig = grid_.get_headl(myid, ilevel - 1);
        while (ig > 0) {
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                grid_.flag1[idc] = 0;
                if (ilevel <= lmin) {
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

    if (ilevel >= lmin) {
        smooth_fine(ilevel);
    }

    if (ilevel >= lmin && (ed > 0.0 || ep > 0.0 || ev > 0.0 || eb2 > 0.0)) {
        if (ilevel == 1) {
            for (int i = 1; i <= grid_.ncoarse; ++i) {
                if (grid_.flag1[i-1] == 0) {
                    bool ok = false;
                    CellState state_c = get_cell_state(i);
                    for (int idim = 0; idim < NDIM; ++idim) {
                        int nx_dim = (idim == 0) ? grid_.nx : ((idim == 1) ? grid_.ny : grid_.nz);
                        int nskip = (idim == 0) ? 1 : ((idim == 1) ? grid_.nx : grid_.nx * grid_.ny);
                        int id_l = i - nskip; if (id_l < 1) id_l += nx_dim * nskip;
                        int id_r = i + nskip; if (id_r > grid_.ncoarse) id_r -= nx_dim * nskip;
                        
                        CellState state_l = get_cell_state(id_l);
                        CellState state_r = get_cell_state(id_r);
                        
                        if (ed > 0.0) {
                            real_t error = get_err_grad(state_l.d, state_c.d, state_r.d, 1e-10);
                            if (error > ed) { ok = true; break; }
                        }
                        if (ep > 0.0) {
                            real_t error = get_err_grad(state_l.p, state_c.p, state_r.p, 1e-10);
                            if (error > ep) { ok = true; break; }
                        }
                        if (ev > 0.0) {
                            real_t error = get_err_grad_u(state_l.v[idim], state_c.v[idim], state_r.v[idim],
                                                          state_l.p, state_c.p, state_r.p,
                                                          state_l.d, state_c.d, state_r.d);
                            if (error > ev) { ok = true; break; }
                        }
                        if (eb2 > 0.0) {
                            real_t error = get_err_grad(state_l.b2, state_c.b2, state_r.b2, 1e-10);
                            if (error > eb2) { ok = true; break; }
                        }
                    }
                    if (ok) grid_.flag1[i-1] = 1;
                }
            }
        } else {
            int ig = grid_.get_headl(myid, ilevel - 1);
            while (ig > 0) {
                int ign[7]; grid_.get_nbor_grids(ig, ign);
                for (int ic = 1; ic <= n2d; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                    if (grid_.flag1[idc] == 0) {
                        bool ok = false;
                        CellState state_c = get_cell_state(idc + 1);
                        int icn_nb[6]; grid_.get_nbor_cells(ign, ic, icn_nb, ig);
                        for (int idim = 0; idim < NDIM; ++idim) {
                            int id_l = icn_nb[idim*2];
                            int id_r = icn_nb[idim*2+1];
                            CellState state_l = (id_l > 0 && id_l <= grid_.ncell) ? get_cell_state(id_l) : state_c;
                            CellState state_r = (id_r > 0 && id_r <= grid_.ncell) ? get_cell_state(id_r) : state_c;
                            
                            if (ed > 0.0) {
                                real_t error = get_err_grad(state_l.d, state_c.d, state_r.d, 1e-10);
                                if (error > ed) { ok = true; break; }
                            }
                            if (ep > 0.0) {
                                real_t error = get_err_grad(state_l.p, state_c.p, state_r.p, 1e-10);
                                if (error > ep) { ok = true; break; }
                            }
                            if (ev > 0.0) {
                                real_t error = get_err_grad_u(state_l.v[idim], state_c.v[idim], state_r.v[idim],
                                                              state_l.p, state_c.p, state_r.p,
                                                              state_l.d, state_c.d, state_r.d);
                                if (error > ev) { ok = true; break; }
                            }
                            if (eb2 > 0.0) {
                                real_t error = get_err_grad(state_l.b2, state_c.b2, state_r.b2, 1e-10);
                                if (error > eb2) { ok = true; break; }
                            }
                        }
                        if (ok) grid_.flag1[idc] = 1;
                    }
                }
                ig = grid_.next[ig - 1];
            }
        }
    }

    for (int i = 0; i < nexp; ++i) {
        smooth_fine(ilevel);
    }

    authorize_fine(ilevel);
    if (ilevel > lmin + 1 && icount < nsubcycle_val) {
        ensure_ref_rules(ilevel);
    }
    
    if (verbose) {
        if (ilevel > 1) RAMSES_INFO("Complete flag");
    }

    if (config_.get_bool("run_params", "verbose", false)) {
        int nflagged = 0;
        if (ilevel == 1) {
            for (int i = 0; i < grid_.ncoarse; ++i) {
                if (grid_.flag1[i] == 1) nflagged++;
            }
        } else {
            int ig = grid_.get_headl(myid, ilevel - 1);
            while (ig > 0) {
                for (int ic = 1; ic <= n2d; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                    if (grid_.flag1[idc] == 1) nflagged++;
                }
                ig = grid_.next[ig - 1];
            }
        }
        RAMSES_INFO("   ==> Flag {} cells (C++ level {})", nflagged, ilevel - 1);
    }
}


void TreeUpdater::authorize_fine(int ilevel) {
    if (ilevel >= grid_.nlevelmax) return;
    int myid = MpiManager::instance().rank() + 1;
    for (int ic = 0; ic < grid_.ncell; ++ic) {
        grid_.flag2[ic] = (grid_.cpu_map[ic] == myid) ? 1 : 0;
    }
}

void TreeUpdater::ensure_ref_rules(int ilevel) {
    int myid = MpiManager::instance().rank() + 1;
    int ngridmax = grid_.ngridmax;
    int ig = grid_.get_headl(myid, ilevel - 1);
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
            int n2d = (1 << NDIM);
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * ngridmax + ig - 1;
                grid_.flag1[idc] = 0;
            }
        }
        ig = grid_.next[ig - 1];
    }
}

} // namespace ramses
