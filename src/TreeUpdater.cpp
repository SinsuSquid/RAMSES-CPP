#include "ramses/TreeUpdater.hpp"
#include "ramses/AmrGrid.hpp"
#include "ramses/Constants.hpp"
#include "ramses/MpiManager.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

namespace ramses {

TreeUpdater::TreeUpdater(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

static int get_nbor_of_coarse(const AmrGrid& grid, int icell, int idim, int side) {
    int idx = icell - 1, nx = grid.nx, ny = grid.ny, nz = grid.nz;
    int iz = idx / (nx * ny), iy = (idx % (nx * ny)) / nx, ix = idx % nx;
    int i = ix, j = iy, k = iz, di = 0, dj = 0, dk = 0;
    if (idim == 0) di = (side == 0) ? -1 : 1;
    else if (idim == 1) dj = (side == 0) ? -1 : 1;
    else if (idim == 2) dk = (side == 0) ? -1 : 1;
    int ni = i + di, nj = j + dj, nk = k + dk;
    bool out = (ni < 0 || ni >= nx || nj < 0 || nj >= ny || nk < 0 || nk >= nz);
    if (out && grid.nboundary > 0) {
        int ri = (ni < 0) ? -1 : (ni >= nx) ? 1 : 0, rj = (nj < 0) ? -1 : (nj >= ny) ? 1 : 0, rk = (nk < 0) ? -1 : (nk >= nz) ? 1 : 0;
        for (int ib = 1; ib <= grid.nboundary; ++ib) {
            bool match = true;
            if (grid.ibound_min[ib-1] != 0 && grid.ibound_min[ib-1] != ri) match = false;
            if (grid.ibound_max[ib-1] != 0 && grid.ibound_max[ib-1] != ri) match = false;
            if (NDIM > 1) { if (grid.jbound_min[ib-1] != 0 && grid.jbound_min[ib-1] != rj) match = false; if (grid.jbound_max[ib-1] != 0 && grid.jbound_max[ib-1] != rj) match = false; }
            if (NDIM > 2) { if (grid.kbound_min[ib-1] != 0 && grid.kbound_min[ib-1] != rk) match = false; if (grid.kbound_max[ib-1] != 0 && grid.kbound_max[ib-1] != rk) match = false; }
            if (match) return -ib;
        }
    }
    ni = (ni + nx) % nx; nj = (nj + ny) % ny; nk = (nk + nz) % nz;
    return 1 + ni + nj * nx + nk * nx * ny;
}

void TreeUpdater::make_grid_fine(int ilevel) {
    if (ilevel >= grid_.nlevelmax) return;
    int n2d = (1 << NDIM);
    std::vector<int> grids_to_refine, coarse_cells_to_refine;
    if (ilevel == 1) { for (int i = 1; i <= grid_.ncoarse; ++i) if (grid_.flag1[i-1] > 0 && grid_.son[i-1] == 0) coarse_cells_to_refine.push_back(i); }
    else {
        int myid = MpiManager::instance().rank() + 1, igrid = grid_.get_headl(myid, ilevel);
        while (igrid > 0) {
            bool refine = false; for (int ic = 1; ic <= n2d; ++ic) { if (grid_.flag1[grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid - 1] > 0 && grid_.son[grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid - 1] == 0) { refine = true; break; } }
            if (refine) grids_to_refine.push_back(igrid);
            igrid = grid_.next[igrid - 1];
        }
    }
    for (int ic_coarse : coarse_cells_to_refine) {
        int new_ig = grid_.get_free_grid(); if (new_ig == 0) return;
        grid_.son[ic_coarse - 1] = new_ig; grid_.father[new_ig - 1] = ic_coarse;
        for (int idim = 0; idim < 3; ++idim) for (int side = 0; side < 2; ++side) grid_.nbor[(idim * 2 + side) * grid_.ngridmax + new_ig - 1] = get_nbor_of_coarse(grid_, ic_coarse, idim, side);
        for (int d = 1; d <= NDIM; ++d) {
            int ixyz[3], idx = ic_coarse - 1; ixyz[2] = idx / (grid_.nx * grid_.ny); idx %= (grid_.nx * grid_.ny); ixyz[1] = idx / grid_.nx; ixyz[0] = idx % grid_.nx;
            real_t dx_coarse = grid_.boxlen / std::max({grid_.nx, grid_.ny, grid_.nz}); grid_.xg[(d-1)*grid_.ngridmax + new_ig - 1] = (ixyz[d-1] + 0.5) * dx_coarse;
        }
        for (int isc = 1; isc <= n2d; ++isc) {
            int id_child = grid_.ncoarse + (isc - 1) * grid_.ngridmax + new_ig - 1;
            for (int iv = 1; iv <= grid_.nvar; ++iv) { grid_.uold_vec[(iv-1)*grid_.ncell + id_child] = grid_.uold_vec[(iv-1)*grid_.ncell + ic_coarse - 1]; grid_.unew_vec[(iv-1)*grid_.ncell + id_child] = grid_.uold_vec[(iv-1)*grid_.ncell + ic_coarse - 1]; }
        }
        grid_.add_to_level_list(new_ig, ilevel + 1);
    }
    for (int ig : grids_to_refine) {
        int ign[7]; grid_.get_nbor_grids(ig, ign);
        for (int ic = 1; ic <= n2d; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
            if (grid_.flag1[idc] > 0 && grid_.son[idc] == 0) {
                int new_ig = grid_.get_free_grid(); if (new_ig == 0) return;
                grid_.son[idc] = new_ig; grid_.father[new_ig - 1] = idc + 1;
                int icn_nb[6]; grid_.get_nbor_cells(ign, ic, icn_nb, ig);
                for (int i = 0; i < 6; ++i) grid_.nbor[i * grid_.ngridmax + new_ig - 1] = icn_nb[i];
                for (int d = 1; d <= NDIM; ++d) {
                    real_t dx_level = grid_.boxlen / (real_t)(grid_.nx * (1 << (ilevel - 1)));
                    int ix = (ic - 1) & 1, iy = ((ic - 1) & 2) >> 1, iz = ((ic - 1) & 4) >> 2; int ixyz[3] = {ix, iy, iz};
                    grid_.xg[(d - 1) * grid_.ngridmax + (new_ig - 1)] = grid_.xg[(d - 1) * grid_.ngridmax + (ig - 1)] + (ixyz[d-1] - 0.5) * dx_level * 0.5;
                }
                real_t u1[7][20] = {0}, u2[8][20] = {0}; for(int iv=1; iv<=grid_.nvar; ++iv) u1[0][iv-1] = grid_.uold(idc+1, iv);
                for(int idim=0; idim<NDIM; ++idim) {
                    int id_l = icn_nb[idim*2], id_r = icn_nb[idim*2+1];
                    if (id_l > 0) for(int iv=1; iv<=grid_.nvar; ++iv) u1[2*idim+1][iv-1] = grid_.uold(id_l, iv); else for(int iv=1; iv<=grid_.nvar; ++iv) u1[2*idim+1][iv-1] = u1[0][iv-1];
                    if (id_r > 0) for(int iv=1; iv<=grid_.nvar; ++iv) u1[2*idim+2][iv-1] = grid_.uold(id_r, iv); else for(int iv=1; iv<=grid_.nvar; ++iv) u1[2*idim+2][iv-1] = u1[0][iv-1];
                }
                if (interpol_hook_) interpol_hook_(u1, u2); else { for(int i=0; i<8; ++i) for(int iv=0; iv<grid_.nvar; ++iv) u2[i][iv] = u1[0][iv]; }
                for (int isc = 1; isc <= n2d; ++isc) {
                    int id_child = grid_.ncoarse + (isc - 1) * grid_.ngridmax + new_ig - 1;
                    for (int iv = 1; iv <= grid_.nvar; ++iv) { grid_.uold_vec[(iv-1)*grid_.ncell + id_child] = u2[isc-1][iv-1]; grid_.unew_vec[(iv-1)*grid_.ncell + id_child] = u2[isc-1][iv-1]; }
                }
                grid_.add_to_level_list(new_ig, ilevel + 1);
            }
        }
    }
}

static void kill_recursive(AmrGrid& grid, int ig_to_kill, int ilevel) {
    if (ig_to_kill <= 0) return;
    int myid = MpiManager::instance().rank() + 1, n2d = (1 << NDIM);
    for (int ic = 1; ic <= n2d; ++ic) {
        int idc = grid.ncoarse + (ic - 1) * grid.ngridmax + ig_to_kill - 1;
        if (grid.son[idc] > 0) { kill_recursive(grid, grid.son[idc], ilevel + 1); grid.son[idc] = 0; }
    }
    int p = grid.prev[ig_to_kill - 1], n = grid.next[ig_to_kill - 1];
    if (p > 0) grid.next[p - 1] = n; else grid.headl_vec[(ilevel-1)*grid.ncpu + (myid-1)] = n;
    if (n > 0) grid.prev[n - 1] = p; else grid.taill_vec[(ilevel-1)*grid.ncpu + (myid-1)] = p;
    grid.numbl_vec[(ilevel-1)*grid.ncpu + (myid-1)]--;
    grid.father[ig_to_kill - 1] = 0; for(int i=0; i<6; ++i) grid.nbor[i * grid.ngridmax + ig_to_kill - 1] = 0;
    grid.free_grid(ig_to_kill);
}

void TreeUpdater::remove_grid_fine(int ilevel) {
    if (ilevel >= grid_.nlevelmax) return;
    int myid = MpiManager::instance().rank() + 1, igrid = grid_.get_headl(myid, ilevel), n2d = (1 << NDIM);
    std::vector<int> gk;
    while (igrid > 0) {
        for (int ic = 1; ic <= n2d; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid - 1;
            if (grid_.flag1[idc] == 0 && grid_.son[idc] > 0) { gk.push_back(grid_.son[idc]); grid_.son[idc] = 0; }
        }
        igrid = grid_.next[igrid - 1];
    }
    for (int ig : gk) kill_recursive(grid_, ig, ilevel + 1);
}

void TreeUpdater::restrict_fine(int ilevel) {
    if (ilevel >= grid_.nlevelmax) return;
    int myid = MpiManager::instance().rank() + 1, igrid = grid_.get_headl(myid, ilevel), n2d = (1 << NDIM);
    while (igrid > 0) {
        for (int ic = 1; ic <= n2d; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid - 1, ig_son = grid_.son[idc];
            if (ig_son > 0) {
                for (int iv = 1; iv <= grid_.nvar; ++iv) {
                    real_t sum = 0; for (int is = 1; is <= n2d; ++is) sum += grid_.unew(grid_.ncoarse + (is - 1) * grid_.ngridmax + ig_son, iv);
                    grid_.unew(idc + 1, iv) = sum / n2d;
                }
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

void TreeUpdater::flag_fine(int ilevel, real_t ed, real_t ep, real_t ev, real_t eb2, const std::vector<real_t>& evar) {
    int nexp = config_.get_int("amr_params", "nexpand", 0), myid = MpiManager::instance().rank() + 1, n2d = (1 << NDIM);
    int lmin = config_.get_int("amr_params", "levelmin", 1);
    
    if (ilevel < lmin) {
        if (ilevel == 1) { for (int i = 1; i <= grid_.ncoarse; ++i) grid_.flag1[i-1] = 1; }
        else {
            int ig = grid_.get_headl(myid, ilevel);
            while (ig > 0) { for (int ic = 1; ic <= n2d; ++ic) grid_.flag1[grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1] = 1; ig = grid_.next[ig - 1]; }
        }
    } else {
        int ig = grid_.get_headl(myid, ilevel);
        while (ig > 0) {
            int ign[7]; grid_.get_nbor_grids(ig, ign);
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1; grid_.flag1[idc] = 0;
                int icn[6]; grid_.get_nbor_cells(ign, ic, icn, ig);
                real_t dm = grid_.uold(idc + 1, 1), fd = 1e-10, g = grid_.gamma, v2m = 0;
                for(int j=1; j<=NDIM; ++j) { real_t v = grid_.uold(idc+1, 1+j)/dm; v2m += v*v; }
                real_t em_m = 0;
#ifdef MHD
                real_t Ax=0.5*(grid_.uold(idc+1,6)+grid_.uold(idc+1,grid_.nvar-2)), Ay=0.5*(grid_.uold(idc+1,7)+grid_.uold(idc+1,grid_.nvar-1)), Az=0.5*(grid_.uold(idc+1,8)+grid_.uold(idc+1,grid_.nvar));
                em_m = 0.5*(Ax*Ax+Ay*Ay+Az*Az);
#endif
                real_t pm = std::max((grid_.uold(idc+1, 5) - 0.5*dm*v2m - em_m)*(g-1.0), dm*fd), cm = std::sqrt(std::max(g*pm/dm, fd*fd));
                
                for (int idim = 0; idim < NDIM; ++idim) {
                    int il = icn[idim * 2], ir = icn[idim * 2 + 1];
                    real_t dl = (il > 0) ? grid_.uold(il, 1) : dm, dr = (ir > 0) ? grid_.uold(ir, 1) : dm, v2l=0, v2r=0, em_l=0, em_r=0;
                    for(int j=1; j<=NDIM; ++j) { v2l += std::pow(grid_.uold(il > 0 ? il : idc+1, 1+j)/dl, 2); v2r += std::pow(grid_.uold(ir > 0 ? ir : idc+1, 1+j)/dr, 2); }
#ifdef MHD
                    if (il > 0) { real_t Alx=0.5*(grid_.uold(il,6)+grid_.uold(il,grid_.nvar-2)), Aly=0.5*(grid_.uold(il,7)+grid_.uold(il,grid_.nvar-1)), Alz=0.5*(grid_.uold(il,8)+grid_.uold(il,grid_.nvar)); em_l = 0.5*(Alx*Alx+Aly*Aly+Alz*Alz); } else em_l = em_m;
                    if (ir > 0) { real_t Arx=0.5*(grid_.uold(ir,6)+grid_.uold(ir,grid_.nvar-2)), Ary=0.5*(grid_.uold(ir,7)+grid_.uold(ir,grid_.nvar-1)), Arz=0.5*(grid_.uold(ir,8)+grid_.uold(ir,grid_.nvar)); em_r = 0.5*(Arx*Arx+Ary*Ary+Arz*Arz); } else em_r = em_m;
#endif
                    if (ed > 0) { real_t gr = 2.0 * std::max(std::abs(dr - dm) / (dr + dm + fd), std::abs(dm - dl) / (dm + dl + fd)); if (gr > ed) { grid_.flag1[idc] = 1; break; } }
                    if (ep > 0) {
                        real_t pl = std::max((grid_.uold(il > 0 ? il : idc+1, 5) - 0.5*dl*v2l - em_l)*(g-1.0), dl*fd), pr = std::max((grid_.uold(ir > 0 ? ir : idc+1, 5) - 0.5*dr*v2r - em_r)*(g-1.0), dr*fd);
                        real_t gr = 2.0 * std::max(std::abs(pr - pm) / (pr + pm + fd), std::abs(pm - pl) / (pm + pl + fd)); if (gr > ep) { grid_.flag1[idc] = 1; break; }
                    }
                    if (eb2 > 0) { real_t gr = 2.0 * std::max(std::abs(em_r - em_m) / (em_r + em_m + fd), std::abs(em_m - em_l) / (em_m + em_l + fd)); if (gr > eb2) { grid_.flag1[idc] = 1; break; } }
                    if (ev > 0) {
                        for (int jd = 1; jd <= NDIM; ++jd) {
                            real_t vm = grid_.uold(idc+1, 1+jd)/dm, vl, vr;
                            if (il > 0) vl = grid_.uold(il, 1+jd)/dl; else { vl = vm; int ib = -il, bt = 1; if(ib>0 && ib<=(int)grid_.bound_type.size()) bt=grid_.bound_type[ib-1]; if(bt==1 && jd==idim+1) vl = -vm; }
                            if (ir > 0) vr = grid_.uold(ir, 1+jd)/dr; else { vr = vm; int ib = -ir, bt = 1; if(ib>0 && ib<=(int)grid_.bound_type.size()) bt=grid_.bound_type[ib-1]; if(bt==1 && jd==idim+1) vr = -vm; }
                            real_t cll = (il > 0) ? std::sqrt(std::max(g*(grid_.uold(il, 5)-0.5*dl*vl*vl-em_l)*(g-1.0)/dl, fd*fd)) : cm, crr = (ir > 0) ? std::sqrt(std::max(g*(grid_.uold(ir, 5)-0.5*dr*vr*vr-em_r)*(g-1.0)/dr, fd*fd)) : cm;
                            real_t gr = 2.0 * std::max(std::abs(vr - vm) / (crr + cm + std::abs(vr) + std::abs(vm) + fd), std::abs(vm - vl) / (cm + cll + std::abs(vm) + std::abs(vl) + fd)); if (gr > ev) { grid_.flag1[idc] = 1; break; }
                        }
                    }
                    if (!evar.empty()) {
                        int nh = 5;
#ifdef MHD
                        nh = 11;
#endif
                        int nener = config_.get_int("hydro_params", "nener", 0);
                        for (int iv = 0; iv < (int)evar.size(); ++iv) {
                            if (evar[iv] <= 0) continue;
                            int idx_v = nh + nener + iv + 1;
                            real_t qm = grid_.uold(idc+1, idx_v)/dm;
                            real_t ql = (il > 0) ? grid_.uold(il, idx_v)/dl : qm;
                            real_t qr = (ir > 0) ? grid_.uold(ir, idx_v)/dr : qm;
                            real_t gr = 2.0 * std::max(std::abs(qr - qm) / (std::abs(qr) + std::abs(qm) + fd), std::abs(qm - ql) / (std::abs(qm) + std::abs(ql) + fd));
                            if (gr > evar[iv]) { grid_.flag1[idc] = 1; break; }
                        }
                    }
                    if (grid_.flag1[idc] == 1) break;
                }
            }
            ig = grid_.next[ig - 1];
        }
    }
    auto smooth = [&]() {
        std::vector<int> tf;
        int head = grid_.get_headl(myid, ilevel);
        int ig_s = head;
        while (ig_s > 0) {
            int ign[7]; grid_.get_nbor_grids(ig_s, ign);
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig_s - 1;
                if (grid_.flag1[idc] == 0) {
                    int icn[6]; grid_.get_nbor_cells(ign, ic, icn, ig_s);
                    bool nbor_flagged = false;
                    for (int in = 0; in < 2 * NDIM; ++in) {
                        if (icn[in] > 0) {
                            if (grid_.flag1[icn[in]-1] == 1) { nbor_flagged = true; break; }
                        }
                    }
                    if (nbor_flagged) tf.push_back(idc);
                }
            }
            ig_s = grid_.next[ig_s - 1];
        }
        for (int idx : tf) grid_.flag1[idx] = 1;
    };
    
    // Perform expansion steps as specified in nexpand (defaulting to 1 if 0)
    int actual_nexp = (nexp > 0) ? nexp : 1;
    for (int ie = 0; ie < actual_nexp; ++ie) smooth();
}

} // namespace ramses
