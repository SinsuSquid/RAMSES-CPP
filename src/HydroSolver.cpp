#include "ramses/HydroSolver.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/RiemannSolver.hpp"
#include "ramses/Parameters.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace ramses {

void HydroSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    std::vector<int> active_octs;
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        active_octs.push_back(igrid);
        igrid = grid_.next[igrid - 1];
    }
    if (!active_octs.empty()) godfine1(active_octs, ilevel, dt, dx);
}

void HydroSolver::set_unew(int ilevel) {
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            // DEBUG: Check idc bounds
            if (idc < 1 || idc > grid_.ncell) std::cerr<<"[set_unew] idc out of bounds: "<<idc<<" ncell="<<grid_.ncell<<std::endl;
            for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.unew(idc, iv) = grid_.uold(idc, iv);
        }
        igrid = grid_.next[igrid - 1];
    }
}

void HydroSolver::set_uold(int ilevel) {
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.uold(idc, iv) = grid_.unew(idc, iv);
        }
        igrid = grid_.next[igrid - 1];
    }
}

void HydroSolver::godfine1(const std::vector<int>& octs, int ilevel, real_t dt, real_t dx) {
    real_t gamma = grid_.gamma;
    real_t dt_dx = dt / dx;

    for (int igrid : octs) {
        gather_stencil(igrid, ilevel, *stencil_ptr_);
        
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son[idc] != 0) continue;

            int cz = (NDIM > 2) ? (ic - 1) / 4 : 0;
            int cy = (NDIM > 1) ? ((ic - 1) % 4) / 2 : 0;
            int cx = (ic - 1) % 2;
            int sz = (NDIM > 2 ? 2 : 0) + cz;
            int sy = (NDIM > 1 ? 2 : 0) + cy;
            int sx = 2 + cx;

            real_t un[5]; 
            for(int iv=0; iv<5; ++iv) un[iv] = grid_.uold(idc, iv+1);

            for (int idim = 0; idim < NDIM; ++idim) {
                for (int iside = 0; iside < 2; ++iside) {
                    int sxl = sx, syl = sy, szl = sz;
                    int sxr = sx, syr = sy, szr = sz;
                    if (idim == 0) { if (iside == 0) sxl--; else sxr++; }
                    if (idim == 1) { if (iside == 0) syl--; else syr++; }
                    if (idim == 2) { if (iside == 0) szl--; else szr++; }

                    int szL = szl, syL = syl, sxL = sxl;
                    int szR = szr, syR = syr, sxR = sxr;
                    int szLL = szL, syLL = syL, sxLL = sxL;
                    int szRR = szR, syRR = syR, sxRR = sxR;

                    if (idim == 0) { sxLL--; sxRR++; }
                    if (idim == 1) { syLL--; syRR++; }
                    if (idim == 2) { szLL--; szRR++; }

                    real_t qL[5], qR[5], qLL[5], qRR[5];
                    ctoprim(stencil_ptr_->uloc[szL][syL][sxL], qL, gamma);
                    ctoprim(stencil_ptr_->uloc[szR][syR][sxR], qR, gamma);
                    ctoprim(stencil_ptr_->uloc[szLL][syLL][sxLL], qLL, gamma);
                    ctoprim(stencil_ptr_->uloc[szRR][syRR][sxRR], qRR, gamma);

                    auto minmod = [](real_t a, real_t b) {
                        if (a * b <= 0.0) return 0.0;
                        return (std::abs(a) < std::abs(b)) ? a : b;
                    };

                    real_t ql[5], qr[5], flux[5];
                    for (int v = 0; v < 5; ++v) {
                        real_t dqL = minmod(qL[v] - qLL[v], qR[v] - qL[v]);
                        real_t dqR = minmod(qR[v] - qL[v], qRR[v] - qR[v]);
                        ql[v] = qL[v] + 0.5 * dqL;
                        qr[v] = qR[v] - 0.5 * dqR;
                    }

                    if (idim == 1) { std::swap(ql[1], ql[2]); std::swap(qr[1], qr[2]); }
                    if (idim == 2) { std::swap(ql[1], ql[3]); std::swap(qr[1], qr[3]); }

                    RiemannSolver::solve_hllc(ql, qr, flux, gamma);

                    if (idim == 1) std::swap(flux[1], flux[2]);
                    if (idim == 2) std::swap(flux[1], flux[3]);

                    real_t sign = (iside == 0) ? 1.0 : -1.0;
                    for(int iv=0; iv<5; ++iv) un[iv] += sign * flux[iv] * dt_dx;
                }
            }
            for (int iv = 1; iv <= 5; ++iv) grid_.unew(idc, iv) = un[iv-1];
        }
    }
}

void HydroSolver::ctoprim(const real_t u[], real_t q[], real_t gamma) {
    const real_t smallr = 1e-10;
    real_t d = std::max(u[0], smallr);
    q[0] = d;
    real_t vel2 = 0.0;
    for (int i = 1; i <= 3; ++i) {
        q[i] = (i <= NDIM) ? u[i] / d : 0.0;
        vel2 += q[i] * q[i];
    }
    real_t e_int = u[4] - 0.5 * d * vel2;
    q[4] = std::max(e_int * (gamma - 1.0), d * 1e-10);
}

void HydroSolver::gather_stencil(int igrid, int ilevel, LocalStencil& stencil) {
    int nbors_father[27] = {0};
    grid_.get_3x3x3_father(igrid, nbors_father);
    int myid = MpiManager::instance().rank() + 1;

    for (int i = 0; i < constants::threetondim; ++i) {
        int ifather = nbors_father[i];
        
        int fz = (NDIM > 2) ? i / 9 : 0;
        int fy = (NDIM > 1) ? (i % (NDIM > 2 ? 9 : 9)) / 3 : 0; 
        if (NDIM == 2) fy = i / 3;
        int fx = i % 3;
        
        if (ifather > 0 && grid_.cpu_map[ifather] == myid) {
            if (grid_.son[ifather] > 0) {
                int ig = grid_.son[ifather];
                for (int ic = 1; ic <= constants::twotondim; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
                    int cz = (NDIM > 2) ? (ic - 1) / 4 : 0;
                    int cy = (NDIM > 1) ? ((ic - 1) % 4) / 2 : 0;
                    int cx = (ic - 1) % 2;
                    int sx = fx * 2 + cx, sy = fy * 2 + cy, sz = fz * 2 + cz;
                    for (int iv = 1; iv <= grid_.nvar; ++iv) stencil.uloc[sz][sy][sx][iv - 1] = grid_.uold(idc, iv);
                    stencil.refined[sz][sy][sx] = true;
                }
            } else {
                for (int cz = 0; cz < (NDIM > 2 ? 2 : 1); ++cz) 
                for (int cy = 0; cy < (NDIM > 1 ? 2 : 1); ++cy) 
                for (int cx = 0; cx < 2; ++cx) {
                    int sx = fx * 2 + cx, sy = fy * 2 + cy, sz = fz * 2 + cz;
                    for (int iv = 1; iv <= grid_.nvar; ++iv) stencil.uloc[sz][sy][sx][iv - 1] = grid_.uold(ifather, iv);
                    stencil.refined[sz][sy][sx] = false;
                }
            }
        }
    }
}

real_t HydroSolver::compute_courant_step(int ilevel_ignored, real_t dx_coarse, real_t gamma, real_t courant_factor) {
    real_t dt_max = 1e30; const real_t smallr = 1e-10;
    int myid = MpiManager::instance().rank() + 1;
    
    for (int il = 1; il <= grid_.nlevelmax; ++il) {
        real_t dx = dx_coarse / static_cast<real_t>(1 << (il - 1));
        int igrid = grid_.headl(myid, il);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int i = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                if (grid_.son[i] != 0) continue; 
                real_t d = std::max(grid_.uold(i, 1), smallr);
                real_t vel2 = 0.0, vel_max = 0.0;
                for (int idim = 1; idim <= NDIM; ++idim) {
                    real_t v = grid_.uold(i, 1 + idim) / d;
                    vel2 += v * v; vel_max = std::max(vel_max, std::abs(v));
                }
                real_t e_int = grid_.uold(i, 5) - 0.5 * d * vel2;
                real_t p = std::max(e_int * (gamma - 1.0), d * 1e-10);
                dt_max = std::min(dt_max, courant_factor * dx / (vel_max + std::sqrt(gamma * p / d)));
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    return dt_max;
}

void HydroSolver::get_diagnostics(int ilevel, real_t dx, real_t& min_d, real_t& max_v, real_t& min_t, real_t& max_t) {
    min_d = 1e30; max_v = 0.0; min_t = 1e30; max_t = 0.0;
    int myid = MpiManager::instance().rank() + 1;
    real_t gamma = grid_.gamma;
    real_t mu = 1.4, kB = 1.3806e-16, mH = 1.67e-24;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int id = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son[id] != 0) continue;
            real_t d = grid_.uold(id, 1); min_d = std::min(min_d, d);
            real_t v2 = 0; for(int d_idx=1; d_idx<=NDIM; ++d_idx) v2 += std::pow(grid_.uold(id, 1+d_idx)/d, 2);
            max_v = std::max(max_v, std::sqrt(v2));
            real_t eint = grid_.uold(id, 5) - 0.5 * d * v2;
            real_t T = eint * params::units_pressure / (d * params::units_density) * (mu * mH / kB) * (gamma - 1.0);
            min_t = std::min(min_t, T); max_t = std::max(max_t, T);
        }
        igrid = grid_.next[igrid - 1];
    }
}

void HydroSolver::add_gravity_source_terms(int ilevel, real_t dt) {
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int id = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            real_t d = grid_.uold(id, 1);
            real_t e_kin_old = 0; for(int i=1; i<=NDIM; ++i) e_kin_old += 0.5 * std::pow(grid_.uold(id, 1+i), 2) / d;
            for (int idim = 1; idim <= NDIM; ++idim) grid_.uold(id, 1 + idim) += d * grid_.f(id, idim) * dt;
            real_t e_kin_new = 0; for(int i=1; i<=NDIM; ++i) e_kin_new += 0.5 * std::pow(grid_.uold(id, 1+i), 2) / d;
            grid_.uold(id, 5) += (e_kin_new - e_kin_old);
        }
        igrid = grid_.next[igrid - 1];
    }
}

} // namespace ramses
