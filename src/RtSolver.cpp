#include "ramses/RtSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/MpiManager.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace ramses {

#define SSIZE 6
#define FSIZE 4

void RtSolver::initialize() {
    nGroups = 0;
#ifdef RAMSES_NGROUPS
    nGroups = RAMSES_NGROUPS;
#endif
    nGroups = config_.get_int("rt_params", "nGroups", nGroups);
    rt_c_speed = config_.get_double("rt_params", "rt_c", 1.0);
    rt_use_hll = config_.get("rt_params", "rt_riemann", "hll") == "hll";

    if (nGroups > 0) {
        load_hll_eigenvalues();
        chem_ = std::make_unique<RtChemistry>(nGroups);
    }
}

void RtSolver::load_hll_eigenvalues() {
    std::string path = config_.get("rt_params", "hll_evals_file", "legacy/rt/hll_evals.list");
    std::ifstream file(path);
    if (!file.is_open()) {
        // Try other common relative paths
        std::vector<std::string> alt_paths = {"../legacy/rt/hll_evals.list", "../../legacy/rt/hll_evals.list", "../../../legacy/rt/hll_evals.list"};
        for (const auto& alt : alt_paths) {
            file.open(alt);
            if (file.is_open()) {
                path = alt;
                break;
            }
            file.clear();
        }
    }

    if (!file.is_open()) {
        std::cerr << "[RtSolver] Error: Could not open HLL eigenvalues file: " << path << " (checked common relative paths)" << std::endl;
        return;
    }

    int n; file >> n;
    lambda1.assign(101, std::vector<real_t>(101, 0.0));
    lambda4.assign(101, std::vector<real_t>(101, 0.0));

    for (int i = 0; i < 101; ++i) {
        for (int j = 0; j < 101; ++j) {
            int ii, jj;
            real_t d1, d2;
            file >> ii >> jj >> lambda1[ii][jj] >> d1 >> d2 >> lambda4[ii][jj];
        }
    }
    std::cout << "[RtSolver] Loaded HLL eigenvalues from " << path << std::endl;
}

void RtSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    if (nGroups <= 0) return;
    std::vector<int> active_octs;
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        active_octs.push_back(igrid);
        igrid = grid_.next[igrid - 1];
    }
    if (!active_octs.empty()) rt_godfine1(active_octs, ilevel, dt, dx);
}

void RtSolver::rt_godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx) {
    int ncache = ind_grid.size();
    int nvar_hydro = 5; 
#ifdef MHD
    nvar_hydro = 11;
#endif
    int nrtvar = nGroups * (1 + NDIM) + nIons_;

    int s_vol = 1;
    for(int d=0; d<NDIM; ++d) s_vol *= SSIZE;
    int f_vol = 1;
    for(int d=0; d<NDIM; ++d) f_vol *= FSIZE;

    uloc.assign((size_t)ncache * s_vol * nrtvar, 0.0);
    cFlx.assign((size_t)ncache * s_vol * (NDIM + 1) * NDIM, 0.0);
    lmin.assign((size_t)ncache * s_vol * NDIM, 0.0);
    lmax.assign((size_t)ncache * s_vol * NDIM, 0.0);
    flux.assign((size_t)ncache * f_vol * nrtvar * NDIM, 0.0);
    ok_stencil.assign((size_t)ncache * s_vol, false);

    auto get_u_idx = [&](int n, int i, int j, int k, int iv) {
        int idx = i;
        if (NDIM >= 2) idx += j * SSIZE;
        if (NDIM >= 3) idx += k * SSIZE * SSIZE;
        return ((size_t)n * s_vol + idx) * nrtvar + (iv - 1);
    };

    auto get_s_idx = [&](int n, int i, int j, int k) {
        int idx = i;
        if (NDIM >= 2) idx += j * SSIZE;
        if (NDIM >= 3) idx += k * SSIZE * SSIZE;
        return (size_t)n * s_vol + idx;
    };

    // 1. Gather stencil
    for (int n = 0; n < ncache; ++n) {
        int ig = ind_grid[n];
        int nbors[27];
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
            grid_.get_27_cell_neighbors(idc, nbors);
            
            int ix0 = (ic-1) % 2 + 2;
            int iy0 = (NDIM >= 2) ? ((ic-1) / 2) % 2 + 2 : 0;
            int iz0 = (NDIM >= 3) ? (ic-1) / 4 + 2 : 0;

            for (int kn = 0; kn < (NDIM >= 3 ? 3 : 1); ++kn) {
                for (int jn = 0; jn < (NDIM >= 2 ? 3 : 1); ++jn) {
                    for (int in = 0; in < 3; ++in) {
                        int n_idx = (kn * 3 + jn) * 3 + in;
                        int idn = nbors[n_idx];
                        if (idn > 0) {
                            int ix = ix0 + in - 1;
                            int iy = (NDIM >= 2) ? iy0 + jn - 1 : 0;
                            int iz = (NDIM >= 3) ? iz0 + kn - 1 : 0;
                            for (int iv = 1; iv <= nrtvar; ++iv) {
                                uloc[get_u_idx(n, ix, iy, iz, iv)] = grid_.uold(idn, nvar_hydro + iv);
                            }
                            if (grid_.son[idn] != 0) ok_stencil[get_s_idx(n, ix, iy, iz)] = true;
                        }
                    }
                }
            }
        }
    }

    // 2. Compute fluxes for each group
    for (int igr = 0; igr < nGroups; ++igr) {
        int iP0 = nIons_ + 1 + igr * (1 + NDIM);
        cmp_rt_faces(ncache, iP0, ilevel, dt, dx);
    }

    // 3. Update state
    for (int n = 0; n < ncache; ++n) {
        int igrid = ind_grid[n];
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son[idc] != 0) continue;

            int ix = (ic-1) % 2 + 2;
            int iy = (NDIM >= 2) ? ((ic-1) / 2) % 2 + 2 : 0;
            int iz = (NDIM >= 3) ? (ic-1) / 4 + 2 : 0;

            for (int iv = 1; iv <= nrtvar; ++iv) {
                real_t dflx = 0;
                for (int idim = 0; idim < NDIM; ++idim) {
                    int i0=0, j0=0, k0=0;
                    if(idim==0) i0=1; else if(idim==1) j0=1; else k0=1;
                    
                    auto get_f_idx = [&](int nn, int ii, int jj, int kk, int ivv, int idim_idx) {
                        int fi = ii - 2, fj = (NDIM >= 2 ? jj - 2 : 0), fk = (NDIM >= 3 ? kk - 2 : 0);
                        int i_off = fi;
                        if (NDIM >= 2) i_off += fj * FSIZE;
                        if (NDIM >= 3) i_off += fk * FSIZE * FSIZE;
                        return ((size_t)nn * f_vol * nrtvar + (size_t)i_off * nrtvar + (ivv - 1)) * NDIM + idim_idx;
                    };
                    
                    dflx += flux[get_f_idx(n, ix, iy, iz, iv, idim)] - flux[get_f_idx(n, ix+i0, iy+j0, iz+k0, iv, idim)];
                }
                grid_.unew(idc, nvar_hydro + iv) += dflx;
            }
        }
    }
}

void RtSolver::cmp_flux_tensors(int ncache, int iP0, int ilevel) {
    int s_vol = 1; for(int d=0; d<NDIM; ++d) s_vol *= SSIZE;
    int nrtvar = nGroups * (1 + NDIM) + nIons_;
    real_t c2 = rt_c_speed * rt_c_speed;

    for (int n = 0; n < ncache; ++n) {
        for (int k = 0; k < (NDIM >= 3 ? SSIZE : 1); ++k) {
            for (int j = 0; j < (NDIM >= 2 ? SSIZE : 1); ++j) {
                for (int i = 0; i < SSIZE; ++i) {
                    int idx_stencil = i;
                    if (NDIM >= 2) idx_stencil += j * SSIZE;
                    if (NDIM >= 3) idx_stencil += k * SSIZE * SSIZE;
                    size_t idx = ((size_t)n * s_vol + idx_stencil);
                    
                    real_t Np = uloc[idx * nrtvar + (iP0 - 1)];
                    if (Np < 0) Np = smallNp;
                    
                    real_t Fp[3] = {0,0,0};
                    real_t Fp_sq = 0;
                    for(int d=0; d<NDIM; ++d) {
                        Fp[d] = uloc[idx * nrtvar + (iP0 + d)];
                        Fp_sq += Fp[d] * Fp[d];
                    }

                    size_t c_idx = idx * (NDIM + 1) * NDIM;
                    
                    if (Np > smallNp) {
                        real_t red_f_sq = Fp_sq / (rt_c_speed * rt_c_speed * Np * Np);
                        red_f_sq = std::min(red_f_sq, 1.0);
                        real_t chi = (3.0 + 4.0 * red_f_sq) / (5.0 + 2.0 * std::sqrt(std::max(0.0, 4.0 - 3.0 * red_f_sq)));
                        
                        real_t iterm = (1.0 - chi) / 2.0;
                        real_t oterm = (3.0 * chi - 1.0) / 2.0;
                        
                        real_t u_vec[3] = {0,0,0};
                        if (Fp_sq > 0) {
                            real_t Fp_mag = std::sqrt(Fp_sq);
                            for(int d=0; d<NDIM; ++d) u_vec[d] = Fp[d] / Fp_mag;
                        }

                        for (int id = 0; id < NDIM; ++id) {
                            size_t base = c_idx + id * (NDIM + 1);
                            cFlx[base] = Fp[id];
                            for (int kk = 0; kk < NDIM; ++kk) {
                                real_t tensor_val = oterm * u_vec[kk] * u_vec[id];
                                if (kk == id) tensor_val += iterm;
                                cFlx[base + 1 + kk] = tensor_val * c2 * Np;
                            }
                        }
                    } else {
                        for (int id = 0; id < NDIM; ++id) {
                            size_t base = c_idx + id * (NDIM + 1);
                            cFlx[base] = 0.0;
                            for (int kk = 0; kk < NDIM; ++kk) {
                                cFlx[base + 1 + kk] = (kk == id) ? (c2 * Np / 3.0) : 0.0;
                            }
                        }
                    }
                }
            }
        }
    }
}

void RtSolver::cmp_eigenvals(int ncache, int iP0, int ilevel) {
    int s_vol = 1; for(int d=0; d<NDIM; ++d) s_vol *= SSIZE;
    int nrtvar = nGroups * (1 + NDIM) + nIons_;
    for (int n = 0; n < ncache; ++n) {
        for (int k = 0; k < (NDIM >= 3 ? SSIZE : 1); ++k) {
            for (int j = 0; j < (NDIM >= 2 ? SSIZE : 1); ++j) {
                for (int i = 0; i < SSIZE; ++i) {
                    int idx_stencil = i;
                    if (NDIM >= 2) idx_stencil += j * SSIZE;
                    if (NDIM >= 3) idx_stencil += k * SSIZE * SSIZE;
                    size_t idx = ((size_t)n * s_vol + idx_stencil);
                    
                    real_t Np = uloc[idx * nrtvar + (iP0 - 1)];
                    real_t Fp[3] = {0,0,0};
                    real_t Fp_mag_sq = 0;
                    for(int d=0; d<NDIM; ++d) {
                        Fp[d] = uloc[idx * nrtvar + (iP0 + d)];
                        Fp_mag_sq += Fp[d] * Fp[d];
                    }
                    real_t Fp_mag = std::sqrt(Fp_mag_sq);
                    
                    real_t ff = (Np > smallNp) ? Fp_mag / (rt_c_speed * Np) : 0.0;
                    ff = std::max(0.0, std::min(1.0, ff));
                    
                    for (int id = 0; id < NDIM; ++id) {
                        real_t omega = (Fp_mag > 0) ? Fp[id] / Fp_mag : 0.0;
                        omega = std::max(-1.0, std::min(1.0, omega));
                        inp_eigenvals(ff, omega, lmin[idx * NDIM + id], lmax[idx * NDIM + id]);
                    }
                }
            }
        }
    }
}

void RtSolver::inp_eigenvals(real_t ff, real_t omega, real_t& lm, real_t& lp) {
    real_t theta = std::acos(omega);
    real_t lff = ff * 100.0;
    real_t ltt = theta / M_PI * 100.0;
    
    int ii = std::min((int)lff, 99);
    int jj = std::min((int)ltt, 99);
    real_t dd1 = lff - ii;
    real_t dd2 = ltt - jj;
    real_t de1 = 1.0 - dd1;
    real_t de2 = 1.0 - dd2;
    
    lm = de1*de2*lambda1[ii][jj] + dd1*de2*lambda1[ii+1][jj] + de1*dd2*lambda1[ii][jj+1] + dd1*dd2*lambda1[ii+1][jj+1];
    lp = de1*de2*lambda4[ii][jj] + dd1*de2*lambda4[ii+1][jj] + de1*dd2*lambda4[ii][jj+1] + dd1*dd2*lambda4[ii+1][jj+1];
}

void RtSolver::cmp_rt_faces(int ncache, int iP0, int ilevel, real_t dt, real_t dx) {
    cmp_flux_tensors(ncache, iP0, ilevel);
    if (rt_use_hll) cmp_eigenvals(ncache, iP0, ilevel);
    
    real_t dtdx = dt / dx;
    int s_vol = 1; for(int d=0; d<NDIM; ++d) s_vol *= SSIZE;
    int f_vol = 1; for(int d=0; d<NDIM; ++d) f_vol *= FSIZE;
    int nrtvar = nGroups * (1 + NDIM) + nIons_;

    for (int idim = 0; idim < NDIM; ++idim) {
        for (int k = (NDIM >= 3 ? 2 : 0); k < (NDIM >= 3 ? 6 : 1); ++k) {
            for (int j = (NDIM >= 2 ? 2 : 0); j < (NDIM >= 2 ? 6 : 1); ++j) {
                for (int i = 2; i < 6; ++i) {
                    for (int n = 0; n < ncache; ++n) {
                        int i0=0, j0=0, k0=0;
                        if(idim==0) i0=1; else if(idim==1) j0=1; else k0=1;
                        
                        auto get_s_idx_local = [&](int nn, int ii, int jj, int kk) {
                            int offset = ii;
                            if (NDIM >= 2) offset += jj * SSIZE;
                            if (NDIM >= 3) offset += kk * SSIZE * SSIZE;
                            return (size_t)nn * s_vol + offset;
                        };

                        size_t idx_dn = get_s_idx_local(n, i - i0, j - j0, k - k0);
                        size_t idx_up = get_s_idx_local(n, i, j, k);

                        real_t fdn[4], fup[4], udn[4], uup[4]; // Max 3+1 = 4
                        size_t c_idx_dn = idx_dn * (NDIM + 1) * NDIM;
                        size_t c_idx_up = idx_up * (NDIM + 1) * NDIM;
                        
                        for(int d=0; d<NDIM+1; ++d) {
                            fdn[d] = cFlx[c_idx_dn + idim * (NDIM + 1) + d];
                            fup[d] = cFlx[c_idx_up + idim * (NDIM + 1) + d];
                            udn[d] = uloc[idx_dn * nrtvar + (iP0 - 1 + d)];
                            uup[d] = uloc[idx_up * nrtvar + (iP0 - 1 + d)];
                        }

                        real_t lminus = 0, lplus = 0;
                        if (rt_use_hll) {
                            lminus = std::min({lmin[idx_dn * NDIM + idim], lmin[idx_up * NDIM + idim], 0.0});
                            lplus = std::max({lmax[idx_dn * NDIM + idim], lmax[idx_up * NDIM + idim], 0.0});
                        }

                        real_t face_flx[4];
                        if (rt_use_hll) {
                            real_t div = lplus - lminus;
                            if (std::abs(div) < 1e-20) {
                                for(int d=0; d<NDIM+1; ++d) face_flx[d] = 0;
                            } else {
                                for(int d=0; d<NDIM+1; ++d) 
                                    face_flx[d] = (lplus * fdn[d] - lminus * fup[d] + lplus * lminus * rt_c_speed * (uup[d] - udn[d])) / div;
                            }
                        } else {
                            for(int d=0; d<NDIM+1; ++d)
                                face_flx[d] = 0.5 * (fdn[d] + fup[d] - rt_c_speed * (uup[d] - udn[d]));
                        }

                        if (i >= 2 && i <= 5 && j >= (NDIM>=2?2:0) && j <= (NDIM>=2?5:0) && k >= (NDIM>=3?2:0) && k <= (NDIM>=3?5:0)) {
                            int fi = i - 2, fj = (NDIM >= 2 ? j - 2 : 0), fk = (NDIM >= 3 ? k - 2 : 0);
                            int i_off = fi;
                            if (NDIM >= 2) i_off += fj * FSIZE;
                            if (NDIM >= 3) i_off += fk * FSIZE * FSIZE;
                            for(int d=0; d<NDIM+1; ++d) {
                                flux[((size_t)n * f_vol * nrtvar + (size_t)i_off * nrtvar + (iP0 - 1 + d)) * NDIM + idim] = face_flx[d] * dtdx;
                            }
                        }
                    }
                }
            }
        }
    }
}

void RtSolver::apply_source_terms(int ilevel, real_t dt) {
    if (nGroups <= 0 || !chem_) return;
    int myid = MpiManager::instance().rank() + 1;
    int nvar_hydro = 5; 
#ifdef MHD
    nvar_hydro = 11;
#endif
    int iIons = nvar_hydro + 1;

    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son[idc] != 0) continue;

            real_t nH = grid_.uold(idc, 1) * params::units_density / 1.67e-24;
            real_t T2 = (grid_.uold(idc, 5) - 0.5 * grid_.uold(idc, 2)*grid_.uold(idc, 2)/grid_.uold(idc, 1)) * (grid_.gamma - 1.0) / grid_.uold(idc, 1) * params::units_pressure / params::units_density * 1.67e-24 / 1.38e-16;

            real_t xion[3] = {0,0,0};
            for(int i=0; i<std::min(nIons_, 3); ++i) xion[i] = grid_.uold(idc, iIons + i) / grid_.uold(idc, 1);
            
            real_t Np[10], Fp[10][3];
            for (int ig = 0; ig < std::min(nGroups, 10); ++ig) {
                int iP0_abs = iIons + nIons_ + ig * (1+NDIM);
                Np[ig] = grid_.uold(idc, iP0_abs + 1);
                for(int d=0; d<NDIM; ++d) Fp[ig][d] = grid_.uold(idc, iP0_abs + 2 + d);
            }

            chem_->solve_chemistry(T2, xion, Np, Fp, nH, dt * params::units_time, 1.0);

            for(int i=0; i<std::min(nIons_, 3); ++i) grid_.uold(idc, iIons + i) = xion[i] * grid_.uold(idc, 1);
            for (int ig = 0; ig < std::min(nGroups, 10); ++ig) {
                int iP0_abs = iIons + nIons_ + ig * (1+NDIM);
                grid_.uold(idc, iP0_abs + 1) = Np[ig];
                for(int d=0; d<NDIM; ++d) grid_.uold(idc, iP0_abs + 2 + d) = Fp[ig][d];
            }
            real_t mu = 1.0 / (1.0 + xion[0]);
            real_t p_new = nH * 1.38e-16 * T2 / mu / params::units_pressure;
            grid_.uold(idc, 5) = p_new / (grid_.gamma - 1.0) + 0.5 * grid_.uold(idc, 2)*grid_.uold(idc, 2)/grid_.uold(idc, 1);
        }
        igrid = grid_.next[igrid - 1];
    }
}

void RtSolver::set_unew(int ilevel) {
    if (nGroups <= 0) return;
    int nvar_hydro = 5;
#ifdef MHD
    nvar_hydro = 11;
#endif
    int nrtvar = nGroups * (1 + NDIM) + nIons_;
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            for (int iv = 1; iv <= nrtvar; ++iv) {
                grid_.unew(idc, nvar_hydro + iv) = grid_.uold(idc, nvar_hydro + iv);
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

void RtSolver::set_uold(int ilevel) {
    if (nGroups <= 0) return;
    int nvar_hydro = 5;
#ifdef MHD
    nvar_hydro = 11;
#endif
    int nrtvar = nGroups * (1 + NDIM) + nIons_;
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            for (int iv = 1; iv <= nrtvar; ++iv) {
                grid_.uold(idc, nvar_hydro + iv) = grid_.unew(idc, nvar_hydro + iv);
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

} // namespace ramses
