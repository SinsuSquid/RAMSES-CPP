#include "ramses/solvers/mhd/MhdSolver.hpp"
#include "ramses/solvers/physics/EquationOfState.hpp"
#include "ramses/core/MpiManager.hpp"
#include "ramses/core/Parameters.hpp"
#include "ramses/core/Constants.hpp"
#include "ramses/solvers/hydro/SlopeLimiter.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

namespace ramses {
MhdSolver::MhdSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {
    stencil_ptr_ = std::make_unique<LocalStencil>();
}

MhdSolver::~MhdSolver() = default;


real_t MhdSolver::compute_max_div_b(int ilevel, real_t dx) {
    real_t max_div = 0.0;
    int nvar_pure = grid_.nvar - 3;
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= (1 << NDIM); ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                real_t div = (grid_.uold(idc, nvar_pure + 1) - grid_.uold(idc, 6));
                if (NDIM > 1) div += (grid_.uold(idc, nvar_pure + 2) - grid_.uold(idc, 7));
                if (NDIM > 2) div += (grid_.uold(idc, nvar_pure + 3) - grid_.uold(idc, 8));
                max_div = std::max(max_div, std::abs(div) / dx);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    return max_div;
}

void MhdSolver::set_unew(int ilevel) {
    int n2d = (1 << NDIM);
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                for (int ivar = 1; ivar <= grid_.nvar; ++ivar) grid_.unew(idc, ivar) = grid_.uold(idc, ivar);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void MhdSolver::set_uold(int ilevel) {
    int n2d = (1 << NDIM);
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                for (int ivar = 1; ivar <= grid_.nvar; ++ivar) grid_.uold(idc, ivar) = grid_.unew(idc, ivar);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

real_t MhdSolver::compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor) {
    real_t dt_max = params::boxlen / 1e-10, smallr = 1e-10;
    int n2d = (1 << NDIM);
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                real_t d = std::max(grid_.uold(idc, 1), smallr), u = grid_.uold(idc, 2)/d, v = grid_.uold(idc, 3)/d, w = grid_.uold(idc, 4)/d;
                real_t A = grid_.uold(idc, 6), B = grid_.uold(idc, 7), C = grid_.uold(idc, 8), etot = grid_.uold(idc, 5);
                real_t ekin = 0.5*d*(u*u+v*v+w*w), emag = 0.125*std::pow(A+grid_.uold(idc,grid_.nvar-2),2) + 0.125*std::pow(B+grid_.uold(idc,grid_.nvar-1),2) + 0.125*std::pow(C+grid_.uold(idc,grid_.nvar),2);
                real_t e_nonthermal = 0.0;
                int iener = 5;
                for (int ie = 0; ie < nener_; ++ie) {
                    e_nonthermal += grid_.uold(idc, iener + 1 + ie);
                }
                real_t p = EquationOfState::get_pressure(d, etot - ekin - emag - e_nonthermal, gamma);
                real_t q_mhd[64] = {0};
                q_mhd[0] = d; q_mhd[1] = p; q_mhd[2] = u;
                q_mhd[3] = 0.5*(A+grid_.uold(idc,grid_.nvar-2));
                q_mhd[4] = v;
                q_mhd[5] = 0.5*(B+grid_.uold(idc,grid_.nvar-1));
                q_mhd[6] = w;
                q_mhd[7] = 0.5*(C+grid_.uold(idc,grid_.nvar));
                for (int ie = 0; ie < nener_; ++ie) {
                    q_mhd[8 + ie] = grid_.uold(idc, iener + 1 + ie) * (grid_.gamma_rad[ie] - 1.0);
                }
                real_t vel_fast;
                find_speed_fast(q_mhd, vel_fast, gamma);
                dt_max = std::min(dt_max, courant_factor * dx / (std::sqrt(u*u+v*v+w*w) + vel_fast));
                for(int idim=1; idim<=NDIM; ++idim) {
                    real_t acc = std::abs(grid_.f(idc, idim));
                    if (acc > 0) dt_max = std::min(dt_max, courant_factor * std::sqrt(dx / acc));
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    return dt_max;
}

void MhdSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    set_unew(ilevel);
    std::vector<int> active_octs;
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) { active_octs.push_back(igrid); igrid = grid_.next[igrid - 1]; }
    }
    if (!active_octs.empty()) godfine1(active_octs, ilevel, dt, dx);
}

void MhdSolver::gather_stencil(int igrid, int ilevel, LocalStencil& stencil) {
    int ind_father = grid_.father[igrid - 1], nbors_27[27]; grid_.get_27_cell_neighbors(ind_father, nbors_27);
    int pos = (ind_father > grid_.ncoarse) ? ((ind_father - grid_.ncoarse - 1) / grid_.ngridmax) + 1 : 1;
    int myid = MpiManager::instance().rank() + 1, nvar = grid_.nvar;
    // For 1D/2D, the hardcoded 3D arrays lll/mmm are not applicable; use fallback logic
#if NDIM == 1
    const int nstencil = 2;  // 2 neighbors in 1D
#elif NDIM == 2
    const int nstencil = 9;  // 9 neighbors in 2D
#else
    const int nstencil = 27; // 27 neighbors in 3D
#endif
    for (int j = 0; j < nstencil; ++j) {
        int f_idx = (NDIM == 3) ? constants::lll[pos-1][j] : j+1;
        int c_pos = (NDIM == 3) ? constants::mmm[pos-1][j] : ((j % (1<<NDIM)) + 1);
        int if_n = nbors_27[f_idx - 1];
        int sz = j / 9; int sy = (j % 9) / 3; int sx = j % 3;
        for (int iv = 0; iv < 20; ++iv) stencil.uloc[sz][sy][sx][iv] = 0.0;
        if (if_n > 0 && if_n <= grid_.ncell && grid_.cpu_map[if_n-1] == myid) {
            if (grid_.son[if_n-1] > 0) {
                // If neighbor is refined, try to find child cell; otherwise use coarse cell value
                int child_grid = grid_.son[if_n-1];
                int idc = grid_.ncoarse + (c_pos - 1) * grid_.ngridmax + child_grid;
                if (idc >= 1 && idc <= grid_.ncell) {
                    for (int iv = 1; iv <= nvar; ++iv) stencil.uloc[sz][sy][sx][iv - 1] = grid_.uold(idc, iv);
                    stencil.refined[sz][sy][sx] = (grid_.son[idc-1] > 0);
                } else {
                    // Formula produced out-of-bounds; use parent cell instead
                    for (int iv = 1; iv <= nvar; ++iv) stencil.uloc[sz][sy][sx][iv - 1] = grid_.uold(if_n, iv);
                    stencil.refined[sz][sy][sx] = false;
                }
            } else {
                for (int iv = 1; iv <= nvar; ++iv) stencil.uloc[sz][sy][sx][iv - 1] = grid_.uold(if_n, iv);
                stencil.refined[sz][sy][sx] = false;
            }
        } else if (if_n <= 0) {
            for (int iv = 1; iv <= nvar; ++iv) stencil.uloc[sz][sy][sx][iv - 1] = grid_.uold(grid_.ncoarse + (c_pos - 1) * grid_.ngridmax + igrid, iv);
            int ib = -if_n; int bt = (ib > 0 && ib <= (int)grid_.bound_type.size()) ? grid_.bound_type[ib - 1] : 1;
            if (bt == 1) { if (f_idx == 11 || f_idx == 17) stencil.uloc[sz][sy][sx][1] *= -1.0; if (f_idx == 13 || f_idx == 15) stencil.uloc[sz][sy][sx][2] *= -1.0; }
        }
    }
}

void MhdSolver::ctoprim(const real_t u[64], real_t q[64], real_t bf[3][2], real_t gamma) {
    real_t d = std::max(u[0], 1e-10); q[0]=d; q[1]=u[1]/d; q[2]=u[2]/d; q[3]=u[3]/d;
    bf[0][0]=u[5]; bf[1][0]=u[6]; bf[2][0]=u[7];
    int nvp = grid_.nvar - 3; bf[0][1]=u[nvp+0]; bf[1][1]=u[nvp+1]; bf[2][1]=u[nvp+2];
    real_t ek = 0.5*d*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]), em = 0.125*(std::pow(bf[0][0]+bf[0][1],2)+std::pow(bf[1][0]+bf[1][1],2)+std::pow(bf[2][0]+bf[2][1],2));
    real_t e_nonthermal = 0.0;
    for (int ie = 0; ie < nener_; ++ie) {
        int iv = 8 + ie;
        q[iv] = u[iv] * (grid_.gamma_rad[ie] - 1.0);
        e_nonthermal += u[iv];
    }
    q[4] = EquationOfState::get_pressure(d, u[4] - ek - em - e_nonthermal, gamma);
    q[5]=0.5*(bf[0][0]+bf[0][1]); q[6]=0.5*(bf[1][0]+bf[1][1]); q[7]=0.5*(bf[2][0]+bf[2][1]);
    for (int iv = 8 + nener_; iv < nvp; ++iv) q[iv] = u[iv] / d;
}

void MhdSolver::trace(const real_t qloc[6][6][6][64], const real_t bfloc[6][6][6][3][2], const real_t dq[6][6][6][3][64], const real_t dbf[6][6][6][3][2], real_t dt, real_t dx, real_t qm[6][6][6][3][64], real_t qp[6][6][6][3][64]) {
    real_t dtdx = dt/dx;
    for (int k = 1; k < 5; ++k) for (int j = 1; j < 5; ++j) for (int i = 1; i < 5; ++i) {
        real_t r=qloc[i][j][k][0], p=qloc[i][j][k][4], u=qloc[i][j][k][1], v=qloc[i][j][k][2], w=qloc[i][j][k][3];
        for(int id=0; id<3; ++id) {
            real_t dr=dq[i][j][k][id][0], dp=dq[i][j][k][id][4], du=dq[i][j][k][id][1];
            real_t sr = -0.5*dtdx*(u*dr + r*du), su = -0.5*dtdx*(u*du + dp/r);
            for(int iv=0; iv<grid_.nvar; ++iv) { 
                real_t s = (iv==0)?sr:(iv==1)?su:0.0;
                qp[i][j][k][id][iv] = qloc[i][j][k][iv] + s - 0.5*dq[i][j][k][id][iv]; 
                qm[i][j][k][id][iv] = qloc[i][j][k][iv] + s + 0.5*dq[i][j][k][id][iv]; 
            }
        }
    }
}

void MhdSolver::cmpflxm(const real_t qm[6][6][6][3][64], const real_t qp[6][6][6][3][64], const real_t bfloc[6][6][6][3][2], int idim, real_t gamma, real_t flux[6][6][6][64]) {
    for (int k=1; k<5; ++k) for (int j=1; j<5; ++j) for (int i=1; i<5; ++i) {
        int ni=i+(idim==0), nj=j+(idim==1), nk=k+(idim==2);
        real_t qL[64] = {0}, qR[64] = {0}, f[9];
        auto fl = [&](real_t* q, const real_t s[64], real_t bn) { q[0]=s[0]; q[1]=s[4]; q[2]=s[1+idim]; q[3]=bn; q[4]=s[1+(idim+1)%3]; q[5]=s[5+(idim+1)%3]; q[6]=s[1+(idim+2)%3]; q[7]=s[5+(idim+2)%3]; };
        fl(qL, qm[i][j][k][idim], bfloc[i][j][k][idim][1]); fl(qR, qp[ni][nj][nk][idim], bfloc[ni][nj][nk][idim][0]);
        for (int ie = 0; ie < nener_; ++ie) {
            qL[8 + ie] = qm[i][j][k][idim][8 + ie];
            qR[8 + ie] = qp[ni][nj][nk][idim][8 + ie];
        }
        if (config_.get("hydro_params", "riemann", "hlld") == "llf") llf(qL, qR, f, gamma); else hlld(qL, qR, f, gamma);
        flux[i][j][k][0]=f[0]; flux[i][j][k][1]=f[1]; flux[i][j][k][2]=f[2]; flux[i][j][k][4]=f[4]; flux[i][j][k][6]=f[6]; flux[i][j][k][5]=f[5]; flux[i][j][k][7]=f[7];
        
        real_t mass_flux = f[0];
        real_t u_interface = mass_flux / std::max((mass_flux > 0) ? qL[0] : qR[0], (real_t)1e-10);
        for (int ie = 0; ie < nener_; ++ie) {
            int iv = 8 + ie;
            real_t p_rad = (mass_flux > 0) ? qm[i][j][k][idim][iv] : qp[ni][nj][nk][idim][iv];
            flux[i][j][k][iv] = u_interface * (p_rad / (grid_.gamma_rad[ie] - 1.0));
        }
        int nvp = grid_.nvar - 3;
        for (int iv = 8 + nener_; iv < nvp; ++iv) {
            real_t concentration = (mass_flux > 0) ? qm[i][j][k][idim][iv] : qp[ni][nj][nk][idim][iv];
            flux[i][j][k][iv] = mass_flux * concentration;
        }
    }
}

void MhdSolver::godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx) {
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4), dt_dx = dt/dx;
    for (int igrid : ind_grid) {
        gather_stencil(igrid, ilevel, *stencil_ptr_);
        real_t qloc[6][6][6][64], bfloc[6][6][6][3][2], dq[6][6][6][3][64]={0}, dbf[6][6][6][3][2]={0};
        for(int k=0; k<6; ++k) for(int j=0; j<6; ++j) for(int i=0; i<6; ++i) ctoprim(stencil_ptr_->uloc[i][j][k], qloc[i][j][k], bfloc[i][j][k], gamma);
        int st = config_.get_int("hydro_params", "slope_type", 1);
        for(int k=1; k<5; ++k) for(int j=1; j<5; ++j) for(int i=1; i<5; ++i) {
            for(int iv=0; iv<grid_.nvar; ++iv) { dq[i][j][k][0][iv]=SlopeLimiter::compute_slope(qloc[i-1][j][k][iv],qloc[i][j][k][iv],qloc[i+1][j][k][iv],st); if(NDIM>1) dq[i][j][k][1][iv]=SlopeLimiter::compute_slope(qloc[i][j-1][k][iv],qloc[i][j][k][iv],qloc[i][j+1][k][iv],st); }
        }
        real_t qm[6][6][6][3][64], qp[6][6][6][3][64], fluxes[3][6][6][6][64]={0}, emfz[6][6][6]={0};
        trace(qloc, bfloc, dq, dbf, dt, dx, qm, qp);
        for(int d=0; d<NDIM; ++d) cmpflxm(qm, qp, bfloc, d, gamma, fluxes[d]);
        for(int k=1; k<6; ++k) for(int j=1; j<6; ++j) for(int i=1; i<6; ++i) if(NDIM>1) emfz[i][j][k]=0.25*(fluxes[0][i][j][k][6]+fluxes[0][i][j-1][k][6]-fluxes[1][i][j][k][5]-fluxes[1][i-1][j][k][5]);
        for(int k2=0; k2<(NDIM>2?2:1); ++k2) for(int j2=0; j2<(NDIM>1?2:1); ++j2) for(int i2=0; i2<2; ++i2) {
            int i=2+i2, j=2+j2, k=2+k2, icp=1+i2+2*j2+4*k2, idc=grid_.ncoarse+(icp-1)*grid_.ngridmax+igrid;
            if(icp<=(1<<NDIM)) {
                for(int d=0; d<NDIM; ++d) {
                    int il=i-(d==0), jl=j-(d==1), kl=k-(d==2);
                    real_t fL[64], fR[64]; for(int iv=0; iv<grid_.nvar; ++iv) { fL[iv]=fluxes[d][il][jl][kl][iv]; fR[iv]=fluxes[d][i][j][k][iv]; }
                    grid_.unew(idc, 1) += (fL[0]-fR[0])*dt_dx; grid_.unew(idc, 5) += (fL[1]-fR[1])*dt_dx;
                    grid_.unew(idc, 1+d+1) += (fL[2]-fR[2])*dt_dx; grid_.unew(idc, 1+(d+1)%3+1) += (fL[4]-fR[4])*dt_dx; grid_.unew(idc, 1+(d+2)%3+1) += (fL[6]-fR[6])*dt_dx;
                    for (int ie = 0; ie < nener_; ++ie) {
                        grid_.unew(idc, 9 + ie) += (fL[8 + ie] - fR[8 + ie]) * dt_dx;
                    }
                    int nvp = grid_.nvar - 3;
                    for (int iv = 8 + nener_; iv < nvp; ++iv) {
                        grid_.unew(idc, 1 + iv) += (fL[iv] - fR[iv]) * dt_dx;
                    }
                }
                if (nener_ > 0) {
                    int ign[7]; grid_.get_nbor_grids(igrid, ign);
                    int icn[6]; grid_.get_nbor_cells(ign, icp, icn, igrid);
                    real_t divu = 0.0;
                    real_t d_c = std::max(grid_.uold(idc, 1), (real_t)1e-10);
                    for (int idim = 0; idim < NDIM; ++idim) {
                        int id_l = icn[idim * 2];
                        int id_r = icn[idim * 2 + 1];
                        
                        real_t v_l = 0.0;
                        real_t dx_l = dx;
                        if (id_l > 0) {
                            real_t d_l = std::max(grid_.uold(id_l, 1), (real_t)1e-10);
                            v_l = grid_.uold(id_l, 2 + idim) / d_l;
                        } else {
                            v_l = grid_.uold(idc, 2 + idim) / d_c;
                            int ib = -id_l;
                            if (ib > 0 && ib <= (int)grid_.bound_type.size() && grid_.bound_type[ib - 1] == 1) {
                                v_l *= -1.0;
                            }
                            dx_l = dx * 1.5;
                        }
                        
                        real_t v_r = 0.0;
                        real_t dx_r = dx;
                        if (id_r > 0) {
                            real_t d_r = std::max(grid_.uold(id_r, 1), (real_t)1e-10);
                            v_r = grid_.uold(id_r, 2 + idim) / d_r;
                        } else {
                            v_r = grid_.uold(idc, 2 + idim) / d_c;
                            int ib = -id_r;
                            if (ib > 0 && ib <= (int)grid_.bound_type.size() && grid_.bound_type[ib - 1] == 1) {
                                v_r *= -1.0;
                            }
                            dx_r = dx * 1.5;
                        }
                        
                        divu += (v_r - v_l) / (dx_l + dx_r);
                    }
                    for (int ie = 0; ie < nener_; ++ie) {
                        grid_.unew(idc, 9 + ie) -= (grid_.gamma_rad[ie] - 1.0) * grid_.uold(idc, 9 + ie) * divu * dt;
                    }
                }
                if(NDIM>1){
                    grid_.unew(idc, 6) += (emfz[i-1][j][k]-emfz[i-1][j+1][k])*dt_dx; grid_.unew(idc, 7) += (emfz[i][j-1][k]-emfz[i+1][j-1][k])*dt_dx;
                    grid_.unew(idc, grid_.nvar - 2) += (emfz[i][j][k]-emfz[i][j+1][k])*dt_dx; grid_.unew(idc, grid_.nvar - 1) += (emfz[i][j][k]-emfz[i+1][j][k])*dt_dx;
                }
            }
        }
        // MHD Refluxing (EMF Averaging at Coarse-Fine interfaces)
        if(ilevel > 1 && NDIM > 1) {
            int ind_father_cell = grid_.father[igrid-1];   // Returns cell index (1..ncell)
            // Find parent grid: cell index = ncoarse + (ic-1)*ngridmax + igrid_parent for each ic
            // Try all ic values to extract the matching parent grid index
            for(int ic=1; ic<=constants::twotondim; ic++) {
                int igrid_parent = ind_father_cell - grid_.ncoarse - (ic-1)*grid_.ngridmax;
                if(igrid_parent >= 1 && igrid_parent <= grid_.ngridmax) {
                    int idc = grid_.ncoarse+(ic-1)*grid_.ngridmax+igrid_parent;
                    if(idc >= 1 && idc <= grid_.ncell && grid_.son[idc-1] == igrid) {
                        // Found the correct parent grid
                        // TODO: Implement full EMF refluxing to parent grid boundaries
                        // For now, just break after identifying the parent grid
                        break;
                    }
                }
            }
        }
    }
}

void MhdSolver::llf(const real_t* ql, const real_t* qr, real_t* f, real_t gamma) {
    real_t fl[9], fr[9];
    real_t cl, cr;
    real_t et = 1.0 / (gamma - 1.0);

    // Compute MHD fluxes and fast wave speeds
    find_mhd_flux(ql, nullptr, fl, gamma);
    find_mhd_flux(qr, nullptr, fr, gamma);
    find_speed_fast(ql, cl, gamma);
    find_speed_fast(qr, cr, gamma);

    // Maximum signal speed
    real_t sm = std::max(std::abs(ql[2]) + cl, std::abs(qr[2]) + cr);

    // Helper to calculate conserved variable state from primitive state
    auto primitive_to_conserved = [&](const real_t* q, real_t* c) {
        c[0] = q[0];                              // density (rho)
        c[1] = q[0] * q[2];                       // momentum x (rho*u)
        c[2] = q[0] * q[4];                       // momentum y (rho*v)
        c[3] = q[0] * q[6];                       // momentum z (rho*w)
        
        // Total energy: internal + kinetic + magnetic
        c[4] = q[1] * et 
             + 0.5 * q[0] * (q[2]*q[2] + q[4]*q[4] + q[6]*q[6]) 
             + 0.5 * (q[3]*q[3] + q[5]*q[5] + q[7]*q[7]);
             
        c[5] = q[3];                              // magnetic field component x (Bx)
        c[6] = q[5];                              // magnetic field component y (By)
        c[7] = q[7];                              // magnetic field component z (Bz)
        c[8] = q[1] * et;                         // internal energy
    };

    real_t cll[9], crr[9];
    primitive_to_conserved(ql, cll);
    primitive_to_conserved(qr, crr);

    // Local Lax-Friedrichs flux formula
    for (int i = 0; i < 9; ++i) {
        f[i] = 0.5 * (fl[i] + fr[i] - sm * (crr[i] - cll[i]));
    }
}

void MhdSolver::hlld(const real_t* ql_in, const real_t* qr_in, real_t* fgdnv, real_t gamma) {
    real_t entho = 1.0 / (gamma - 1.0);
    real_t half = 0.5;
    real_t A = half * (ql_in[3] + qr_in[3]);
    real_t sgnm = (A >= 0) ? 1.0 : -1.0;

    real_t ql[8], qr[8];
    for (int i = 0; i < 8; ++i) {
        ql[i] = ql_in[i];
        qr[i] = qr_in[i];
    }
    ql[3] = A;
    qr[3] = A;

    // Extract primitive variables for left state
    real_t rl = ql[0], Pl = ql[1], ul = ql[2], vl = ql[4], Bl = ql[5], wl = ql[6], Cl = ql[7];
    // Extract primitive variables for right state
    real_t rr = qr[0], Pr = qr[1], ur = qr[2], vr = qr[4], Br = qr[5], wr = qr[6], Cr = qr[7];

    // Compute left state kinetic, magnetic, total energy, total pressure, and B-dot-v
    real_t el = half * rl * (ul * ul + vl * vl + wl * wl);
    real_t emagl = half * (A * A + Bl * Bl + Cl * Cl);
    real_t etl = Pl * entho + el + emagl;
    real_t Ptl = Pl + emagl;
    real_t vBl = ul * A + vl * Bl + wl * Cl;

    // Compute right state kinetic, magnetic, total energy, total pressure, and B-dot-v
    real_t er = half * rr * (ur * ur + vr * vr + wr * wr);
    real_t emagr = half * (A * A + Br * Br + Cr * Cr);
    real_t etr = Pr * entho + er + emagr;
    real_t Ptr = Pr + emagr;
    real_t vBr = ur * A + vr * Br + wr * Cr;

    // Wave speeds
    real_t cl, cr;
    find_speed_fast(ql_in, cl, gamma);
    find_speed_fast(qr_in, cr, gamma);
    real_t SL = std::min(ul, ur) - std::max(cl, cr);
    real_t SR = std::max(ul, ur) + std::max(cl, cr);

    // Compute star state velocity and total pressure
    real_t rcl = rl * (ul - SL);
    real_t rcr = rr * (SR - ur);
    real_t ust = (rcr * ur + rcl * ul + (Ptl - Ptr)) / (rcr + rcl);
    real_t Ptst = (rcr * Ptl + rcl * Ptr + rcl * rcr * (ul - ur)) / (rcr + rcl);

    // Compute left star state density and field components
    real_t rstl = rl * (SL - ul) / (SL - ust);
    real_t estl = rl * (SL - ul) * (SL - ust) - A * A;
    real_t ell = rl * (SL - ul) * (SL - ul) - A * A;
    real_t vstl, Bstl, wstl, Cstl;
    if (std::abs(estl) < 1e-4 * A * A) {
        vstl = vl;
        Bstl = Bl;
        wstl = wl;
        Cstl = Cl;
    } else {
        vstl = vl - A * Bl * (ust - ul) / estl;
        Bstl = Bl * ell / estl;
        wstl = wl - A * Cl * (ust - ul) / estl;
        Cstl = Cl * ell / estl;
    }

    // Left star state energies and Alfven speed wave checks
    real_t vdbstl = ust * A + vstl * Bstl + wstl * Cstl;
    real_t etstl = ((SL - ul) * etl - Ptl * ul + Ptst * ust + A * (vBl - vdbstl)) / (SL - ust);
    real_t sqrl = std::sqrt(rstl);
    real_t calfl = std::abs(A) / sqrl;
    real_t SAL = ust - calfl;

    // Compute right star state density and field components
    real_t rstr = rr * (SR - ur) / (SR - ust);
    real_t estr = rr * (SR - ur) * (SR - ust) - A * A;
    real_t err = rr * (SR - ur) * (SR - ur) - A * A;
    real_t vstr, Bstr, wstr, Cstr;
    if (std::abs(estr) < 1e-4 * A * A) {
        vstr = vr;
        Bstr = Br;
        wstr = wr;
        Cstr = Cr;
    } else {
        vstr = vr - A * Br * (ust - ur) / estr;
        Bstr = Br * err / estr;
        wstr = wr - A * Cr * (ust - ur) / estr;
        Cstr = Cr * err / estr;
    }

    // Right star state energies and Alfven speed wave checks
    real_t vdbstr = ust * A + vstr * Bstr + wstr * Cstr;
    real_t etstr = ((SR - ur) * etr - Ptr * ur + Ptst * ust + A * (vBr - vdbstr)) / (SR - ust);
    real_t sqrr = std::sqrt(rstr);
    real_t calfr = std::abs(A) / sqrr;
    real_t SAR = ust + calfr;

    // Double-star state velocities and magnetic fields
    real_t vss = (sqrl * vstl + sqrr * vstr + sgnm * (Bstr - Bstl)) / (sqrl + sqrr);
    real_t wss = (sqrl * wstl + sqrr * wstr + sgnm * (Cstr - Cstl)) / (sqrl + sqrr);
    real_t Bss = (sqrl * Bstl + sqrr * Bstr + sgnm * sqrl * sqrr * (vstr - vstl)) / (sqrl + sqrr);
    real_t Css = (sqrl * Cstl + sqrr * Cstr + sgnm * sqrl * sqrr * (wstr - wstl)) / (sqrl + sqrr);
    real_t vdbss = ust * A + vss * Bss + wss * Css;
    real_t etssl = etstl - sgnm * sqrl * (vdbstl - vdbss);
    real_t etssr = etstr + sgnm * sqrr * (vdbstr - vdbss);

    // Resolve final state based on SL, SAL, ust, SAR, SR wave speeds
    real_t ro, uo, vo, wo, Bo, Co, Ptoto, etoto, vdb, ei;
    real_t eintl = Pl * entho;
    real_t eintr = Pr * entho;

    if (SL > 0) {
        ro = rl; uo = ul; vo = vl; wo = wl; Bo = Bl; Co = Cl; Ptoto = Ptl; etoto = etl; vdb = vBl; ei = eintl;
    } else if (SAL > 0) {
        ro = rstl; uo = ust; vo = vstl; wo = wstl; Bo = Bstl; Co = Cstl; Ptoto = Ptst; etoto = etstl; vdb = vdbstl; ei = eintl * (SL - ul) / (SL - ust);
    } else if (ust > 0) {
        ro = rstl; uo = ust; vo = vss; wo = wss; Bo = Bss; Co = Css; Ptoto = Ptst; etoto = etssl; vdb = vdbss; ei = eintl * (SL - ul) / (SL - ust);
    } else if (SAR > 0) {
        ro = rstr; uo = ust; vo = vss; wo = wss; Bo = Bss; Co = Css; Ptoto = Ptst; etoto = etssr; vdb = vdbss; ei = eintr * (SR - ur) / (SR - ust);
    } else if (SR > 0) {
        ro = rstr; uo = ust; vo = vstr; wo = wstr; Bo = Bstr; Co = Cstr; Ptoto = Ptst; etoto = etstr; vdb = vdbstr; ei = eintr * (SR - ur) / (SR - ust);
    } else {
        ro = rr; uo = ur; vo = vr; wo = wr; Bo = Br; Co = Cr; Ptoto = Ptr; etoto = etr; vdb = vBr; ei = eintr;
    }

    // Set Godunov flux
    fgdnv[0] = ro * uo;
    fgdnv[1] = (etoto + Ptoto) * uo - A * vdb;
    fgdnv[2] = ro * uo * uo + Ptoto - A * A;
    fgdnv[3] = 0;
    fgdnv[4] = ro * uo * vo - A * Bo;
    fgdnv[5] = Bo * uo - A * vo;
    fgdnv[6] = ro * uo * wo - A * Co;
    fgdnv[7] = Co * uo - A * wo;
    fgdnv[8] = uo * ei;
}

void MhdSolver::find_mhd_flux(const real_t* q, real_t* c, real_t* f, real_t gamma) {
    real_t et = 1.0 / (gamma - 1.0);
    real_t d = q[0], p = q[1], u = q[2], A = q[3], v = q[4], B = q[5], w = q[6], C = q[7];
    real_t ek = 0.5 * d * (u * u + v * v + w * w);
    real_t em = 0.5 * (A * A + B * B + C * C);
    real_t etot = p * et + ek + em;
    real_t pt = p + em;
    
    f[0] = d * u;
    f[1] = (etot + pt) * u - A * (A * u + B * v + C * w);
    f[2] = d * u * u + pt - A * A;
    f[3] = 0;
    f[4] = d * u * v - A * B;
    f[5] = B * u - A * v;
    f[6] = d * u * w - A * C;
    f[7] = C * u - A * w;
    f[8] = p * et * u;
}

void MhdSolver::find_speed_fast(const real_t* q, real_t& v, real_t gamma) {
    real_t d = q[0], p = q[1], A = q[3], B = q[5], C = q[7];
    real_t b2 = A * A + B * B + C * C;
    real_t c2 = gamma * p / d;
    for (int ie = 0; ie < nener_; ++ie) {
        c2 += grid_.gamma_rad[ie] * q[8 + ie] / d;
    }
    real_t d2 = 0.5 * (b2 / d + c2);
    v = std::sqrt(d2 + std::sqrt(std::max(0.0, d2 * d2 - c2 * A * A / d)));
}

void MhdSolver::get_diagnostics(int il, real_t dx, real_t& mi, real_t& mv, real_t& mb) {}

} // namespace ramses
