#include "ramses/MhdSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Muscl.hpp"
#include "ramses/SlopeLimiter.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

namespace ramses {

real_t MhdSolver::compute_max_div_b(int ilevel, real_t dx) {
    real_t max_div = 0.0;
    int nvar_pure = grid_.nvar - 3;
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
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
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                for (int ivar = 1; ivar <= grid_.nvar; ++ivar) grid_.unew(idc, ivar) = grid_.uold(idc, ivar);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void MhdSolver::set_uold(int ilevel) {
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                for (int ivar = 1; ivar <= grid_.nvar; ++ivar) grid_.uold(idc, ivar) = grid_.unew(idc, ivar);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

real_t MhdSolver::compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor) {
    real_t dt_max = 1e30, smallr = 1e-10;
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                real_t d = std::max(grid_.uold(idc, 1), smallr), u = grid_.uold(idc, 2)/d, v = grid_.uold(idc, 3)/d, w = grid_.uold(idc, 4)/d;
                real_t A = grid_.uold(idc, 6), B = grid_.uold(idc, 7), C = grid_.uold(idc, 8), etot = grid_.uold(idc, 5);
                real_t ekin = 0.5*d*(u*u+v*v+w*w), emag = 0.5*(A*A+B*B+C*C), p = std::max((etot-ekin-emag)*(gamma-1.0), d*1e-10);
                real_t q_mhd[8] = {d, p, u, A, v, B, w, C}, vel_fast;
                find_speed_fast(q_mhd, vel_fast, gamma);
                dt_max = std::min(dt_max, dx / (std::sqrt(u*u+v*v+w*w) + vel_fast));
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    return dt_max * courant_factor;
}

void MhdSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    std::vector<int> active_octs;
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) { active_octs.push_back(igrid); igrid = grid_.next[igrid - 1]; }
    }
    if (!active_octs.empty()) godfine1(active_octs, ilevel, dt, dx);
}

void MhdSolver::gather_stencil(int igrid, int ilevel, LocalStencil& stencil) {
    int nbors_father[27]; grid_.get_3x3x3_father(igrid, nbors_father);
    int nvar = grid_.nvar;
    for (int k1 = 0; k1 < 3; ++k1) for (int j1 = 0; j1 < 3; ++j1) for (int i1 = 0; i1 < 3; ++i1) {
        int ifather = nbors_father[i1 + 3*j1 + 9*k1];
        int ison = (ifather > 0 && ifather <= (int)grid_.son.size() - 1) ? grid_.son[ifather] : 0;
        for (int k2 = 0; k2 < 2; ++k2) for (int j2 = 0; j2 < 2; ++j2) for (int i2 = 0; i2 < 2; ++i2) {
            int i3 = i1*2+i2, j3 = j1*2+j2, k3 = k1*2+k2, icp = 1+i2+2*j2+4*k2;
            for (int iv = 0; iv < 20; ++iv) stencil.uloc[i3][j3][k3][iv] = 0.0;
            if (ison > 0 && icp <= constants::twotondim) {
                int idc = grid_.ncoarse + (icp-1)*grid_.ngridmax + ison;
                if (idc > 0 && idc <= grid_.ncell) {
                    for (int iv = 1; iv <= nvar; ++iv) stencil.uloc[i3][j3][k3][iv-1] = grid_.uold(idc, iv);
                    stencil.refined[i3][j3][k3] = (grid_.son[idc] > 0);
                }
            } else if (ifather > 0 && ifather <= grid_.ncell) {
                for (int iv = 1; iv <= nvar; ++iv) stencil.uloc[i3][j3][k3][iv-1] = grid_.uold(ifather, iv);
                stencil.refined[i3][j3][k3] = false;
            } else {
                int ai2 = i2, aj2 = (NDIM > 1) ? j2 : 0, ak2 = (NDIM > 2) ? k2 : 0, aicp = 1+ai2+2*aj2+4*ak2;
                if (ison > 0 && aicp <= constants::twotondim) {
                    int idc = grid_.ncoarse + (aicp-1)*grid_.ngridmax + ison;
                    if (idc > 0 && idc <= grid_.ncell) for (int iv = 1; iv <= nvar; ++iv) stencil.uloc[i3][j3][k3][iv-1] = grid_.uold(idc, iv);
                }
            }
        }
    }
}

void MhdSolver::ctoprim(const real_t u[20], real_t q[20], real_t bf[3][2], real_t gamma) {
    real_t d = std::max(u[0], 1e-10); q[0]=d; q[1]=u[1]/d; q[2]=u[2]/d; q[3]=u[3]/d;
    bf[0][0]=u[5]; bf[1][0]=u[6]; bf[2][0]=u[7];
    int nvp = grid_.nvar - 3; bf[0][1]=u[nvp+0]; bf[1][1]=u[nvp+1]; bf[2][1]=u[nvp+2];
    q[5]=0.5*(bf[0][0]+bf[0][1]); q[6]=0.5*(bf[1][0]+bf[1][1]); q[7]=0.5*(bf[2][0]+bf[2][1]);
    real_t ek=0.5*d*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]), em=0.5*(q[5]*q[5]+q[6]*q[6]+q[7]*q[7]);
    q[4]=std::max((u[4]-ek-em)*(gamma-1.0), d*1e-10);
    for (int iv=8; iv<nvp; ++iv) q[iv]=u[iv]/d;
}

void MhdSolver::trace(const real_t qloc[6][6][6][20], const real_t bfloc[6][6][6][3][2], const real_t dq[6][6][6][3][20], const real_t dbf[6][6][6][3][2], real_t dt, real_t dx, real_t qm[6][6][6][3][20], real_t qp[6][6][6][3][20]) {
    real_t dtdx = dt/dx, gamma = config_.get_double("hydro_params", "gamma", 1.4), sp = 1e-10;
    real_t ez[6][6][6]={0}, ey[6][6][6]={0}, ex[6][6][6]={0};
    for(int k=1; k<6; ++k) for(int j=1; j<6; ++j) for(int i=1; i<6; ++i) {
        real_t u=0.25*(qloc[i-1][j-1][k][1]+qloc[i-1][j][k][1]+qloc[i][j-1][k][1]+qloc[i][j][k][1]), v=0.25*(qloc[i-1][j-1][k][2]+qloc[i-1][j][k][2]+qloc[i][j-1][k][2]+qloc[i][j][k][2]);
        real_t Bx=0.5*(bfloc[i][j-1][k][0][0]+bfloc[i][j][k][0][0]), By=0.5*(bfloc[i-1][j][k][1][0]+bfloc[i][j][k][1][0]);
        ez[i][j][k]=u*By-v*Bx;
    }
    for (int k = 1; k < 5; ++k) for (int j = 1; j < 5; ++j) for (int i = 1; i < 5; ++i) {
        real_t r=qloc[i][j][k][0], u=qloc[i][j][k][1], v=qloc[i][j][k][2], w=qloc[i][j][k][3];
        real_t p=qloc[i][j][k][4], A=qloc[i][j][k][5], B=qloc[i][j][k][6], C=qloc[i][j][k][7];

        real_t drx=0.5*dq[i][j][k][0][0], dry=0.5*dq[i][j][k][1][0], drz=0.5*dq[i][j][k][2][0];
        real_t dux=0.5*dq[i][j][k][0][1], duy=0.5*dq[i][j][k][1][1], duz=0.5*dq[i][j][k][2][1];
        real_t dvx=0.5*dq[i][j][k][0][2], dvy=0.5*dq[i][j][k][1][2], dvz=0.5*dq[i][j][k][2][2];
        real_t dwx=0.5*dq[i][j][k][0][3], dwy=0.5*dq[i][j][k][1][3], dwz=0.5*dq[i][j][k][2][3];
        real_t dpx=0.5*dq[i][j][k][0][4], dpy=0.5*dq[i][j][k][1][4], dpz=0.5*dq[i][j][k][2][4];
        real_t dAx=0.5*dq[i][j][k][0][5], dAy=0.5*dq[i][j][k][1][5], dAz=0.5*dq[i][j][k][2][5];
        real_t dBx=0.5*dq[i][j][k][0][6], dBy=0.5*dq[i][j][k][1][6], dBz=0.5*dq[i][j][k][2][6];
        real_t dCx=0.5*dq[i][j][k][0][7], dCy=0.5*dq[i][j][k][1][7], dCz=0.5*dq[i][j][k][2][7];

        // 1D time prediction source terms
        real_t sr = (-u*drx - r*dux) * dtdx;
        real_t su = (-u*dux - (dpx + B*dBx + C*dCx)/r) * dtdx;
        real_t sv = (-u*dvx + A*dBx/r) * dtdx;
        real_t sw = (-u*dwx + A*dCx/r) * dtdx;
        real_t sp_pred = (-u*dpx - gamma*p*dux) * dtdx;
        real_t sB = (-u*dBy + A*dvy - B*dux) * dtdx;
        real_t sC = (-u*dCx + A*dwx - C*dux) * dtdx;

        if (NDIM > 1) {
            sr += (-v*dry - r*dvy) * dtdx;
            su += (-v*duy + B*dAy/r) * dtdx;
            sv += (-v*dvy - (dpy + A*dAy + C*dCy)/r) * dtdx;
            sw += (-v*dwy + B*dCy/r) * dtdx;
            sp_pred += (-v*dpy - gamma*p*dvy) * dtdx;
        }

        // Face-centered B prediction using EMFs
        real_t BLx = bfloc[i][j][k][0][0], BRx = bfloc[i+1][j][k][0][0];
        real_t BLy = bfloc[i][j][k][1][0], BRy = bfloc[i][j+1][k][1][0];
        real_t BLz = bfloc[i][j][k][2][0], BRz = bfloc[i][j][k+1][2][0];

        if (NDIM > 1) {
            BLx += (ez[i][j+1][k] - ez[i][j][k]) * dtdx * 0.5;
            BRx += (ez[i+1][j+1][k] - ez[i+1][j][k]) * dtdx * 0.5;
            BLy -= (ez[i+1][j][k] - ez[i][j][k]) * dtdx * 0.5;
            BRy -= (ez[i+1][j+1][k] - ez[i+1][j+1][k]) * dtdx * 0.5;
        }
        
        auto set_st = [&](int id, real_t dr, real_t du, real_t dv, real_t dw, real_t dp, real_t dA, real_t dB, real_t dC, real_t BL, real_t BR) {
            real_t r_pred = r + sr, u_pred = u + su, v_pred = v + sv, w_pred = w + sw, p_pred = p + sp_pred;
            qp[i][j][k][id][0]=std::max(1e-10, r_pred - dr); qm[i][j][k][id][0]=std::max(1e-10, r_pred + dr);
            qp[i][j][k][id][1]=u_pred - du; qm[i][j][k][id][1]=u_pred + du;
            qp[i][j][k][id][2]=v_pred - dv; qm[i][j][k][id][2]=v_pred + dv;
            qp[i][j][k][id][3]=w_pred - dw; qm[i][j][k][id][3]=w_pred + dw;
            qp[i][j][k][id][4]=std::max(1e-10, p_pred - dp); qm[i][j][k][id][4]=std::max(1e-10, p_pred + dp);
            
            if(id==0){ qp[i][j][k][id][5]=BL; qm[i][j][k][id][5]=BR; qp[i][j][k][id][6]=B+sB-dB; qm[i][j][k][id][6]=B+sB+dB; qp[i][j][k][id][7]=C+sC-dC; qm[i][j][k][id][7]=C+sC+dC; }
            else if(id==1){ qp[i][j][k][id][5]=A+sB-dA; qm[i][j][k][id][5]=A+sB+dA; qp[i][j][k][id][6]=BL; qm[i][j][k][id][6]=BR; qp[i][j][k][id][7]=C+sC-dC; qm[i][j][k][id][7]=C+sC+dC; }
            else { qp[i][j][k][id][5]=A-dA; qm[i][j][k][id][5]=A+dA; qp[i][j][k][id][6]=B-dB; qm[i][j][k][id][6]=B+dB; qp[i][j][k][id][7]=BL; qm[i][j][k][id][7]=BR; }
        };
        set_st(0, drx, dux, dvx, dwx, dpx, dAx, dBx, dCx, BLx, BRx);
        if (NDIM > 1) set_st(1, dry, duy, dvy, dwy, dpy, dAy, dBy, dCy, BLy, BRy);
        if (NDIM > 2) set_st(2, drz, duz, dvz, dwz, dpz, dAz, dBz, dCz, BLz, BRz);
    }
}


void MhdSolver::cmpflxm(const real_t qm[6][6][6][3][20], const real_t qp[6][6][6][3][20], const real_t bfloc[6][6][6][3][2], int idim, real_t gamma, real_t flux[6][6][6][20]) {
    int nvp = grid_.nvar - 3;
    for (int k=1; k<5; ++k) for (int j=1; j<5; ++j) for (int i=1; i<5; ++i) {
        int ni=i+(idim==0), nj=j+(idim==1), nk=k+(idim==2);
        real_t qL[8], qR[8], f[9];
        auto fl = [&](real_t* q, const real_t s[20], real_t bn) { q[0]=s[0]; q[1]=s[4]; q[2]=s[1+idim]; q[3]=bn; q[4]=s[1+(idim+1)%3]; q[5]=s[5+(idim+1)%3]; q[6]=s[1+(idim+2)%3]; q[7]=s[5+(idim+2)%3]; };
        fl(qL, qm[i][j][k][idim], bfloc[i][j][k][idim][1]); fl(qR, qp[ni][nj][nk][idim], bfloc[ni][nj][nk][idim][0]);
        if (config_.get("hydro_params", "riemann", "hlld") == "llf") llf(qL, qR, f, gamma); else hlld(qL, qR, f, gamma);
        flux[i][j][k][0]=f[0]; flux[i][j][k][4]=f[1]; flux[i][j][k][1+idim]=f[2]; flux[i][j][k][1+(idim+1)%3]=f[4]; flux[i][j][k][5+(idim+1)%3]=f[5]; flux[i][j][k][1+(idim+2)%3]=f[6]; flux[i][j][k][5+(idim+2)%3]=f[7];
        for(int iv=8; iv<nvp; ++iv) flux[i][j][k][iv] = (f[0]>0) ? f[0]*qm[i][j][k][idim][iv] : f[0]*qp[ni][nj][nk][idim][iv];
    }
}

void MhdSolver::godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx) {
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4), dt_dx = dt/dx;
    int nvp = grid_.nvar - 3;
    for (int igrid : ind_grid) {
        gather_stencil(igrid, ilevel, *stencil_ptr_);
        real_t qloc[6][6][6][20], bfloc[6][6][6][3][2], dq[6][6][6][3][20]={0}, dbf[6][6][6][3][2]={0};
        for(int k=0; k<6; ++k) for(int j=0; j<6; ++j) for(int i=0; i<6; ++i) ctoprim(stencil_ptr_->uloc[i][j][k], qloc[i][j][k], bfloc[i][j][k], gamma);
        int st = config_.get_int("hydro_params", "slope_type", 1);
        for(int k=1; k<5; ++k) for(int j=1; j<5; ++j) for(int i=1; i<5; ++i) {
            for(int iv=0; iv<nvp; ++iv) { dq[i][j][k][0][iv]=SlopeLimiter::compute_slope(qloc[i-1][j][k][iv],qloc[i][j][k][iv],qloc[i+1][j][k][iv],st); if(NDIM>1) dq[i][j][k][1][iv]=SlopeLimiter::compute_slope(qloc[i][j-1][k][iv],qloc[i][j][k][iv],qloc[i][j+1][k][iv],st); }
            for(int f=i; f<=i+1; ++f) if(NDIM>1) dbf[f][j][k][0][0]=SlopeLimiter::compute_slope(bfloc[f][j-1][k][0][0],bfloc[f][j][k][0][0],bfloc[f][j+1][k][0][0],st);
            for(int f=j; f<=j+1; ++f) if(NDIM>1) dbf[i][f][k][1][0]=SlopeLimiter::compute_slope(bfloc[i-1][f][k][1][0],bfloc[i][f][k][1][0],bfloc[i+1][f][k][1][0],st);
        }
        real_t qm[6][6][6][3][20], qp[6][6][6][3][20], fluxes[3][6][6][6][20]={0}, emfz[6][6][6]={0};
        trace(qloc, bfloc, dq, dbf, dt, dx, qm, qp);
        for(int d=0; d<NDIM; ++d) cmpflxm(qm, qp, bfloc, d, gamma, fluxes[d]);
        for(int k=1; k<6; ++k) for(int j=1; j<6; ++j) for(int i=1; i<6; ++i) if(NDIM>1) emfz[i][j][k]=0.25*(fluxes[0][i][j][k][6]+fluxes[0][i][j-1][k][6]-fluxes[1][i][j][k][5]-fluxes[1][i-1][j][k][5]);
        for(int k2=0; k2<(NDIM>2?2:1); ++k2) for(int j2=0; j2<(NDIM>1?2:1); ++j2) for(int i2=0; i2<2; ++i2) {
            int i=2+i2, j=2+j2, k=2+k2, icp=1+i2+2*j2+4*k2, idc=grid_.ncoarse+(icp-1)*grid_.ngridmax+igrid;
            if(icp<=constants::twotondim) {
                for(int iv=1; iv<=5; ++iv) for(int d=0; d<NDIM; ++d) grid_.unew(idc, iv) += (fluxes[d][i-(d==0)][j-(d==1)][k-(d==2)][iv-1] - fluxes[d][i][j][k][iv-1])*dt_dx;
                if(NDIM>1){ grid_.unew(idc, 6) += (emfz[i][j][k]-emfz[i][j+1][k])*dt_dx; grid_.unew(idc, 7) += (emfz[i+1][j][k]-emfz[i][j][k])*dt_dx; grid_.unew(idc, 9) += (emfz[i+1][j][k]-emfz[i+1][j+1][k])*dt_dx; grid_.unew(idc, 10) += (emfz[i+1][j+1][k]-emfz[i][j+1][k])*dt_dx; }
            }
        }
    }
}

void MhdSolver::llf(const real_t* ql, const real_t* qr, real_t* f, real_t gamma) {
    real_t fl[9], fr[9], cl, cr; find_mhd_flux(ql, nullptr, fl, gamma); find_mhd_flux(qr, nullptr, fr, gamma); find_speed_fast(ql, cl, gamma); find_speed_fast(qr, cr, gamma);
    real_t sm = std::max(std::abs(ql[2])+cl, std::abs(qr[2])+cr), et = 1.0/(gamma-1.0);
    auto gc = [&](const real_t* q, real_t* c) { c[0]=q[0]; c[1]=q[0]*q[2]; c[2]=q[0]*q[4]; c[3]=q[0]*q[6]; c[4]=q[1]*et+0.5*q[0]*(q[2]*q[2]+q[4]*q[4]+q[6]*q[6])+0.5*(q[3]*q[3]+q[5]*q[5]+q[7]*q[7]); c[5]=q[3]; c[6]=q[5]; c[7]=q[7]; c[8]=q[1]*et; };
    real_t cll[9], crr[9]; gc(ql, cll); gc(qr, crr);
    for(int i=0; i<9; ++i) f[i] = 0.5*(fl[i]+fr[i]-sm*(crr[i]-cll[i]));
}

void MhdSolver::hlld(const real_t* ql_in, const real_t* qr_in, real_t* fgdnv, real_t gamma) {
    real_t entho = 1.0/(gamma-1.0), half = 0.5;
    real_t A = half*(ql_in[3]+qr_in[3]), sgnm = (A>=0)?1.0:-1.0;
    real_t ql[8], qr[8]; for(int i=0; i<8; ++i){ ql[i]=ql_in[i]; qr[i]=qr_in[i]; } ql[3]=A; qr[3]=A;
    real_t rl=ql[0],Pl=ql[1],ul=ql[2],vl=ql[4],Bl=ql[5],wl=ql[6],Cl=ql[7];
    real_t rr=qr[0],Pr=qr[1],ur=qr[2],vr=qr[4],Br=qr[5],wr=qr[6],Cr=qr[7];
    real_t ecinl=half*rl*(ul*ul+vl*vl+wl*wl), emagl=half*(A*A+Bl*Bl+Cl*Cl), etotl=Pl*entho+ecinl+emagl, Ptotl=Pl+emagl, vdotBl=ul*A+vl*Bl+wl*Cl, eintl=Pl*entho;
    real_t ecinr=half*rr*(ur*ur+vr*vr+wr*wr), emagr=half*(A*A+Br*Br+Cr*Cr), etotr=Pr*entho+ecinr+emagr, Ptotr=Pr+emagr, vdotBr=ur*A+vr*Br+wr*Cr, eintr=Pr*entho;
    real_t cl, cr; find_speed_fast(ql, cl, gamma); find_speed_fast(qr, cr, gamma);
    real_t SL=std::min(ul,ur)-std::max(cl,cr), SR=std::max(ul,ur)+std::max(cl,cr);
    real_t rcl=rl*(ul-SL), rcr=rr*(SR-ur), ustar=(rcr*ur+rcl*ul+(Ptotl-Ptotr))/(rcr+rcl), Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl);
    real_t rstarl=rl*(SL-ul)/(SL-ustar), estl=rl*(SL-ul)*(SL-ustar)-A*A, ell=rl*(SL-ul)*(SL-ul)-A*A, vstl, Bstl, wstl, Cstl;
    if(std::abs(estl)<1e-4*A*A){ vstl=vl; Bstl=Bl; wstl=wl; Cstl=Cl; } else { vstl=vl-A*Bl*(ustar-ul)/estl; Bstl=Bl*ell/estl; wstl=wl-A*Cl*(ustar-ul)/estl; Cstl=Cl*ell/estl; }
    real_t vdbstl=ustar*A+vstl*Bstl+wstl*Cstl, etstl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar+A*(vdotBl-vdbstl))/(SL-ustar), sqrl=std::sqrt(rstarl), calfl=std::abs(A)/sqrl, SAL=ustar-calfl;
    real_t rstarr=rr*(SR-ur)/(SR-ustar), estr=rr*(SR-ur)*(SR-ustar)-A*A, err=rr*(SR-ur)*(SR-ur)-A*A, vstr, Bstr, wstr, Cstr;
    if(std::abs(estr)<1e-4*A*A){ vstr=vr; Bstr=Br; wstr=wr; Cstr=Cr; } else { vstr=vr-A*Br*(ustar-ur)/estr; Bstr=Br*err/estr; wstr=wr-A*Cr*(ustar-ur)/estr; Cstr=Cr*err/estr; }
    real_t vdbstr=ustar*A+vstr*Bstr+wstr*Cstr, etstr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar+A*(vdotBr-vdbstr))/(SR-ustar), sqrr=std::sqrt(rstarr), calfr=std::abs(A)/sqrr, SAR=ustar+calfr;
    real_t vss=(sqrl*vstl+sqrr*vstr+sgnm*(Bstr-Bstl))/(sqrl+sqrr), wss=(sqrl*wstl+sqrr*wstr+sgnm*(Cstr-Cstl))/(sqrl+sqrr), Bss=(sqrl*Bstl+sqrr*Bstr+sgnm*sqrl*sqrr*(vstr-vstl))/(sqrl+sqrr), Css=(sqrl*Cstl+sqrr*Cstr+sgnm*sqrl*sqrr*(wstr-wstl))/(sqrl+sqrr);
    real_t vdbss=ustar*A+vss*Bss+wss*Css, etssl=etstl-sgnm*sqrl*(vdbstl-vdbss), etssr=etstr+sgnm*sqrr*(vdbstr-vdbss);
    real_t ro, uo, vo, wo, Bo, Co, Ptoto, etoto, vdb, ei;
    if(SL>0){ ro=rl; uo=ul; vo=vl; wo=wl; Bo=Bl; Co=Cl; Ptoto=Ptotl; etoto=etotl; vdb=vdotBl; ei=eintl; }
    else if(SAL>0){ ro=rstarl; uo=ustar; vo=vstl; wo=wstl; Bo=Bstl; Co=Cstl; Ptoto=Ptotstar; etoto=etstl; vdb=vdbstl; ei=eintl*(SL-ul)/(SL-ustar); }
    else if(ustar>0){ ro=rstarl; uo=ustar; vo=vss; wo=wss; Bo=Bss; Co=Css; Ptoto=Ptotstar; etoto=etssl; vdb=vdbss; ei=eintl*(SL-ul)/(SL-ustar); }
    else if(SAR>0){ ro=rstarr; uo=ustar; vo=vss; wo=wss; Bo=Bss; Co=Css; Ptoto=Ptotstar; etoto=etssr; vdb=vdbss; ei=eintr*(SR-ur)/(SR-ustar); }
    else if(SR>0){ ro=rstarr; uo=ustar; vo=vstr; wo=wstr; Bo=Bstr; Co=Cstr; Ptoto=Ptotstar; etoto=etstr; vdb=vdbstr; ei=eintr*(SR-ur)/(SR-ustar); }
    else { ro=rr; uo=ur; vo=vr; wo=wr; Bo=Br; Co=Cr; Ptoto=Ptotr; etoto=etotr; vdb=vdotBr; ei=eintr; }
    fgdnv[0]=ro*uo; fgdnv[1]=(etoto+Ptoto)*uo-A*vdb; fgdnv[2]=ro*uo*uo+Ptoto-A*A; fgdnv[3]=0; fgdnv[4]=ro*uo*vo-A*Bo; fgdnv[5]=Bo*uo-A*vo; fgdnv[6]=ro*uo*wo-A*Co; fgdnv[7]=Co*uo-A*wo; fgdnv[8]=uo*ei;
}

void MhdSolver::find_mhd_flux(const real_t* q, real_t* c, real_t* f, real_t gamma) {
    real_t et = 1.0/(gamma-1.0), d=q[0], p=q[1], u=q[2], A=q[3], v=q[4], B=q[5], w=q[6], C=q[7];
    real_t ek=0.5*d*(u*u+v*v+w*w), em=0.5*(A*A+B*B+C*C), etot=p*et+ek+em, pt=p+em;
    f[0]=d*u; f[1]=(etot+pt)*u-A*(A*u+B*v+C*w); f[2]=d*u*u+pt-A*A; f[3]=0; f[4]=d*u*v-A*B; f[5]=B*u-A*v; f[6]=d*u*w-A*C; f[7]=C*u-A*w; f[8]=p*et*u;
}

void MhdSolver::find_speed_fast(const real_t* q, real_t& v, real_t gamma) {
    real_t d=q[0], p=q[1], A=q[3], B=q[5], C=q[7], b2=A*A+B*B+C*C, c2=gamma*p/d, d2=0.5*(b2/d+c2);
    v = std::sqrt(d2 + std::sqrt(std::max(0.0, d2*d2-c2*A*A/d)));
}

void MhdSolver::get_diagnostics(int il, real_t dx, real_t& mi, real_t& mv, real_t& mb) {
    mi=1e30; mv=0; for(int c=1; c<=grid_.ncpu; ++c) { int g=grid_.headl(c,il); while(g>0){ for(int i=1; i<=constants::twotondim; ++i){ int id=grid_.ncoarse+(i-1)*grid_.ngridmax+g; if(grid_.son[id]!=0) continue; real_t d=grid_.unew(id,1); mi=std::min(mi,d); real_t u=grid_.unew(id,2)/std::max(1e-10,d),v=grid_.unew(id,3)/std::max(1e-10,d),w=grid_.unew(id,4)/std::max(1e-10,d); mv=std::max(mv,std::sqrt(u*u+v*v+w*w)); } g=grid_.next[g-1]; } } mb=compute_max_div_b(il,dx);
}

} // namespace ramses
