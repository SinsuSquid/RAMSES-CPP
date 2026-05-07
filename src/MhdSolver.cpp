#include "ramses/MhdSolver.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Muscl.hpp"
#include "ramses/SlopeLimiter.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

namespace ramses {

MhdSolver::MhdSolver(AmrGrid& grid, Config& config) 
    : grid_(grid), config_(config), stencil_ptr_(std::make_unique<LocalStencil>()) {}

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
    real_t dt_max = 1e30, smallr = 1e-10;
    int n2d = (1 << NDIM);
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= n2d; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                real_t d = std::max(grid_.uold(idc, 1), smallr), u = grid_.uold(idc, 2)/d, v = grid_.uold(idc, 3)/d, w = grid_.uold(idc, 4)/d;
                real_t A = grid_.uold(idc, 6), B = grid_.uold(idc, 7), C = grid_.uold(idc, 8), etot = grid_.uold(idc, 5);
                real_t ekin = 0.5*d*(u*u+v*v+w*w), emag = 0.125*std::pow(A+grid_.uold(idc,grid_.nvar-2),2) + 0.125*std::pow(B+grid_.uold(idc,grid_.nvar-1),2) + 0.125*std::pow(C+grid_.uold(idc,grid_.nvar),2);
                real_t p = std::max((etot-ekin-emag)*(gamma-1.0), d*1e-10);
                real_t q_mhd[8] = {d, p, u, 0.5*(A+grid_.uold(idc,grid_.nvar-2)), v, 0.5*(B+grid_.uold(idc,grid_.nvar-1)), w, 0.5*(C+grid_.uold(idc,grid_.nvar))}, vel_fast;
                find_speed_fast(q_mhd, vel_fast, gamma);
                dt_max = std::min(dt_max, dx / (std::sqrt(u*u+v*v+w*w) + vel_fast));
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    return dt_max * courant_factor;
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
    for (int j = 0; j < 27; ++j) {
        int f_idx = constants::lll[pos-1][j], c_pos = constants::mmm[pos-1][j], if_n = nbors_27[f_idx - 1];
        int sz = j / 9; int sy = (j % 9) / 3; int sx = j % 3;
        for (int iv = 0; iv < 20; ++iv) stencil.uloc[sz][sy][sx][iv] = 0.0;
        if (if_n > 0 && grid_.cpu_map[if_n-1] == myid) {
            if (grid_.son[if_n-1] > 0) {
                int idc = grid_.ncoarse + (c_pos - 1) * grid_.ngridmax + grid_.son[if_n-1];
                for (int iv = 1; iv <= nvar; ++iv) stencil.uloc[sz][sy][sx][iv - 1] = grid_.uold(idc, iv);
                stencil.refined[sz][sy][sx] = (grid_.son[idc-1] > 0);
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
    q[4] = std::max((u[4] - ek - em) * (gamma - 1.0), d * 1e-10);
    q[5]=0.5*(bf[0][0]+bf[0][1]); q[6]=0.5*(bf[1][0]+bf[1][1]); q[7]=0.5*(bf[2][0]+bf[2][1]);
    for (int iv=8; iv<nvp; ++iv) q[iv]=u[iv]/d;
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
        int ni=i+(idim==0), nj=j+(idim==1), nk=k+(idim==2); real_t qL[8], qR[8], f[9];
        auto fl = [&](real_t* q, const real_t s[64], real_t bn) { q[0]=s[0]; q[1]=s[4]; q[2]=s[1+idim]; q[3]=bn; q[4]=s[1+(idim+1)%3]; q[5]=s[5+(idim+1)%3]; q[6]=s[1+(idim+2)%3]; q[7]=s[5+(idim+2)%3]; };
        fl(qL, qm[i][j][k][idim], bfloc[i][j][k][idim][1]); fl(qR, qp[ni][nj][nk][idim], bfloc[ni][nj][nk][idim][0]);
        if (config_.get("hydro_params", "riemann", "hlld") == "llf") llf(qL, qR, f, gamma); else hlld(qL, qR, f, gamma);
        flux[i][j][k][0]=f[0]; flux[i][j][k][1]=f[1]; flux[i][j][k][2]=f[2]; flux[i][j][k][4]=f[4]; flux[i][j][k][6]=f[6]; flux[i][j][k][5]=f[5]; flux[i][j][k][7]=f[7];
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
                    real_t fL[64], fR[64]; for(int iv=0; iv<20; ++iv) { fL[iv]=fluxes[d][il][jl][kl][iv]; fR[iv]=fluxes[d][i][j][k][iv]; }
                    grid_.unew(idc, 1) += (fL[0]-fR[0])*dt_dx; grid_.unew(idc, 5) += (fL[1]-fR[1])*dt_dx;
                    grid_.unew(idc, 1+d+1) += (fL[2]-fR[2])*dt_dx; grid_.unew(idc, 1+(d+1)%3+1) += (fL[4]-fR[4])*dt_dx; grid_.unew(idc, 1+(d+2)%3+1) += (fL[6]-fR[6])*dt_dx;
                }
                if(NDIM>1){
                    grid_.unew(idc, 6) += (emfz[i-1][j][k]-emfz[i-1][j+1][k])*dt_dx; grid_.unew(idc, 7) += (emfz[i][j-1][k]-emfz[i+1][j-1][k])*dt_dx;
                    grid_.unew(idc, 9) += (emfz[i][j][k]-emfz[i][j+1][k])*dt_dx; grid_.unew(idc, 10) += (emfz[i][j][k]-emfz[i+1][j][k])*dt_dx;
                }
            }
        }
        // MHD Refluxing (EMF Averaging at Coarse-Fine interfaces)
        if(ilevel > 1 && NDIM > 1) {
            int ind_father = grid_.father[igrid-1];
            int ign[7]; grid_.get_nbor_grids(ind_father, ign);
            for(int ic=1; ic<=constants::twotondim; ic++) {
                int idc = grid_.ncoarse+(ic-1)*grid_.ngridmax+ind_father;
                if(grid_.son[idc-1] == igrid) {
                    // This oct is a child of idc. Apply EMF averages to coarse neighbors of idc.
                    // Simplified: for 2D, we only have emfz.
                    real_t weight = 0.25 * 1.0; // Assume leaf-leaf interface for now
                    real_t df = (emfz[2][2][2] + emfz[2][2][3]) * 0.25 * dt_dx; // Bottom-left corner EMF
                    // ... (Full EMF refluxing logic requires nbors_father_cells map)
                }
            }
        }
    }
}

void MhdSolver::llf(const real_t* ql, const real_t* qr, real_t* f, real_t gamma) {
    real_t fl[9], fr[9], cl, cr, et=1.0/(gamma-1.0); find_mhd_flux(ql, nullptr, fl, gamma); find_mhd_flux(qr, nullptr, fr, gamma); find_speed_fast(ql, cl, gamma); find_speed_fast(qr, cr, gamma);
    real_t sm = std::max(std::abs(ql[2])+cl, std::abs(qr[2])+cr);
    auto gc = [&](const real_t* q, real_t* c) { c[0]=q[0]; c[1]=q[0]*q[2]; c[2]=q[0]*q[4]; c[3]=q[0]*q[6]; c[4]=q[1]*et+0.5*q[0]*(q[2]*q[2]+q[4]*q[4]+q[6]*q[6])+0.5*(q[3]*q[3]+q[5]*q[5]+q[7]*q[7]); c[5]=q[3]; c[6]=q[5]; c[7]=q[7]; c[8]=q[1]*et; };
    real_t cll[9], crr[9]; gc(ql, cll); gc(qr, crr);
    for(int i=0; i<9; ++i) f[i] = 0.5*(fl[i]+fr[i]-sm*(crr[i]-cll[i]));
}

void MhdSolver::hlld(const real_t* ql_in, const real_t* qr_in, real_t* fgdnv, real_t gamma) {
    real_t entho = 1.0/(gamma-1.0), half = 0.5, A = half*(ql_in[3]+qr_in[3]), sgnm = (A>=0)?1.0:-1.0;
    real_t ql[8], qr[8]; for(int i=0; i<8; ++i){ ql[i]=ql_in[i]; qr[i]=qr_in[i]; } ql[3]=A; qr[3]=A;
    real_t rl=ql[0],Pl=ql[1],ul=ql[2],vl=ql[4],Bl=ql[5],wl=ql[6],Cl=ql[7];
    real_t rr=qr[0],Pr=qr[1],ur=qr[2],vr=qr[4],Br=qr[5],wr=qr[6],Cr=qr[7];
    real_t el=half*rl*(ul*ul+vl*vl+wl*wl), emagl=half*(A*A+Bl*Bl+Cl*Cl), etl=Pl*entho+el+emagl, Ptl=Pl+emagl, vBl=ul*A+vl*Bl+wl*Cl;
    real_t er=half*rr*(ur*ur+vr*vr+wr*wr), emagr=half*(A*A+Br*Br+Cr*Cr), etr=Pr*entho+er+emagr, Ptr=Pr+emagr, vBr=ur*A+vr*Br+wr*Cr;
    real_t cl, cr; find_speed_fast(ql, cl, gamma); find_speed_fast(qr, cr, gamma);
    real_t SL=std::min(ul,ur)-std::max(cl,cr), SR=std::max(ul,ur)+std::max(cl,cr);
    real_t rcl=rl*(ul-SL), rcr=rr*(SR-ur), ust=(rcr*ur+rcl*ul+(Ptl-Ptr))/(rcr+rcl), Ptst=(rcr*Ptl+rcl*Ptr+rcl*rcr*(ul-ur))/(rcr+rcl);
    real_t rstl=rl*(SL-ul)/(SL-ust), estl=rl*(SL-ul)*(SL-ust)-A*A, ell=rl*(SL-ul)*(SL-ul)-A*A, vstl, Bstl, wstl, Cstl;
    if(std::abs(estl)<1e-4*A*A){ vstl=vl; Bstl=Bl; wstl=wl; Cstl=Cl; } else { vstl=vl-A*Bl*(ust-ul)/estl; Bstl=Bl*ell/estl; wstl=wl-A*Cl*(ust-ul)/estl; Cstl=Cl*ell/estl; }
    real_t vdbstl=ust*A+vstl*Bstl+wstl*Cstl, etstl=((SL-ul)*etl-Ptl*ul+Ptst*ust+A*(vBl-vdbstl))/(SL-ust), sqrl=std::sqrt(rstl), calfl=std::abs(A)/sqrl, SAL=ust-calfl;
    real_t rstr=rr*(SR-ur)/(SR-ust), estr=rr*(SR-ur)*(SR-ust)-A*A, err=rr*(SR-ur)*(SR-ur)-A*A, vstr, Bstr, wstr, Cstr;
    if(std::abs(estr)<1e-4*A*A){ vstr=vr; Bstr=Br; wstr=wr; Cstr=Cr; } else { vstr=vr-A*Br*(ust-ur)/estr; Bstr=Br*err/estr; wstr=wr-A*Cr*(ust-ur)/estr; Cstr=Cr*err/estr; }
    real_t vdbstr=ust*A+vstr*Bstr+wstr*Cstr, etstr=((SR-ur)*etr-Ptr*ur+Ptst*ust+A*(vBr-vdbstr))/(SR-ust), sqrr=std::sqrt(rstr), calfr=std::abs(A)/sqrr, SAR=ust+calfr;
    real_t vss=(sqrl*vstl+sqrr*vstr+sgnm*(Bstr-Bstl))/(sqrl+sqrr), wss=(sqrl*wstl+sqrr*wstr+sgnm*(Cstr-Cstl))/(sqrl+sqrr), Bss=(sqrl*Bstl+sqrr*Bstr+sgnm*sqrl*sqrr*(vstr-vstl))/(sqrl+sqrr), Css=(sqrl*Cstl+sqrr*Cstr+sgnm*sqrl*sqrr*(wstr-wstl))/(sqrl+sqrr);
    real_t vdbss=ust*A+vss*Bss+wss*Css, etssl=etstl-sgnm*sqrl*(vdbstl-vdbss), etssr=etstr+sgnm*sqrr*(vdbstr-vdbss);
    real_t ro, uo, vo, wo, Bo, Co, Ptoto, etoto, vdb, ei, eintl=Pl*entho, eintr=Pr*entho;
    if(SL>0){ ro=rl; uo=ul; vo=vl; wo=wl; Bo=Bl; Co=Cl; Ptoto=Ptl; etoto=etl; vdb=vBl; ei=eintl; }
    else if(SAL>0){ ro=rstl; uo=ust; vo=vstl; wo=wstl; Bo=Bstl; Co=Cstl; Ptoto=Ptst; etoto=etstl; vdb=vdbstl; ei=eintl*(SL-ul)/(SL-ust); }
    else if(ust>0){ ro=rstl; uo=ust; vo=vss; wo=wss; Bo=Bss; Co=Css; Ptoto=Ptst; etoto=etssl; vdb=vdbss; ei=eintl*(SL-ul)/(SL-ust); }
    else if(SAR>0){ ro=rstr; uo=ust; vo=vss; wo=wss; Bo=Bss; Co=Css; Ptoto=Ptst; etoto=etssr; vdb=vdbss; ei=eintr*(SR-ur)/(SR-ust); }
    else if(SR>0){ ro=rstr; uo=ust; vo=vstr; wo=wstr; Bo=Bstr; Co=Cstr; Ptoto=Ptst; etoto=etstr; vdb=vdbstr; ei=eintr*(SR-ur)/(SR-ust); }
    else { ro=rr; uo=ur; vo=vr; wo=wr; Bo=Br; Co=Cr; Ptoto=Ptr; etoto=etr; vdb=vBr; ei=eintr; }
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

void MhdSolver::get_diagnostics(int il, real_t dx, real_t& mi, real_t& mv, real_t& mb) {}

} // namespace ramses
