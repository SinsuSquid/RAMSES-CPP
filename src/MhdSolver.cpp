#include "ramses/MhdSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Muscl.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

namespace ramses {

void MhdSolver::set_unew(int ilevel) {
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                for (int ivar = 1; ivar <= grid_.nvar; ++ivar) {
                    grid_.unew(ind_cell, ivar) = grid_.uold(ind_cell, ivar);
                }
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
                int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                for (int ivar = 1; ivar <= grid_.nvar; ++ivar) {
                    grid_.uold(ind_cell, ivar) = grid_.unew(ind_cell, ivar);
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

real_t MhdSolver::compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor) {
    real_t dt_max = 1e30;
    real_t smallr = 1e-10;

    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                real_t q[8];
                for(int i=0; i<8; ++i) q[i] = grid_.uold(ind_cell, i+1);
                
                real_t vel_fast;
                find_speed_fast(q, vel_fast, gamma);

                real_t d = std::max(grid_.uold(ind_cell, 1), smallr);
                real_t u = grid_.uold(ind_cell, 2) / d;
                real_t v = grid_.uold(ind_cell, 3) / d;
                real_t w = grid_.uold(ind_cell, 4) / d;
                real_t vel = std::sqrt(u*u + v*v + w*w);
                
                real_t dt_cell = dx / (vel + vel_fast);
                dt_max = std::min(dt_max, dt_cell);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    return dt_max * courant_factor;
}

void MhdSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    int nx = params::nx * (1 << (ilevel - 1));
    int nvar = grid_.nvar;
    
    std::vector<real_t> q(nx * nvar);
    std::vector<real_t> dq(nx * nvar);
    std::vector<real_t> qL(nx * nvar);
    std::vector<real_t> qR(nx * nvar);
    std::vector<real_t> flux(nx * nvar);

    for (int i = 0; i < nx; ++i) {
        int ind_cell = i + 1;
        for (int iv = 0; iv < nvar; ++iv) {
            q[i * nvar + iv] = grid_.uold(ind_cell, iv + 1);
        }
    }

    Muscl::compute_slopes(q.data(), dq.data(), nx, nvar, 1);
    Muscl::reconstruct(q.data(), dq.data(), qL.data(), qR.data(), nx, nvar);

    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    for (int i = 0; i < nx - 1; ++i) {
        hlld(&qR[i * nvar], &qL[(i + 1) * nvar], &flux[i * nvar], gamma);
    }

    for (int i = 1; i < nx - 1; ++i) {
        int ind_cell = i + 1;
        for (int iv = 0; iv < nvar; ++iv) {
            grid_.unew(ind_cell, iv + 1) = q[i * nvar + iv] - (dt / dx) * (flux[i * nvar + iv] - flux[(i - 1) * nvar + iv]);
        }
    }
}

void MhdSolver::hlld(const real_t* qleft, const real_t* qright, real_t* fgdnv, real_t gamma) {
    // Ported HLLD implementation...
    real_t entho = 1.0 / (gamma - 1.0);
    real_t zero = 0.0;
    real_t half = 0.5;
    real_t one = 1.0;

    real_t A = half * (qleft[3] + qright[3]);
    real_t sgnm = (A >= 0) ? 1.0 : -1.0;

    real_t ql[8], qr[8];
    for(int i=0; i<8; ++i) { ql[i]=qleft[i]; qr[i]=qright[i]; }
    ql[3]=A; qr[3]=A;

    real_t rl = ql[0], Pl = ql[1], ul = ql[2];
    real_t vl = ql[4], Bl = ql[5], wl = ql[6], Cl = ql[7];
    real_t ecinl = half * (ul * ul + vl * vl + wl * wl) * rl;
    real_t emagl = half * (A * A + Bl * Bl + Cl * Cl);
    real_t etotl = Pl * entho + ecinl + emagl;
    real_t Ptotl = Pl + emagl;
    real_t vdotBl = ul * A + vl * Bl + wl * Cl;
    real_t eintl = Pl * entho;

    real_t rr = qr[0], Pr = qr[1], ur = qr[2];
    real_t vr = qr[4], Br = qr[5], wr = qr[6], Cr = qr[7];
    real_t ecinr = half * (ur * ur + vr * vr + wr * wr) * rr;
    real_t emagr = half * (A * A + Br * Br + Cr * Cr);
    real_t etotr = Pr * entho + ecinr + emagr;
    real_t Ptotr = Pr + emagr;
    real_t vdotBr = ur * A + vr * Br + wr * Cr;
    real_t eintr = Pr * entho;

    real_t cfastl, cfastr;
    find_speed_fast(ql, cfastl, gamma);
    find_speed_fast(qr, cfastr, gamma);

    real_t SL = std::min(ul, ur) - std::max(cfastl, cfastr);
    real_t SR = std::max(ul, ur) + std::max(cfastl, cfastr);

    real_t rcl = rl * (ul - SL);
    real_t rcr = rr * (SR - ur);
    real_t ustar = (rcr * ur + rcl * ul + (Ptotl - Ptotr)) / (rcr + rcl);
    real_t Ptotstar = (rcr * Ptotl + rcl * Ptotr + rcl * rcr * (ul - ur)) / (rcr + rcl);

    real_t rstarl = rl * (SL - ul) / (SL - ustar);
    real_t estar_l = rl * (SL - ul) * (SL - ustar) - A * A;
    real_t el = rl * (SL - ul) * (SL - ul) - A * A;
    real_t vstarl, Bstarl, wstarl, Cstarl;
    if (std::abs(estar_l) < 1e-6 * A * A) {
        vstarl = vl; Bstarl = Bl; wstarl = wl; Cstarl = Cl;
    } else {
        vstarl = vl - A * Bl * (ustar - ul) / estar_l;
        Bstarl = Bl * el / estar_l;
        wstarl = wl - A * Cl * (ustar - ul) / estar_l;
        Cstarl = Cl * el / estar_l;
    }
    real_t vdotBstarl = ustar * A + vstarl * Bstarl + wstarl * Cstarl;
    real_t etotstarl = ((SL - ul) * etotl - Ptotl * ul + Ptotstar * ustar + A * (vdotBl - vdotBstarl)) / (SL - ustar);
    real_t sqrrstarl = std::sqrt(rstarl);
    real_t calfvenl = std::abs(A) / sqrrstarl;
    real_t SAL = ustar - calfvenl;

    real_t rstarr = rr * (SR - ur) / (SR - ustar);
    real_t estar_r = rr * (SR - ur) * (SR - ustar) - A * A;
    real_t er = rr * (SR - ur) * (SR - ur) - A * A;
    real_t vstarr, Bstarr, wstarr, Cstarr;
    if (std::abs(estar_r) < 1e-6 * A * A) {
        vstarr = vr; Bstarr = Br; wstarr = wr; Cstarr = Cr;
    } else {
        vstarr = vr - A * Br * (ustar - ur) / estar_r;
        Bstarr = Br * er / estar_r;
        wstarr = wr - A * Cr * (ustar - ur) / estar_r;
        Cstarr = Cr * er / estar_r;
    }
    real_t vdotBstarr = ustar * A + vstarr * Bstarr + wstarr * Cstarr;
    real_t etotstarr = ((SR - ur) * etotr - Ptotr * ur + Ptotstar * ustar + A * (vdotBr - vdotBstarr)) / (SR - ustar);
    real_t sqrrstarr = std::sqrt(rstarr);
    real_t calfvenr = std::abs(A) / sqrrstarr;
    real_t SAR = ustar + calfvenr;

    real_t vstarstar = (sqrrstarl * vstarl + sqrrstarr * vstarr + sgnm * (Bstarr - Bstarl)) / (sqrrstarl + sqrrstarr);
    real_t wstarstar = (sqrrstarl * wstarl + sqrrstarr * wstarr + sgnm * (Cstarr - Cstarl)) / (sqrrstarl + sqrrstarr);
    real_t Bstarstar = (sqrrstarl * Bstarr + sqrrstarr * Bstarl + sgnm * sqrrstarl * sqrrstarr * (vstarr - vstarl)) / (sqrrstarl + sqrrstarr);
    real_t Cstarstar = (sqrrstarl * Cstarr + sqrrstarr * Cstarl + sgnm * sqrrstarl * sqrrstarr * (wstarr - wstarl)) / (sqrrstarl + sqrrstarr);
    real_t vdotBstarstar = ustar * A + vstarstar * Bstarstar + wstarstar * Cstarstar;
    real_t etotstarstarl = etotstarl - sgnm * sqrrstarl * (vdotBstarl - vdotBstarstar);
    real_t etotstarstarr = etotstarr + sgnm * sqrrstarr * (vdotBstarr - vdotBstarstar);

    real_t ro, uo, vo, wo, Bo, Co, Ptoto, etoto, vdotBo, einto;
    if (SL > 0.0) {
        ro = rl; uo = ul; vo = vl; wo = wl; Bo = Bl; Co = Cl; Ptoto = Ptotl; etoto = etotl; vdotBo = vdotBl; einto = eintl;
    } else if (SAL > 0.0) {
        ro = rstarl; uo = ustar; vo = vstarl; wo = wstarl; Bo = Bstarl; Co = Cstarl; Ptoto = Ptotstar; etoto = etotstarstarl; vdotBo = vdotBstarl; einto = eintl*(SL-ul)/(SL-ustar);
    } else if (ustar > 0.0) {
        ro = rstarl; uo = ustar; vo = vstarstar; wo = wstarstar; Bo = Bstarstar; Co = Cstarstar; Ptoto = Ptotstar; etoto = etotstarstarl; vdotBo = vdotBstarstar; einto = eintl*(SL-ul)/(SL-ustar);
    } else if (SAR > 0.0) {
        ro = rstarr; uo = ustar; vo = vstarstar; wo = wstarstar; Bo = Bstarstar; Co = Cstarstar; Ptoto = Ptotstar; etoto = etotstarstarr; vdotBo = vdotBstarstar; einto = eintr*(SR-ur)/(SR-ustar);
    } else if (SR > 0.0) {
        ro = rstarr; uo = ustar; vo = vstarr; wo = wstarr; Bo = Bstarr; Co = Cstarr; Ptoto = Ptotstar; etoto = etotstarr; vdotBo = vdotBstarr; einto = eintr*(SR-ur)/(SR-ustar);
    } else {
        ro = rr; uo = ur; vo = vr; wo = wr; Bo = Br; Co = Cr; Ptoto = Ptotr; etoto = etotr; vdotBo = vdotBr; einto = eintr;
    }

    fgdnv[0] = ro * uo;
    fgdnv[1] = (etoto + Ptoto) * uo - A * vdotBo;
    fgdnv[2] = ro * uo * uo + Ptoto - A * A;
    fgdnv[3] = 0.0;
    fgdnv[4] = ro * uo * vo - A * Bo;
    fgdnv[5] = Bo * uo - A * vo;
    fgdnv[6] = ro * uo * wo - A * Co;
    fgdnv[7] = Co * uo - A * wo;
    fgdnv[8] = uo * einto;
}

void MhdSolver::find_mhd_flux(const real_t* qvar, real_t* cvar, real_t* ff, real_t gamma) {
    real_t entho = 1.0 / (gamma - 1.0);
    real_t d = qvar[0], P = qvar[1], u = qvar[2], A = qvar[3];
    real_t v = qvar[4], B = qvar[5], w = qvar[6], C = qvar[7];
    real_t ecin = 0.5 * (u * u + v * v + w * w) * d;
    real_t emag = 0.5 * (A * A + B * B + C * C);
    real_t etot = P * entho + ecin + emag;
    real_t Ptot = P + emag;

    ff[0] = d * u;
    ff[1] = (etot + Ptot) * u - A * (A * u + B * v + C * w);
    ff[2] = d * u * u + Ptot - A * A;
    ff[3] = 0.0;
    ff[4] = d * u * v - A * B;
    ff[5] = B * u - A * v;
    ff[6] = d * u * w - A * C;
    ff[7] = C * u - A * w;
    ff[8] = P * entho * u;
}

void MhdSolver::find_speed_fast(const real_t* qvar, real_t& vel_info, real_t gamma) {
    real_t d = qvar[0], P = qvar[1], u = qvar[2], A = qvar[3];
    real_t v = qvar[4], B = qvar[5], w = qvar[6], C = qvar[7];
    real_t B2 = A * A + B * B + C * C;
    real_t c2 = gamma * P / d;
    real_t d2 = 0.5 * (B2 / d + c2);
    vel_info = std::sqrt(d2 + std::sqrt(std::max(0.0, d2 * d2 - c2 * A * A / d)));
}

} // namespace ramses
