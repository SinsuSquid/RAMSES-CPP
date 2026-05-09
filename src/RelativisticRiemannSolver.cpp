#include "ramses/RelativisticRiemannSolver.hpp"
#include <cmath>
#include <algorithm>

namespace ramses {

void RelativisticRiemannSolver::solve_llf(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, const std::string& eos) {
    // Basic LLF for RHD
    real_t fl[20], fr[20], ul[20], ur[20], cfl, cfr;
    find_rhd_flux(ql, ul, fl, gamma, eos);
    find_rhd_flux(qr, ur, fr, gamma, eos);
    find_speed_fast(ql, cfl, gamma, eos);
    find_speed_fast(qr, cfr, gamma, eos);
    real_t s_max = std::max(std::abs(ql[2]) + cfl, std::abs(qr[2]) + cfr); // Crude estimate
    for (int i = 0; i < 20; ++i) {
        flux[i] = 0.5 * (fl[i] + fr[i]) - 0.5 * s_max * (ur[i] - ul[i]);
    }
}

void RelativisticRiemannSolver::solve_hll(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, const std::string& eos) {
    real_t fl[20], fr[20], ul[20], ur[20], cfl, cfr;
    find_rhd_flux(ql, ul, fl, gamma, eos);
    find_rhd_flux(qr, ur, fr, gamma, eos);
    find_speed_fast(ql, cfl, gamma, eos);
    find_speed_fast(qr, cfr, gamma, eos);

    real_t vnl = ql[2], vnr = qr[2];
    real_t vlsq = ql[2]*ql[2] + ql[3]*ql[3] + ql[4]*ql[4];
    real_t vrsq = qr[2]*qr[2] + qr[3]*qr[3] + qr[4]*qr[4];
    real_t vplsq = ql[3]*ql[3] + ql[4]*ql[4];
    real_t vprsq = qr[3]*qr[3] + qr[4]*qr[4];

    real_t lorl = 1.0 / std::sqrt(1.0 - vlsq);
    real_t lorr = 1.0 / std::sqrt(1.0 - vrsq);

    real_t factorl = cfl * std::sqrt(1.0 - vnl*vnl - cfl*cfl * vplsq) / lorl;
    real_t factorr = cfr * std::sqrt(1.0 - vnr*vnr - cfr*cfr * vprsq) / lorr;

    real_t lleftp = (vnl * (1.0 - cfl*cfl) + factorl) / (1.0 - cfl*cfl * vlsq);
    real_t lleftm = (vnl * (1.0 - cfl*cfl) - factorl) / (1.0 - cfl*cfl * vlsq);
    real_t lrightp = (vnr * (1.0 - cfr*cfr) + factorr) / (1.0 - cfr*cfr * vrsq);
    real_t lrightm = (vnr * (1.0 - cfr*cfr) - factorr) / (1.0 - cfr*cfr * vrsq);

    real_t sr = std::max({0.0, lleftp, lrightp});
    real_t sl = std::min({0.0, lleftm, lrightm});

    if (sl >= 0.0) {
        for (int i = 0; i < 20; ++i) flux[i] = fl[i];
    } else if (sr <= 0.0) {
        for (int i = 0; i < 20; ++i) flux[i] = fr[i];
    } else {
        real_t inv_srs = 1.0 / (sr - sl);
        for (int i = 0; i < 20; ++i) {
            flux[i] = (sr * fl[i] - sl * fr[i] + sr * sl * (ur[i] - ul[i])) * inv_srs;
        }
    }
}

void RelativisticRiemannSolver::solve_hllc(const real_t ql[], const real_t qr[], real_t flux[], real_t gamma, const std::string& eos) {
    real_t fl[20], fr[20], ul[20], ur[20], cfl, cfr;
    find_rhd_flux(ql, ul, fl, gamma, eos);
    find_rhd_flux(qr, ur, fr, gamma, eos);
    find_speed_fast(ql, cfl, gamma, eos);
    find_speed_fast(qr, cfr, gamma, eos);

    real_t vnl = ql[2], vnr = qr[2];
    real_t vlsq = ql[2]*ql[2] + ql[3]*ql[3] + ql[4]*ql[4];
    real_t vrsq = qr[2]*qr[2] + qr[3]*qr[3] + qr[4]*qr[4];
    real_t vplsq = ql[3]*ql[3] + ql[4]*ql[4];
    real_t vprsq = qr[3]*qr[3] + qr[4]*qr[4];

    real_t lorl = 1.0 / std::sqrt(1.0 - vlsq);
    real_t lorr = 1.0 / std::sqrt(1.0 - vrsq);

    real_t factorl = cfl * std::sqrt(1.0 - vnl*vnl - cfl*cfl * vplsq) / lorl;
    real_t factorr = cfr * std::sqrt(1.0 - vnr*vnr - cfr*cfr * vprsq) / lorr;

    real_t lleftp = (vnl * (1.0 - cfl*cfl) + factorl) / (1.0 - cfl*cfl * vlsq);
    real_t lleftm = (vnl * (1.0 - cfl*cfl) - factorl) / (1.0 - cfl*cfl * vlsq);
    real_t lrightp = (vnr * (1.0 - cfr*cfr) + factorr) / (1.0 - cfr*cfr * vrsq);
    real_t lrightm = (vnr * (1.0 - cfr*cfr) - factorr) / (1.0 - cfr*cfr * vrsq);

    real_t sr = std::max({lleftp, lrightp, 0.0});
    real_t sl = std::min({lleftm, lrightm, 0.0});

    if (sl >= 0.0) {
        for (int i = 0; i < 20; ++i) flux[i] = fl[i];
    } else if (sr <= 0.0) {
        for (int i = 0; i < 20; ++i) flux[i] = fr[i];
    } else {
        real_t inv_srs = 1.0 / (sr - sl);
        real_t fhll[20], uhll[20];
        for (int i = 0; i < 20; ++i) {
            fhll[i] = (sr * fl[i] - sl * fr[i] + sr * sl * (ur[i] - ul[i])) * inv_srs;
            uhll[i] = (sr * ur[i] - sl * ul[i] + fl[i] - fr[i]) * inv_srs;
        }

        real_t a = fhll[1]; // Pressure is related to momentum flux in RHD HLLC
        real_t b = -(uhll[1] + fhll[2]);
        real_t c = uhll[2];
        real_t sstar;
        if (b > 0) {
            real_t quad = -0.5 * (b + std::sqrt(std::max(0.0, b*b - 4.0*a*c)));
            sstar = c / quad;
        } else {
            real_t quad = -0.5 * (b - std::sqrt(std::max(0.0, b*b - 4.0*a*c)));
            sstar = c / quad;
        }

        if (sstar >= 0.0) {
            real_t ps = -fhll[1] * sstar + fhll[2];
            real_t den = 1.0 / (sl - sstar);
            real_t usl[20];
            usl[0] = ul[0] * (sl - ql[2]) * den;
            usl[2] = (ul[2] * (sl - ql[2]) + ps - ql[1]) * den;
            usl[3] = ul[3] * (sl - ql[2]) * den;
            usl[4] = ul[4] * (sl - ql[2]) * den;
            usl[1] = (ul[1] * (sl - ql[2]) + ps * sstar - ql[1]*ql[2]) * den;
            for (int i = 0; i < 5; ++i) flux[i] = sl * (usl[i] - ul[i]) + fl[i];
            // Passive scalars
            for (int i = 5; i < 20; ++i) flux[i] = (flux[0] > 0) ? flux[0] * ql[i] : flux[0] * qr[i];
        } else {
            real_t ps = -fhll[1] * sstar + fhll[2];
            real_t den = 1.0 / (sr - sstar);
            real_t usr[20];
            usr[0] = ur[0] * (sr - qr[2]) * den;
            usr[2] = (ur[2] * (sr - qr[2]) + ps - qr[1]) * den;
            usr[3] = ur[3] * (sr - qr[2]) * den;
            usr[4] = ur[4] * (sr - qr[2]) * den;
            usr[1] = (ur[1] * (sr - qr[2]) + ps * sstar - qr[1]*qr[2]) * den;
            for (int i = 0; i < 5; ++i) flux[i] = sr * (usr[i] - ur[i]) + fr[i];
            for (int i = 5; i < 20; ++i) flux[i] = (flux[0] > 0) ? flux[0] * ql[i] : flux[0] * qr[i];
        }
    }
}

void RelativisticRiemannSolver::find_rhd_flux(const real_t q[], real_t u[], real_t f[], real_t gamma, const std::string& eos) {
    real_t d = q[0], p = q[4], vx = q[1], vy = q[2], vz = q[3];
    real_t v2 = vx*vx + vy*vy + vz*vz;
    real_t lor = 1.0 / std::sqrt(1.0 - v2);
    real_t enth;
    if (eos == "TM") {
        real_t tau = p / d;
        enth = d * (2.5 * tau + 1.5 * std::sqrt(tau*tau + 4.0/9.0));
    } else {
        enth = d + p / (gamma - 1.0) + p;
    }

    u[0] = d * lor;
    u[1] = enth * lor*lor - p;
    u[2] = enth * lor*lor * vx;
    u[3] = enth * lor*lor * vy;
    u[4] = enth * lor*lor * vz;
    for (int i = 5; i < 20; ++i) u[i] = d * q[i] * lor;

    f[0] = d * vx * lor;
    f[1] = enth * lor*lor * vx;
    f[2] = enth * lor*lor * vx*vx + p;
    f[3] = enth * lor*lor * vx * vy;
    f[4] = enth * lor*lor * vx * vz;
    for (int i = 5; i < 20; ++i) f[i] = d * vx * q[i] * lor;
}

void RelativisticRiemannSolver::find_speed_fast(const real_t q[], real_t& cs, real_t gamma, const std::string& eos) {
    real_t d = q[0], p = q[4];
    if (eos == "TM") {
        real_t tau = p / d;
        real_t enth = 2.5 * tau + 1.5 * std::sqrt(tau*tau + 4.0/9.0);
        cs = std::sqrt(tau / (3.0 * enth) * (5.0 * enth - 8.0 * tau) / (enth - tau));
    } else {
        real_t enth = d + p / (gamma - 1.0) + p;
        cs = std::sqrt(gamma * p / enth);
    }
}

} // namespace ramses