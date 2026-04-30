#include "ramses/Initializer.hpp"
#include "ramses/Parameters.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>

namespace ramses {

void Initializer::apply_all() {
    std::string condinit_kind = config_.get("init_params", "condinit_kind", "");
    std::cout << "[Initializer] condinit_kind=\"" << condinit_kind << "\"" << std::endl;
    if (condinit_kind == "ana_disk_potential" || condinit_kind == "'ana_disk_potential'") {
        ana_disk_potential_condinit();
    } else if (condinit_kind == "orzag_tang" || condinit_kind == "'orzag_tang'") {
        orzag_tang_condinit();
    } else {
        region_condinit();
    }
}

void Initializer::orzag_tang_condinit() {
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.6666667);
    real_t pi = std::acos(-1.0);
    real_t B0 = 1.0 / std::sqrt(4.0 * pi);

    for (int i = 1; i <= grid_.ncell; ++i) {
        real_t x_cell[3];
        grid_.get_cell_center(i, x_cell);
        
        int ilevel = 1;
        if (i > grid_.ncoarse) {
            int igrid = (i - grid_.ncoarse - 1) % grid_.ngridmax + 1;
            int ifather = grid_.father[igrid - 1];
            ilevel = 2;
            while (ifather > grid_.ncoarse) {
                ilevel++;
                int igrid_father = (ifather - grid_.ncoarse - 1) % grid_.ngridmax + 1;
                ifather = grid_.father[igrid_father - 1];
            }
        }
        real_t dx = 1.0 / static_cast<real_t>(1 << (ilevel - 1));

        real_t xc = x_cell[0], yc = x_cell[1];
        real_t xl = xc - 0.5 * dx, xr = xc + 0.5 * dx;
        real_t yl = yc - 0.5 * dx, yr = yc + 0.5 * dx;

        real_t d = 25.0 / (36.0 * pi);
        real_t u = -std::sin(2.0 * pi * yc);
        real_t v = std::sin(2.0 * pi * xc);
        real_t w = 0.0;
        real_t p = 5.0 / (12.0 * pi);

        auto get_Ar = [&](real_t xx, real_t yy) {
            return B0 * (std::cos(4.0 * pi * xx) / (4.0 * pi) + std::cos(2.0 * pi * yy) / (2.0 * pi));
        };

        real_t Bx_l = (get_Ar(xl, yr) - get_Ar(xl, yl)) / dx;
        real_t Bx_r = (get_Ar(xr, yr) - get_Ar(xr, yl)) / dx;
        real_t By_l = (get_Ar(xl, yl) - get_Ar(xr, yl)) / dx;
        real_t By_r = (get_Ar(xl, yr) - get_Ar(xr, yr)) / dx;

        grid_.uold(i, 1) = d;
        grid_.uold(i, 2) = d * u;
        grid_.uold(i, 3) = d * v;
        grid_.uold(i, 4) = d * w;
        
        real_t Bx_avg = 0.5 * (Bx_l + Bx_r);
        real_t By_avg = 0.5 * (By_l + By_r);
        real_t Bz_avg = 0.0;

        real_t e_kin = 0.5 * d * (u*u + v*v + w*w);
        real_t e_mag = 0.5 * (Bx_avg*Bx_avg + By_avg*By_avg + Bz_avg*Bz_avg);
        grid_.uold(i, 5) = p / (gamma - 1.0) + e_kin + e_mag;

#ifdef MHD
        grid_.uold(i, 6) = Bx_l;
        grid_.uold(i, 7) = By_l;
        grid_.uold(i, 8) = 0.0;
        grid_.uold(i, grid_.nvar - 2) = Bx_r;
        grid_.uold(i, grid_.nvar - 1) = By_r;
        grid_.uold(i, grid_.nvar) = 0.0;
#endif
        grid_.cpu_map[i] = 1;
    }
    std::cout << "[Initializer] Applied Orszag-Tang ICs." << std::endl;
}

void Initializer::ana_disk_potential_condinit() {
    region_condinit();
    std::string param_str = config_.get("poisson_params", "gravity_params", "");
    real_t a1 = 1.42e-3, a2 = 5.49e-4, z0 = 0.18e3;
    if (!param_str.empty()) {
        std::replace(param_str.begin(), param_str.end(), 'd', 'e');
        std::replace(param_str.begin(), param_str.end(), 'D', 'e');
        std::replace(param_str.begin(), param_str.end(), ',', ' ');
        std::stringstream ss(param_str); ss >> a1 >> a2 >> z0;
    }
    real_t sl = config_.get_double("units_params", "units_length", 1.0), st = config_.get_double("units_params", "units_time", 1.0), sv = sl / st;
    real_t kpc2cm = 3.085677581282e21, pc2cm = 3.085677581282e18, Myr2sec = 3.15576e13;
    a1 = a1 * kpc2cm / (Myr2sec * Myr2sec) / sl * (st * st); a2 = a2 / (Myr2sec * Myr2sec) * (st * st); z0 = z0 * pc2cm / sl;
    real_t t0 = 8000.0, mu = config_.get_double("cooling_params", "mu_gas", 1.4), gam = config_.get_double("hydro_params", "gamma", 1.4);
    real_t kB = 1.3806490e-16, mH = 1.6605390e-24, cs2 = (kB * t0 / (mu * mH)) / (sv * sv);
    real_t a1r = a1 / (4.0 * M_PI * cs2) * z0 * z0, a2r = a2 / (2.0 * M_PI * cs2);

    for (int ic = 1; ic <= grid_.ncpu; ++ic) {
        for (int il = 1; il <= grid_.nlevelmax; ++il) {
            int ig = grid_.headl(ic, il);
            while (ig > 0) {
                for (int c = 1; ic <= constants::twotondim; ++c) {
                    int id = grid_.ncoarse + (c - 1) * grid_.ngridmax + ig;
                    real_t x_cell[3]; grid_.get_cell_center(id, x_cell);
                    real_t x = x_cell[NDIM-1]; // Use last dimension for vertical
                    real_t z = (x - 0.5) * params::boxlen;
                    real_t rho = a1r / std::pow(1.0 + std::pow(z / z0, 2), 1.5) + a2r;
                    grid_.uold(id, 1) = rho;
                    for (int d = 1; d <= NDIM; ++d) grid_.uold(id, 1 + d) = 0.0;
                    grid_.uold(id, 5) = (rho * cs2) / (gam - 1.0);
                }
                ig = grid_.next[ig - 1];
            }
        }
    }
}

void Initializer::region_condinit() {
    int nreg = config_.get_int("init_params", "nregion", 1);
    real_t gam = config_.get_double("hydro_params", "gamma", 1.4);
    int nener = config_.get_int("hydro_params", "nener", 0);
    
    auto parse = [&](const std::string& k) {
        std::vector<real_t> v; std::string s = config_.get("init_params", k, "");
        if (s.empty()) return v;
        std::replace(s.begin(), s.end(), 'd', 'e'); std::replace(s.begin(), s.end(), 'D', 'e'); std::replace(s.begin(), s.end(), ',', ' ');
        std::stringstream ss(s); real_t val; while (ss >> val) v.push_back(val); return v;
    };

    auto drs = parse("d_region"), prs = parse("p_region"), urs = parse("u_region"), vrs = parse("v_region"), wrs = parse("w_region");
    auto Ars = parse("A_region"), Brs = parse("B_region"), Crs = parse("C_region"), xcs = parse("x_center"), lxs = parse("length_x");
    auto ycs = parse("y_center"), lys = parse("length_y"), zcs = parse("z_center"), lzs = parse("length_z");

    if (nener == 0) {
        for (int k = 1; k <= 10; ++k) if (!config_.get("init_params", "prad_region(" + std::to_string(k) + ",1)", "").empty()) nener = std::max(nener, k);
    }
    std::vector<std::vector<real_t>> prads(nreg, std::vector<real_t>(nener, 0.0));
    for (int r = 0; r < nreg; ++r) for (int e = 0; e < nener; ++e) {
        prads[r][e] = config_.get_double("init_params", "prad_region(" + std::to_string(e + 1) + "," + std::to_string(r + 1) + ")", 0.0);
    }

    real_t db = !drs.empty() ? drs[0] : 1.0, pb = !prs.empty() ? prs[0] : 1e-5, ub = !urs.empty() ? urs[0] : 0.0, vb = !vrs.empty() ? vrs[0] : 0.0, wb = !wrs.empty() ? wrs[0] : 0.0;
    real_t Ab = !Ars.empty() ? Ars[0] : 0.0, Bb = !Brs.empty() ? Brs[0] : 0.0, Cb = !Crs.empty() ? Crs[0] : 0.0;

    for (int i = 1; i <= grid_.ncell; ++i) {
        real_t x_cell[3]; grid_.get_cell_center(i, x_cell);
        real_t xp = x_cell[0] * params::boxlen, yp = x_cell[1] * params::boxlen, zp = x_cell[2] * params::boxlen;
        int rm = 0;
        for (int r = 1; r < nreg; ++r) {
            real_t xc = (r < (int)xcs.size()) ? xcs[r] : 0.0, lx = (r < (int)lxs.size()) ? lxs[r] : 1e10;
            real_t yc = (r < (int)ycs.size()) ? ycs[r] : 0.0, ly = (r < (int)lys.size()) ? lys[r] : 1e10;
            real_t zc = (r < (int)zcs.size()) ? zcs[r] : 0.0, lz = (r < (int)lzs.size()) ? lzs[r] : 1e10;
            
            bool match = (std::abs(xp - xc) <= 0.5 * lx);
            if (NDIM > 1) match = match && (std::abs(yp - yc) <= 0.5 * ly);
            if (NDIM > 2) match = match && (std::abs(zp - zc) <= 0.5 * lz);
            
            if (match) rm = r;
        }

        real_t d = (rm < (int)drs.size()) ? drs[rm] : db, p = (rm < (int)prs.size()) ? prs[rm] : pb, u = (rm < (int)urs.size()) ? urs[rm] : ub, v = (rm < (int)vrs.size()) ? vrs[rm] : vb, w = (rm < (int)wrs.size()) ? wrs[rm] : wb;
        real_t A = (rm < (int)Ars.size()) ? Ars[rm] : Ab, B = (rm < (int)Brs.size()) ? Brs[rm] : Bb, C = (rm < (int)Crs.size()) ? Crs[rm] : Cb;

        grid_.uold(i, 1) = d; grid_.uold(i, 2) = d * u; grid_.uold(i, 3) = d * v; grid_.uold(i, 4) = d * w;
        grid_.uold(i, 5) = p / (gam - 1.0) + 0.5 * d * (u*u + v*v + w*w) + 0.5 * (A*A + B*B + C*C);

        int iv = 6;
#ifdef MHD
        grid_.uold(i, iv++) = A; grid_.uold(i, iv++) = B; grid_.uold(i, iv++) = C;
#endif
        for (int e = 0; e < nener; ++e) grid_.uold(i, iv++) = prads[rm][e];
#ifdef MHD
        grid_.uold(i, grid_.nvar - 2) = A; grid_.uold(i, grid_.nvar - 1) = B; grid_.uold(i, grid_.nvar) = C;
#endif
        grid_.cpu_map[i] = 1;
    }
    std::cout << "[Initializer] Applied ICs (nvar=" << grid_.nvar << ")." << std::endl;
}

} // namespace ramses
