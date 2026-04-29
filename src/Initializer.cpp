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
    int nvar_pure = grid_.nvar - 3;

    for (int i = 1; i <= grid_.ncell; ++i) {
        real_t x_cell[3];
        grid_.get_cell_center(i, x_cell);
        
        // Find level of this cell
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
        grid_.uold(i, nvar_pure + 1) = Bx_r;
        grid_.uold(i, nvar_pure + 2) = By_r;
        grid_.uold(i, nvar_pure + 3) = 0.0;
#endif
        grid_.cpu_map[i] = 1;
    }
    std::cout << "[Initializer] Applied Orszag-Tang ICs." << std::endl;
#ifdef MHD
    int ilevel = 1;
    real_t dx = params::boxlen / static_cast<real_t>(params::nx);
    real_t max_div = 0.0;
    int nvar_pure_all = grid_.nvar - 3;
    for (int i = 1; i <= grid_.ncell; ++i) {
        real_t div = (grid_.uold(i, nvar_pure_all + 1) - grid_.uold(i, 6));
        if (NDIM > 1) div += (grid_.uold(i, nvar_pure_all + 2) - grid_.uold(i, 7));
        if (NDIM > 2) div += (grid_.uold(i, nvar_pure_all + 3) - grid_.uold(i, 8));
        max_div = std::max(max_div, std::abs(div) / dx);
    }
    std::cout << "[Initializer] Max Div B = " << max_div << std::endl;
#endif
}

void Initializer::ana_disk_potential_condinit() {
    // Apply default region initialization first (optional, for other variables)
    region_condinit();

    std::string param_str = config_.get("poisson_params", "gravity_params", "");
    real_t a1 = 1.42e-3;
    real_t a2 = 5.49e-4;
    real_t z0 = 0.18e3;
    
    if (!param_str.empty()) {
        std::replace(param_str.begin(), param_str.end(), 'd', 'e');
        std::replace(param_str.begin(), param_str.end(), 'D', 'e');
        std::replace(param_str.begin(), param_str.end(), ',', ' ');
        std::stringstream ss(param_str);
        ss >> a1 >> a2 >> z0;
    }
    
    real_t scale_l = config_.get_double("units_params", "units_length", 1.0);
    real_t scale_t = config_.get_double("units_params", "units_time", 1.0);
    real_t scale_v = scale_l / scale_t;

    real_t kpc2cm = 3.085677581282e21;
    real_t pc2cm = 3.085677581282e18;
    real_t Myr2sec = 3.15576e13;

    // Convert gravity params to code units
    a1 = a1 * kpc2cm / (Myr2sec * Myr2sec) / scale_l * (scale_t * scale_t);
    a2 = a2 / (Myr2sec * Myr2sec) * (scale_t * scale_t);
    z0 = z0 * pc2cm / scale_l;

    real_t temp0 = 8000.0;
    real_t mu_gas = config_.get_double("cooling_params", "mu_gas", 1.4);
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);

    real_t kB = 1.3806490e-16;
    real_t mH = 1.6605390e-24;
    real_t cs2 = (kB * temp0 / (mu_gas * mH)) / (scale_v * scale_v);

    // Compute base density for equilibrium
    // rho(z) = rho0 / (1 + (z/z0)^2)^1.5 (for a1)
    real_t a1_rho = a1 / (4.0 * M_PI * cs2) * z0 * z0;
    real_t a2_rho = a2 / (2.0 * M_PI * cs2);

    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        for (int ilevel = 1; ilevel <= grid_.nlevelmax; ++ilevel) {
            int igrid = grid_.headl(icpu, ilevel);
            while (igrid > 0) {
                for (int ic = 1; ic <= constants::twotondim; ++ic) {
                    int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    real_t x = grid_.xg[(NDIM - 1) * grid_.ngridmax + igrid - 1];
                    int ix = (ic - 1) & 1;
                    int iy = ((ic - 1) & 2) >> 1;
                    int iz = ((ic - 1) & 4) >> 2;
                    real_t dx = 0.5 / static_cast<real_t>(1 << (ilevel - 1));
                    if (NDIM == 1) x += (static_cast<real_t>(ix) - 0.5) * dx;
                    else if (NDIM == 2) x += (static_cast<real_t>(iy) - 0.5) * dx;
                    else if (NDIM == 3) x += (static_cast<real_t>(iz) - 0.5) * dx;

                    real_t z_coord = (x - 0.5) * params::boxlen;
                    real_t rho = a1_rho / std::pow(1.0 + std::pow(z_coord / z0, 2), 1.5) + a2_rho;
                    
                    grid_.uold(ind_cell, 1) = rho;
                    // Set all velocities to 0
                    for (int idim = 1; idim <= NDIM; ++idim) {
                        grid_.uold(ind_cell, 1 + idim) = 0.0;
                    }
                    // Pressure via constant temperature
                    real_t p = rho * cs2;
                    grid_.uold(ind_cell, NDIM + 2) = p / (gamma - 1.0); // total energy (kin is 0)
                }
                igrid = grid_.next[igrid - 1];
            }
        }
    }
}

void Initializer::region_condinit() {
    int nregion = config_.get_int("init_params", "nregion", 1);
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    
    auto parse_list = [&](const std::string& key) {
        std::vector<real_t> vals;
        std::string s = config_.get("init_params", key, "");
        if (s.empty()) return vals;
        std::replace(s.begin(), s.end(), 'd', 'e');
        std::replace(s.begin(), s.end(), 'D', 'e');
        std::replace(s.begin(), s.end(), ',', ' ');
        std::stringstream ss(s);
        real_t v;
        while (ss >> v) vals.push_back(v);
        return vals;
    };

    auto d_regs = parse_list("d_region");
    auto p_regs = parse_list("p_region");
    auto u_regs = parse_list("u_region");
    auto v_regs = parse_list("v_region");
    auto w_regs = parse_list("w_region");
    auto A_regs = parse_list("A_region");
    auto B_regs = parse_list("B_region");
    auto C_regs = parse_list("C_region");
    auto x_centers = parse_list("x_center");
    auto length_xs = parse_list("length_x");

    // Default to region 1 for all cells initially
    real_t d_bg = !d_regs.empty() ? d_regs[0] : 1.0;
    real_t p_bg = !p_regs.empty() ? p_regs[0] : 1e-5;
    real_t u_bg = !u_regs.empty() ? u_regs[0] : 0.0;
    real_t v_bg = !v_regs.empty() ? v_regs[0] : 0.0;
    real_t w_bg = !w_regs.empty() ? w_regs[0] : 0.0;
    real_t A_bg = !A_regs.empty() ? A_regs[0] : 0.0;
    real_t B_bg = !B_regs.empty() ? B_regs[0] : 0.0;
    real_t C_bg = !C_regs.empty() ? C_regs[0] : 0.0;

    for (int i = 1; i <= grid_.ncell; ++i) {
        real_t x[3];
        grid_.get_cell_center(i, x);
        real_t x_phys = x[0] * params::boxlen;

        // Determine which region this cell belongs to (last matching region wins)
        int ireg_match = 0; // 0-based index
        for (int ireg = 1; ireg < nregion; ++ireg) {
            real_t xc = (ireg < (int)x_centers.size()) ? x_centers[ireg] : 0.0;
            real_t lx = (ireg < (int)length_xs.size()) ? length_xs[ireg] : 0.0;
            if (std::abs(x_phys - xc) <= 0.5 * lx) {
                ireg_match = ireg;
            }
        }

        real_t d = (ireg_match < (int)d_regs.size()) ? d_regs[ireg_match] : d_bg;
        real_t p = (ireg_match < (int)p_regs.size()) ? p_regs[ireg_match] : p_bg;
        real_t u = (ireg_match < (int)u_regs.size()) ? u_regs[ireg_match] : u_bg;
        real_t v = (ireg_match < (int)v_regs.size()) ? v_regs[ireg_match] : v_bg;
        real_t w = (ireg_match < (int)w_regs.size()) ? w_regs[ireg_match] : w_bg;
        real_t A = (ireg_match < (int)A_regs.size()) ? A_regs[ireg_match] : A_bg;
        real_t B = (ireg_match < (int)B_regs.size()) ? B_regs[ireg_match] : B_bg;
        real_t C = (ireg_match < (int)C_regs.size()) ? C_regs[ireg_match] : C_bg;

        grid_.uold(i, 1) = d;
        grid_.uold(i, 2) = d * u;
        grid_.uold(i, 3) = d * v;
        grid_.uold(i, 4) = d * w;
        
        real_t e_kin = 0.5 * d * (u*u + v*v + w*w);
        real_t e_mag = 0.5 * (A*A + B*B + C*C);
        grid_.uold(i, 5) = p / (gamma - 1.0) + e_kin + e_mag;

#ifdef MHD
        grid_.uold(i, 6) = A;
        grid_.uold(i, 7) = B;
        grid_.uold(i, 8) = C;
        grid_.uold(i, grid_.nvar - 2) = A;
        grid_.uold(i, grid_.nvar - 1) = B;
        grid_.uold(i, grid_.nvar) = C;
#endif
        grid_.cpu_map[i] = 1;
    }
    
    std::cout << "[Initializer] Applied ICs (nvar=" << grid_.nvar << ")." << std::endl;
}

} // namespace ramses
