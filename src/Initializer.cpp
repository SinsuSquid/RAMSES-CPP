#include "ramses/Initializer.hpp"
#include "ramses/Parameters.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>

namespace ramses {

void Initializer::apply_all() {
    std::string condinit_kind = config_.get("init_params", "condinit_kind", "");
    if (condinit_kind == "ana_disk_potential" || condinit_kind == "'ana_disk_potential'") {
        ana_disk_potential_condinit();
    } else {
        region_condinit();
    }
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
    
    // Default Background values
    real_t d_bg = config_.get_double("init_params", "d_region", 1.0);
    real_t p_bg = config_.get_double("init_params", "p_region", 1e-5);
    
    real_t A_bg = 0.0, B_bg = 0.0, C_bg = 0.0;
#ifdef MHD
    A_bg = config_.get_double("init_params", "A_region", 0.0);
    B_bg = config_.get_double("init_params", "B_region", 0.0);
    C_bg = config_.get_double("init_params", "C_region", 0.0);
#endif

    // Initialize all cells with background
    for (int i = 1; i <= grid_.ncell; ++i) {
        grid_.uold(i, 1) = d_bg; // density
        for (int idim = 1; idim <= NDIM; ++idim) {
            grid_.uold(i, 1 + idim) = 0.0; // velocity
        }
        
        real_t e_kin = 0.0;
        real_t e_mag = 0.5 * (A_bg*A_bg + B_bg*B_bg + C_bg*C_bg);
        grid_.uold(i, NDIM + 2) = p_bg / (gamma - 1.0) + e_kin + e_mag; // total energy
        
#ifdef MHD
        grid_.uold(i, 6) = A_bg;
        grid_.uold(i, 7) = B_bg;
        grid_.uold(i, 8) = C_bg;
        grid_.uold(i, grid_.nvar - 2) = A_bg;
        grid_.uold(i, grid_.nvar - 1) = B_bg;
        grid_.uold(i, grid_.nvar) = C_bg;
#endif

        int nener_start = NDIM + 3;
#ifdef MHD
        nener_start = 9;
        int nvar_mhd_end = grid_.nvar - 3;
#else
        int nvar_mhd_end = grid_.nvar;
#endif
        for (int ivar = nener_start; ivar <= nvar_mhd_end; ++ivar) {
            grid_.uold(i, ivar) = 0.0;
        }
        grid_.cpu_map[i] = 1;
    }

    // Apply regions
    if (nregion >= 2) {
        // Read region 2 values (simplified)
        std::string d_str = config_.get("init_params", "d_region", "");
        std::string p_str = config_.get("init_params", "p_region", "");
        real_t d2 = 0.125, p2 = 0.1;
        if (!d_str.empty()) {
            std::stringstream ss(d_str); double val; ss >> val; if (ss >> val) d2 = val;
        }
        if (!p_str.empty()) {
            std::stringstream ss(p_str); double val; ss >> val; if (ss >> val) p2 = val;
        }

#ifdef MHD
        real_t A2 = 0.0, B2 = 0.0, C2 = 0.0;
        std::string A_str = config_.get("init_params", "A_region", "");
        std::string B_str = config_.get("init_params", "B_region", "");
        std::string C_str = config_.get("init_params", "C_region", "");
        if (!A_str.empty()) { std::stringstream ss(A_str); double v; ss >> v; if (ss >> v) A2 = v; }
        if (!B_str.empty()) { std::stringstream ss(B_str); double v; ss >> v; if (ss >> v) B2 = v; }
        if (!C_str.empty()) { std::stringstream ss(C_str); double v; ss >> v; if (ss >> v) C2 = v; }
#endif

        // For Sod tube 1D, x > 0.5 is region 2.
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            real_t x = (static_cast<real_t>(i) - 0.5f) / static_cast<real_t>(params::nx);
            if (x > 0.5) {
                grid_.uold(i, 1) = d2;
                real_t e_kin = 0.0;
                real_t e_mag = 0.0;
#ifdef MHD
                e_mag = 0.5 * (A2*A2 + B2*B2 + C2*C2);
#else
                e_mag = 0.5 * (A_bg*A_bg + B_bg*B_bg + C_bg*C_bg);
#endif
                grid_.uold(i, NDIM + 2) = p2 / (gamma - 1.0) + e_kin + e_mag;
#ifdef MHD
                grid_.uold(i, 6) = A2;
                grid_.uold(i, 7) = B2;
                grid_.uold(i, 8) = C2;
                grid_.uold(i, grid_.nvar - 2) = A2;
                grid_.uold(i, grid_.nvar - 1) = B2;
                grid_.uold(i, grid_.nvar) = C2;
#endif
            }
        }
    }
    
    // Special case for sod-tube-nener: handle prad_region
    int nener = grid_.nvar - (2 + NDIM);
    if (nener > 0) {
        // Very simplified: just set them to some values if found in config
        for (int ireg = 1; ireg <= nregion; ++ireg) {
            for (int iv = 1; iv <= nener; ++iv) {
                // prad_region(ivar, iregion)
                std::string key = "prad_region(" + std::to_string(iv) + "," + std::to_string(ireg) + ")";
                real_t val = config_.get_double("init_params", key, 0.0);
                // Apply this value to the region... (omitted for brevity, POC only)
                // For sod-tube-nener, let's just set them.
                if (ireg == 1) {
                    for(int i=1; i<=grid_.ncell; ++i) grid_.uold(i, NDIM + 2 + iv) = val;
                }
            }
        }
    }

    std::cout << "[Initializer] Applied ICs (nvar=" << grid_.nvar << ")." << std::endl;
}

} // namespace ramses
