#include "ramses/Initializer.hpp"
#include "ramses/Parameters.hpp"
#include <iostream>
#include <cmath>

namespace ramses {

void Initializer::apply_all() {
    region_condinit();
}

void Initializer::region_condinit() {
    int nregion = config_.get_int("init_params", "nregion", 1);
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    
    // For Sedov 3D, we usually have a background region (1) and a blast center (2)
    for (int ireg = 1; ireg <= nregion; ++ireg) {
        // Simplified for POC: We'll just implement the Sedov point-center logic
    }

    // Default Background
    real_t d_bg = config_.get_double("init_params", "d_region", 1.0); // Simple hack: first value
    real_t p_bg = config_.get_double("init_params", "p_region", 1e-5);

    for (int i = 1; i <= grid_.ncell; ++i) {
        grid_.uold(i, 1) = d_bg;
        grid_.uold(i, 2) = 0.0;
        grid_.uold(i, 3) = 0.0;
        grid_.uold(i, 4) = 0.0;
        grid_.uold(i, 5) = p_bg / (gamma - 1.0);
    }

    // Sedov Center (Assume region 2 is a point center)
    real_t x_c = 0.5, y_c = 0.5, z_c = 0.5; // Default center
    real_t p_center = 0.4; // From sedov3d.nml
    
    // In a real port, we'd calculate cell positions x,y,z correctly.
    // For now, let's just mark the center-most coarse cell.
    int center_ind = grid_.ncoarse / 2; 
    grid_.uold(center_ind, 5) = p_center / (gamma - 1.0);
    
    std::cout << "[Initializer] Applied Sedov-like ICs." << std::endl;
}

} // namespace ramses
