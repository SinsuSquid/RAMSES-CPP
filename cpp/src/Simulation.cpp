#include "ramses/Simulation.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Initializer.hpp"
#include <iostream>

namespace ramses {

void Simulation::initialize(const std::string& nml_path) {
    std::cout << "[Simulation] Initializing from " << nml_path << "..." << std::endl;
    if (!config_.parse(nml_path)) {
        return;
    }

    // Map parameters
    namespace p = ramses::params;
    p::nx = config_.get_int("amr_params", "nx", 2);
    p::ny = config_.get_int("amr_params", "ny", 2);
    p::nz = config_.get_int("amr_params", "nz", 2);
    int ngridmax = config_.get_int("amr_params", "ngridmax", 1000);
    int nvar = config_.get_int("hydro_params", "nvar", 5);
    int nlevelmax = config_.get_int("amr_params", "levelmax", 10);
    nstepmax_ = config_.get_int("run_params", "nstepmax", 10);

    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar, 1, nlevelmax);
    
    // Apply Initial Conditions
    Initializer init(grid_, config_);
    init.apply_all();
    
    std::cout << "[Simulation] Grid allocated and initialized." << std::endl;
}

void Simulation::run() {
    std::cout << "[Simulation] Starting main loop..." << std::endl;
    
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    real_t courant_factor = config_.get_double("hydro_params", "courant_factor", 0.8);
    real_t dx = config_.get_double("amr_params", "boxlen", 1.0) / static_cast<real_t>(params::nx);

    while (t_ < tend_ && nstep_ < nstepmax_) {
        nstep_++;
        
        // Calculate dynamic timestep
        real_t dt = hydro_.compute_courant_step(1, dx, gamma, courant_factor);
        if (dt > 0.1) dt = 0.1; // Cap for stability in initial steps

        std::cout << "  Step " << nstep_ << " t=" << t_ << " dt=" << dt << std::endl;
        
        // Root level step
        amr_step(1);
        
        t_ += dt;
    }
    
    std::cout << "[Simulation] Run complete." << std::endl;
}

void Simulation::amr_step(int ilevel) {
    if (ilevel > grid_.nlevelmax) return;

    // 1. Poisson solver (gravity)
    poisson_.solve(ilevel);

    // 2. Hydro update for this level
    hydro_.godunov_fine(ilevel);
    
    // 2. Recursive step for finer levels (Sub-cycling)
    if (ilevel < grid_.nlevelmax) {
        amr_step(ilevel + 1);
        amr_step(ilevel + 1);
    }
    
    // 3. Restriction (average fine cells to coarse)
    // This is where upload_fine(ilevel) would go
}


} // namespace ramses
