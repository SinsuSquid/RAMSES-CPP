#include "ramses/Simulation.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Initializer.hpp"
#include "ramses/RamsesWriter.hpp"
#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <sstream>

namespace ramses {

void Simulation::initialize(const std::string& nml_path) {
    std::cout << "[Simulation] Initializing from " << nml_path << "..." << std::endl;
    if (!config_.parse(nml_path)) {
        return;
    }

    // Map parameters
    namespace p = ramses::params;
    p::nx = config_.get_int("amr_params", "nx", 1);
    p::ny = config_.get_int("amr_params", "ny", 1);
    p::nz = config_.get_int("amr_params", "nz", 1);
    int ngridmax = config_.get_int("amr_params", "ngridmax", 0);
    if (ngridmax == 0) ngridmax = config_.get_int("amr_params", "ngridtot", 1000000);
    int nvar = config_.get_int("hydro_params", "nvar", 5);
    int nlevelmax = config_.get_int("amr_params", "levelmax", 10);
    if (nlevelmax < 20) nlevelmax = 20; // Ensure headroom for ilevel+1 access
    nstepmax_ = config_.get_int("run_params", "nstepmax", 10);

    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar, 1, nlevelmax);
    
    // Build initial tree up to levelmin
    int levelmin = config_.get_int("amr_params", "levelmin", 1);
    std::cout << "[Simulation] Building initial AMR tree up to level " << levelmin << "..." << std::endl;
    
    // Level 1 refinement
    updater_.mark_all(1);
    updater_.refine_coarse();
    
    // Further levels
    for (int ilevel = 1; ilevel < levelmin; ++ilevel) {
        std::cout << "  [Simulation] Initial Refinement Level " << ilevel << " -> " << ilevel + 1 << std::endl;
        updater_.mark_all(ilevel);
        updater_.refine_fine(ilevel);
    }

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
        real_t dt = hydro_.compute_courant_step(1, dx, gamma, courant_factor);
        if (dt > 0.1) dt = 0.1;
        std::cout << "  Step " << nstep_ << " t=" << t_ << " dt=" << dt << std::endl;
        amr_step(1);
        t_ += dt;
    }
    
    std::cout << "[Simulation] Run complete. Dumping final snapshot..." << std::endl;
    dump_snapshot(1);
}

void Simulation::dump_snapshot(int iout) {
    std::stringstream ss;
    ss << "output_" << std::setw(5) << std::setfill('0') << iout;
    std::string dir = ss.str();
    mkdir(dir.c_str(), 0777);
    
    // 1. AMR file
    std::string amr_file = dir + "/amr_00001.out00001"; // Match CPU numbering
    {
        RamsesWriter writer(amr_file);
        if (writer.is_open()) writer.write_amr(grid_);
    }

    // 2. Hydro file
    std::string hydro_file = dir + "/hydro_00001.out00001";
    {
        RamsesWriter writer(hydro_file);
        if (writer.is_open()) writer.write_hydro(grid_, grid_.nlevelmax);
    }

    // 3. Info file (Plain text)
    std::string info_file = dir + "/info_" + std::string(ss.str().substr(7)) + ".txt";
    std::ofstream info(info_file);
    if (info.is_open()) {
        info << "ncpu         = " << grid_.ncpu << "\n";
        info << "ndim         = " << NDIM << "\n";
        info << "nx           = " << params::nx << "\n";
        info << "ny           = " << params::ny << "\n";
        info << "nz           = " << params::nz << "\n";
        info << "levelmin     = " << config_.get_int("amr_params", "levelmin", 1) << "\n";
        info << "levelmax     = " << grid_.nlevelmax << "\n";
        info << "ngridmax     = " << grid_.ngridmax << "\n";
        info << "boxlen       = 1.0\n";
        info << "time         = " << t_ << "\n";
        info << "unit_l       = 1.0\n";
        info << "unit_d       = 1.0\n";
        info << "unit_t       = 1.0\n";
        info.close();
    }

    // 4. Descriptor file
    std::ofstream desc(dir + "/hydro_file_descriptor.txt");
    if (desc.is_open()) {
        desc << "# version:  1\n";
        desc << "1, density, double\n";
        desc << "2, velocity_x, double\n";
        if (NDIM > 1) desc << "3, velocity_y, double\n";
        if (NDIM > 2) desc << "4, velocity_z, double\n";
        desc << "5, pressure, double\n";
        desc.close();
    }
    
    std::cout << "[Simulation] Snapshot " << iout << " dumped successfully." << std::endl;
}

void Simulation::amr_step(int ilevel) {
    if (ilevel > grid_.nlevelmax) return;
    if (grid_.count_grids_at_level(ilevel) == 0) return;

    // 1. Poisson solver (gravity)
    poisson_.solve(ilevel);

    // 2. Hydro update for this level
    hydro_.godunov_fine(ilevel);
    
    // 3. Recursive step for finer levels (Sub-cycling)
    if (ilevel < grid_.nlevelmax) {
        if (grid_.count_grids_at_level(ilevel + 1) > 0) {
            amr_step(ilevel + 1);
            amr_step(ilevel + 1);
        }
    }
}

} // namespace ramses
