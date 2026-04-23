#include "ramses/Simulation.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Initializer.hpp"
#include "ramses/RamsesWriter.hpp"
#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <sstream>
#include <algorithm>

namespace ramses {

void Simulation::initialize(const std::string& nml_path) {
    if (!config_.parse(nml_path)) return;
    
    namespace p = ramses::params;
    p::nx = config_.get_int("amr_params", "nx", 1);
    p::ny = config_.get_int("amr_params", "ny", 1);
    p::nz = config_.get_int("amr_params", "nz", 1);
    
    int ngridmax = config_.get_int("amr_params", "ngridmax", 0);
    if (ngridmax == 0) ngridmax = config_.get_int("amr_params", "ngridtot", 1000000);
    int nvar = config_.get_int("hydro_params", "nvar", 5);
    int nlevelmax = config_.get_int("amr_params", "levelmax", 10);
    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar, 1, nlevelmax);
    
    // Load Outputs
    noutput_ = config_.get_int("output_params", "noutput", 0);
    std::string tout_str = config_.get("output_params", "tout", "");
    if (!tout_str.empty()) {
        std::stringstream ss(tout_str);
        double val;
        while (ss >> val) {
            tout_.push_back(val);
            if (ss.peek() == ',' || ss.peek() == ' ') ss.ignore();
        }
    }
    tend_ = config_.get_double("output_params", "tend", 1e10);
    if (!tout_.empty()) tend_ = tout_.back();

    // Initial AMR tree
    int levelmin = config_.get_int("amr_params", "levelmin", 1);
    updater_.mark_all(1);
    updater_.refine_coarse();
    for (int ilevel = 1; ilevel < levelmin; ++ilevel) {
        updater_.mark_all(ilevel);
        updater_.refine_fine(ilevel);
    }
    
    Initializer init(grid_, config_);
    init.apply_all();
}

void Simulation::run() {
    std::cout << "[Simulation] Starting main loop..." << std::endl;
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    real_t courant_factor = config_.get_double("hydro_params", "courant_factor", 0.8);
    
    dump_snapshot(1);
    int iout = 0;
    int snapshot_count = 2;

    while (t_ < tend_) {
        nstep_++;
        
        // Compute adaptive dt for level 1
        real_t dx = params::boxlen / static_cast<real_t>(params::nx);
        real_t dt = hydro_.compute_courant_step(1, dx, gamma, courant_factor);
        
        // Check if we reached next tout
        if (iout < (int)tout_.size()) {
            if (t_ + dt >= tout_[iout]) {
                dt = tout_[iout] - t_;
                amr_step(1, dt);
                t_ = tout_[iout];
                dump_snapshot(snapshot_count++);
                iout++;
                continue;
            }
        } else if (t_ + dt >= tend_) {
            dt = tend_ - t_;
            amr_step(1, dt);
            t_ = tend_;
            dump_snapshot(snapshot_count++);
            break;
        }

        amr_step(1, dt);
        t_ += dt;
        
        if (nstep_ % 10 == 0) {
            std::cout << "  Step " << nstep_ << " t=" << t_ << " dt=" << dt << std::endl;
        }
    }
}

void Simulation::amr_step(int ilevel, real_t dt) {
    if (ilevel > grid_.nlevelmax) return;
    if (grid_.count_grids_at_level(ilevel) == 0) return;

    // TODO: Pass dt to solvers for physical scaling
    hydro_.godunov_fine(ilevel);

    // Sub-cycling
    if (ilevel < grid_.nlevelmax) {
        if (grid_.count_grids_at_level(ilevel + 1) > 0) {
            amr_step(ilevel + 1, dt / 2.0);
            amr_step(ilevel + 1, dt / 2.0);
        }
    }
}

void Simulation::dump_snapshot(int iout) {
    std::stringstream ss;
    ss << "output_" << std::setw(5) << std::setfill('0') << iout;
    std::string dir = ss.str();
    mkdir(dir.c_str(), 0777);
    std::string nchar = ss.str().substr(7);
    
    std::string amr_path = dir + "/amr_" + nchar + ".out00001";
    {
        RamsesWriter writer(amr_path);
        if (writer.is_open()) writer.write_amr(grid_);
    }

    std::string hydro_path = dir + "/hydro_" + nchar + ".out00001";
    {
        RamsesWriter writer(hydro_path);
        if (writer.is_open()) writer.write_hydro(grid_, grid_.nlevelmax);
    }

    std::string info_file = dir + "/info_" + nchar + ".txt";
    std::ofstream info(info_file);
    if (info.is_open()) {
        info << "ncpu         = 1\n";
        info << "ndim         = " << NDIM << "\n";
        info << "nx           = " << params::nx << "\n";
        info << "ny           = " << params::ny << "\n";
        info << "nz           = " << params::nz << "\n";
        info << "levelmin     = 1\n";
        info << "levelmax     = " << grid_.nlevelmax << "\n";
        info << "ngridmax     = " << grid_.ngridmax << "\n";
        info << "boxlen       = " << params::boxlen << "\n";
        info << "time         = " << t_ << "\n";
        info << "unit_l       = " << config_.get_double("units_params", "units_length", 1.0) << "\n";
        info << "unit_d       = " << config_.get_double("units_params", "units_density", 1.0) << "\n";
        info << "unit_t       = " << config_.get_double("units_params", "units_time", 1.0) << "\n";
        info.close();
    }

    std::ofstream desc(dir + "/hydro_file_descriptor.txt");
    desc << "1, density, double\n2, velocity_x, double\n3, velocity_y, double\n4, velocity_z, double\n5, pressure, double\n";
    
    std::ofstream header(dir + "/header_" + nchar + ".txt");
    header << "total 0\nlost 0\n";
    
    std::cout << "[Simulation] Snapshot " << iout << " dumped (t=" << t_ << ")." << std::endl;
}

} // namespace ramses
