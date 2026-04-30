#include "ramses/Simulation.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Initializer.hpp"
#include "ramses/RamsesWriter.hpp"
#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <sstream>
#include <algorithm>
#include <cmath>

namespace ramses {

void Simulation::initialize(const std::string& nml_path) {
    if (!config_.parse(nml_path)) return;
    
    namespace p = ramses::params;
    p::nx = config_.get_int("amr_params", "nx", 1);
    p::ny = config_.get_int("amr_params", "ny", 1);
    p::nz = config_.get_int("amr_params", "nz", 1);
    p::boxlen = config_.get_double("amr_params", "boxlen", 1.0);
    
    int ngridmax = config_.get_int("amr_params", "ngridmax", 0);
    if (ngridmax == 0) ngridmax = config_.get_int("amr_params", "ngridtot", 1000000);
    
    nener_ = config_.get_int("hydro_params", "nener", 0);
    // If nener is 0, check if prad_region is present in init_params
    if (nener_ == 0) {
        for (int i = 1; i <= 10; ++i) {
            std::string key = "prad_region(" + std::to_string(i) + ",1)";
            if (config_.get( "init_params", key, "").length() > 0) {
                nener_ = std::max(nener_, i);
            }
        }
    }
    // Standard RAMSES: nvar = 2 + NDIM + nener
    int nvar_default = 2 + NDIM + nener_;
#ifdef MHD
    nvar_default = 5 + 3 + nener_; // nhydro=8 (d,u,v,w,p,B_left)
#endif
#ifdef RT
    int nGroups = config_.get_int("rt_params", "nGroups", 0);
    nvar_default += nGroups * (1 + NDIM);
#endif
    int nvar = config_.get_int("hydro_params", "nvar", nvar_default);
    
#ifdef MHD
    int nvar_all = nvar + 3; // extra 3 for B_right
    int nlevelmax = config_.get_int("amr_params", "levelmax", 10);
    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar_all, 1, nlevelmax);
#else
    int nlevelmax = config_.get_int("amr_params", "levelmax", 10);
    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar, 1, nlevelmax);
#endif
    
    // Load Outputs
    noutput_ = config_.get_int("output_params", "noutput", 0);
    ncontrol_ = config_.get_int("run_params", "ncontrol", 1);

    // Units
    p::units_length = config_.get_double("units_params", "units_length", 1.0);
    p::units_density = config_.get_double("units_params", "units_density", 1.0);
    p::units_time = config_.get_double("units_params", "units_time", 1.0);
    p::units_velocity = p::units_length / p::units_time;
    p::units_mass = p::units_density * std::pow(p::units_length, 3);
    p::units_energy = p::units_mass * std::pow(p::units_velocity, 2);
    p::units_pressure = p::units_density * std::pow(p::units_velocity, 2);
    p::units_number_density = p::units_density / (1.67e-24); // approx mH
    
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
    updater_.mark_all(1); // Marks level 0 cells
    updater_.refine_coarse(); // Refines level 0 -> creates level 1 grids
    for (int ilevel = 1; ilevel < levelmin; ++ilevel) {
        updater_.mark_all(ilevel + 1); // Marks level ilevel cells
        updater_.refine_fine(ilevel); // Refines level ilevel -> creates level ilevel+1 grids
    }
    
    Initializer init(grid_, config_);
    init.apply_all();

    // Initial load balance
    balancer_.calculate_hilbert_keys();
    balancer_.balance();

    // RT initialization
    rt_.initialize();
}

void Simulation::run() {
    std::cout << "[Simulation] Starting main loop..." << std::endl;
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    real_t courant_factor = config_.get_double("hydro_params", "courant_factor", 0.8);
    
    if (config_.get_bool("run_params", "poisson", false)) {
        for (int ilevel = grid_.nlevelmax; ilevel >= 1; --ilevel) {
            poisson_.solve(ilevel);
            poisson_.compute_force(ilevel);
        }
    }

    dump_snapshot(1);
    int iout = 0;
    int snapshot_count = 2;

    while (t_ < tend_) {
        // Output Mesh Structure and Main Step info every ncontrol steps
        if (nstep_ % ncontrol_ == 0) {
            std::cout << "Mesh structure" << std::endl;
            for (int i = 1; i <= grid_.nlevelmax; ++i) {
                int ngrids = grid_.count_grids_at_level(i);
                if (ngrids > 0) {
                    printf(" Level %2d has %10d grids\n", i, ngrids);
                }
            }
            // Simplified conservation output
            printf(" Main step=%7d t=%12.5e dt=%10.3e\n", nstep_, t_, 0.0); // dt is for fine steps
        }

        nstep_++;
        
        // Compute adaptive dt for level 1
        real_t dx = params::boxlen / static_cast<real_t>(params::nx);
        real_t dt;
#ifdef MHD
        dt = mhd_.compute_courant_step(1, dx, gamma, courant_factor);
#else
        dt = hydro_.compute_courant_step(1, dx, gamma, courant_factor);
#endif
        
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
            std::cout << "  Step " << nstep_ << " t=" << t_ << " dt=" << dt;
#ifdef MHD
            real_t dx = params::boxlen / static_cast<real_t>(params::nx);
            std::cout << " max_div_b=" << mhd_.compute_max_div_b(1, dx);
#endif
            std::cout << std::endl;
        }
    }
}

void Simulation::amr_step(int ilevel, real_t dt) {
    if (ilevel > grid_.nlevelmax) return;
    if (grid_.count_grids_at_level(ilevel) == 0) return;

#ifdef MHD
    mhd_.set_unew(ilevel);
#else
    hydro_.set_unew(ilevel);
#endif
    rt_.set_unew(ilevel);

    // Sub-cycling
    if (ilevel < grid_.nlevelmax) {
        if (grid_.count_grids_at_level(ilevel + 1) > 0) {
            amr_step(ilevel + 1, dt / 2.0);
            amr_step(ilevel + 1, dt / 2.0);
        }
    }

    real_t dx = params::boxlen / static_cast<real_t>(params::nx * (1 << (ilevel - 1)));
    
    if (config_.get_bool("run_params", "poisson", false)) {
        poisson_.solve(ilevel);
        poisson_.compute_force(ilevel);
        hydro_.add_gravity_source_terms(ilevel, dt);
    }
    
#ifdef MHD
    mhd_.godunov_fine(ilevel, dt, dx);
#else
    hydro_.godunov_fine(ilevel, dt, dx);
#endif

    // Cooling
    cooling_.apply_cooling(ilevel, dt);

    // RT
    rt_.godunov_fine(ilevel, dt, dx);

    // Diagnostics
    real_t mind = 1e30, maxv = 0, maxdivb = 0;
#ifdef MHD
    mhd_.get_diagnostics(ilevel, dx, mind, maxv, maxdivb);
    printf(" Fine step=%7d t=%12.5e dt=%10.3e level=%2d min_d=%9.2e max_v=%9.2e max_div_b=%9.2e\n", 
           nstep_, t_ + dt, dt, ilevel, mind, maxv, maxdivb);
#else
    hydro_.get_diagnostics(ilevel, dx, mind, maxv);
    printf(" Fine step=%7d t=%12.5e dt=%10.3e level=%2d min_d=%9.2e max_v=%9.2e\n", 
           nstep_, t_ + dt, dt, ilevel, mind, maxv);
#endif

    if (config_.get_bool("run_params", "poisson", false)) {
        hydro_.add_gravity_source_terms(ilevel, dt);
    }

#ifdef MHD
    mhd_.set_uold(ilevel);
#else
    hydro_.set_uold(ilevel);
#endif
    rt_.set_uold(ilevel);

    real_t min_d = 1e30, max_v = 0.0, max_div_b = 0.0;
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                real_t d = grid_.uold(ind_cell, 1);
                min_d = std::min(min_d, d);
                real_t v2 = 0;
                for (int i=1; i<=NDIM; ++i) v2 += std::pow(grid_.uold(ind_cell, 1+i)/d, 2);
                max_v = std::max(max_v, std::sqrt(v2));
#ifdef MHD
                int nvar_pure = grid_.nvar - 3;
                real_t div = (grid_.uold(ind_cell, nvar_pure + 1) - grid_.uold(ind_cell, 6));
                if (NDIM > 1) div += (grid_.uold(ind_cell, nvar_pure + 2) - grid_.uold(ind_cell, 7));
                if (NDIM > 2) div += (grid_.uold(ind_cell, nvar_pure + 3) - grid_.uold(ind_cell, 8));
                max_div_b = std::max(max_div_b, std::abs(div) / dx);
#endif
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    std::cout << "  amr_step(" << ilevel << ", dt=" << dt << "): min_d=" << min_d << ", max_v=" << max_v;
#ifdef MHD
    std::cout << ", max_div_b=" << max_div_b;
#endif
    std::cout << std::endl;
}

void Simulation::dump_snapshot(int iout) {
    std::stringstream ss;
    ss << "output_" << std::setw(5) << std::setfill('0') << iout;
    std::string dir = ss.str();
    mkdir(dir.c_str(), 0777);
    std::string nchar = ss.str().substr(7);
    
    SnapshotInfo info;
    info.t = t_;
    info.nstep = nstep_;
    info.noutput = noutput_;
    info.iout = iout;
    info.tout = tout_;
    info.gamma = config_.get_double("hydro_params", "gamma", 1.4);

    std::string amr_path = dir + "/amr_" + nchar + ".out00001";
    {
        RamsesWriter writer(amr_path);
        if (writer.is_open()) writer.write_amr(grid_, info);
    }

    std::string hydro_path = dir + "/hydro_" + nchar + ".out00001";
    {
        RamsesWriter writer(hydro_path);
        if (writer.is_open()) writer.write_hydro(grid_, info);
    }

    if (config_.get_bool("run_params", "poisson", false)) {
        std::string grav_path = dir + "/grav_" + nchar + ".out00001";
        RamsesWriter writer(grav_path);
        if (writer.is_open()) writer.write_grav(grid_, info);
    }

    std::string info_file = dir + "/info_" + nchar + ".txt";
    std::ofstream infof(info_file);
    if (infof.is_open()) {
        infof << std::scientific << std::setprecision(15);
        infof << "ncpu         = " << grid_.ncpu << "\n";
        infof << "ndim         = " << NDIM << "\n";
        infof << "nx           = " << params::nx << "\n";
        infof << "ny           = " << params::ny << "\n";
        infof << "nz           = " << params::nz << "\n";
        infof << "levelmin     = 1\n";
        infof << "levelmax     = " << grid_.nlevelmax << "\n";
        infof << "ngridmax     = " << grid_.ngridmax << "\n";
        infof << "boxlen       = " << params::boxlen << "\n";
        infof << "time         = " << t_ << "\n";
        infof << "unit_l       = " << config_.get_double("units_params", "units_length", 1.0) << "\n";
        infof << "unit_d       = " << config_.get_double("units_params", "units_density", 1.0) << "\n";
        infof << "unit_t       = " << config_.get_double("units_params", "units_time", 1.0) << "\n";
        infof.close();
    }

    std::ofstream desc(dir + "/hydro_file_descriptor.txt");
    if (desc.is_open()) {
        int ivar = 1;
        desc << ivar++ << ", density, double\n";
#ifdef MHD
        desc << ivar++ << ", velocity_x, double\n";
        desc << ivar++ << ", velocity_y, double\n";
        desc << ivar++ << ", velocity_z, double\n";
        desc << ivar++ << ", B_x_left, double\n";
        desc << ivar++ << ", B_y_left, double\n";
        desc << ivar++ << ", B_z_left, double\n";
        desc << ivar++ << ", B_x_right, double\n";
        desc << ivar++ << ", B_y_right, double\n";
        desc << ivar++ << ", B_z_right, double\n";
#else
        for (int i = 1; i <= NDIM; ++i) {
            char dim_char = (i == 1) ? 'x' : (i == 2 ? 'y' : 'z');
            desc << ivar++ << ", velocity_" << dim_char << ", double\n";
        }
#endif
        desc << ivar++ << ", pressure, double\n";
        
        for (int i = 1; i <= nener_; ++i) {
            desc << ivar++ << ", non_thermal_pressure_" << std::setw(2) << std::setfill('0') << i << ", double\n";
        }
#ifndef MHD
        int nvar_logical = 2 + NDIM + nener_;
        int nvar_file = config_.get_int("hydro_params", "nvar", nvar_logical);
        for (int i = nvar_logical + 1; i <= nvar_file; ++i) {
            desc << ivar++ << ", scalar_" << std::setw(2) << std::setfill('0') << i - nvar_logical << ", double\n";
        }
#else
        int nhydro = 8;
        int nvar_logical = nhydro + nener_;
        int nvar_file = config_.get_int("hydro_params", "nvar", nvar_logical);
        for (int i = nvar_logical + 1; i <= nvar_file; ++i) {
            desc << ivar++ << ", scalar_" << std::setw(2) << std::setfill('0') << i - nvar_logical << ", double\n";
        }
#endif
        desc.close();
    }
    
    std::ofstream header(dir + "/header_" + nchar + ".txt");
    header << "total 0\nlost 0\n";
    
    std::cout << "[Simulation] Snapshot " << iout << " dumped (t=" << t_ << ")." << std::endl;
}

} // namespace ramses
