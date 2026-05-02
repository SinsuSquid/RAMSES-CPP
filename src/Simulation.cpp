#include "ramses/Simulation.hpp"
#include "ramses/Initializer.hpp"
#include "ramses/RamsesWriter.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/MpiManager.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <algorithm>

#ifdef RAMSES_USE_MPI
#include <mpi.h>
#endif

namespace ramses {
namespace p = params;

void Simulation::initialize(const std::string& nml_path) {
    config_.parse(nml_path);

    p::nx = config_.get_int("amr_params", "nx", 1);
    p::ny = config_.get_int("amr_params", "ny", 1);
    p::nz = config_.get_int("amr_params", "nz", 1);
    p::boxlen = config_.get_double("amr_params", "boxlen", 1.0);

    int ngridmax = config_.get_int("amr_params", "ngridmax", 0);
    if (ngridmax == 0) ngridmax = config_.get_int("amr_params", "ngridtot", 1000000);
    p::ngridmax = ngridmax;

    nener_ = config_.get_int("hydro_params", "nener", 0);
    int nvar_default = 5 + nener_;
#ifdef MHD
    nvar_default = 8 + nener_;
#endif
#ifdef RT
    int nGroups = config_.get_int("rt_params", "nGroups", 0);
    nvar_default += nGroups * (1 + NDIM);
#endif
    int nvar = config_.get_int("hydro_params", "nvar", nvar_default);

    int levelmin = config_.get_int("amr_params", "levelmin", 1);
    int levelmax = config_.get_int("amr_params", "levelmax", 1);
    params::levelmin = levelmin;
    params::nlevelmax = levelmax;
    int ncpu = MpiManager::instance().size();

    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    grid_.gamma = gamma;

    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar, ncpu, levelmax);

    // Initial refinement pass (up to levelmin)
    for (int il = 1; il < levelmin; ++il) {
        initializer_.apply_all();
        updater_.mark_cells(il);
        if (il == 1) updater_.refine_coarse();
        else updater_.refine_fine(il);
    }
    initializer_.apply_all(); // Final population of leaf cells

    // Parse boundary_params
    grid_.nboundary = config_.get_int("boundary_params", "nboundary", 0);
    if (grid_.nboundary > 0) {
        std::string bmin_s = config_.get("boundary_params", "ibound_min", "");
        std::string bmax_s = config_.get("boundary_params", "ibound_max", "");
        std::string btype_s = config_.get("boundary_params", "bound_type", "");
        
        auto parse_vec = [](const std::string& s) {
            std::vector<int> v;
            std::string sc = s;
            std::replace(sc.begin(), sc.end(), ',', ' ');
            std::stringstream ss(sc); int val;
            while (ss >> val) v.push_back(val);
            return v;
        };
        grid_.ibound_min = parse_vec(bmin_s);
        grid_.ibound_max = parse_vec(bmax_s);
        grid_.bound_type = parse_vec(btype_s);
    }

    t_ = 0.0;
    tend_ = config_.get_double("output_params", "tend", 1.0);
    nstepmax_ = config_.get_int("run_params", "nstepmax", 1000000);
    ncontrol_ = config_.get_int("run_params", "ncontrol", 1);
    if (ncontrol_ > 100) ncontrol_ = 100; // Prevent timeouts
    
    // Parse tout
    std::string tout_s = config_.get("output_params", "tout", "");
    if (!tout_s.empty()) {
        std::replace(tout_s.begin(), tout_s.end(), ',', ' ');
        std::stringstream ss(tout_s); double val;
        while (ss >> val) tout_.push_back(val);
    }
    if (!tout_.empty()) tend_ = std::max(tend_, tout_.back());

    // Parse nsubcycle
    std::string nsub_s = config_.get("run_params", "nsubcycle", "");
    if (!nsub_s.empty()) {
        std::replace(nsub_s.begin(), nsub_s.end(), ',', ' ');
        std::replace(nsub_s.begin(), nsub_s.end(), '*', ' ');
        std::stringstream ss(nsub_s); int val;
        while (ss >> val) nsubcycle_.push_back(val);
    }
}

void Simulation::run() {
    std::cout << "[Simulation] Starting main loop..." << std::endl;
    real_t gamma = grid_.gamma;
    
    // Step 0 Diagnostics
    {
        real_t mind, maxv, mint, maxt;
        hydro_.get_diagnostics(p::levelmin, 0.0, mind, maxv, mint, maxt);
        printf(" Fine step=      0 t=%12.5e dt=%10.3e level=%2d min_d=%9.2e max_v=%9.2e min_T=%9.2e max_T=%9.2e\n", 
               t_, 0.0, p::levelmin, mind, maxv, mint, maxt);
    }

    dump_snapshot(1);
    int iout = 0, snapshot_count = 2;

    while (t_ < tend_ && nstep_ < nstepmax_) {
        nstep_++;
        real_t dt = 1e30;

        for (int il = p::levelmin; il <= grid_.nlevelmax; ++il) {
            if (grid_.count_grids_at_level(il) == 0) continue;
            real_t dx_cell = p::boxlen / static_cast<real_t>(p::nx * (1 << (il - 1)));
#ifdef MHD
            real_t ldt = mhd_.compute_courant_step(il, dx_cell, gamma, config_.get_double("hydro_params", "courant_factor", 0.8));
#else
            real_t ldt = hydro_.compute_courant_step(il, dx_cell, gamma, config_.get_double("hydro_params", "courant_factor", 0.8));
#endif
            dt = std::min(dt, ldt);
        }
        
        if (iout < (int)tout_.size()) {
            if (t_ + dt >= tout_[iout]) {
                dt = tout_[iout] - t_;
                amr_step(p::levelmin, dt);
                t_ += dt;
                
                // Restriction for levels below levelmin
                for (int il = p::levelmin - 1; il >= 1; --il) {
                    if (il == 1) grid_.restrict_coarse();
                    else grid_.restrict_fine(il);
                }

                dump_snapshot(snapshot_count++);
                iout++;
                continue;
            }
        }
        if (t_ + dt > tend_) dt = tend_ - t_;
        if (dt < 1e-15 * tend_ && nstep_ > 0) {
            std::cout << "[Simulation] dt too small, finishing." << std::endl;
            break;
        }
        amr_step(p::levelmin, dt);
        t_ += dt;

        // Restriction for levels below levelmin
        for (int il = p::levelmin - 1; il >= 1; --il) {
            if (il == 1) grid_.restrict_coarse();
            else grid_.restrict_fine(il);
        }

        // Step summary diagnostics
        if (nstep_ % ncontrol_ == 0) {
            real_t mind, maxv, mint, maxt;
            hydro_.get_diagnostics(p::levelmin, 0.0, mind, maxv, mint, maxt);
            printf(" Fine step=%7d t=%12.5e dt=%10.3e level=%2d min_d=%9.2e max_v=%9.2e min_T=%9.2e max_T=%9.2e\n", 
                   nstep_, t_, dt, p::levelmin, mind, maxv, mint, maxt);
        }

        // Dynamic Load Balancing
        if (MpiManager::instance().size() > 1) {
            balancer_.calculate_hilbert_keys();
            balancer_.balance();
        }
    }
    if (iout < (int)tout_.size() || t_ >= tend_) dump_snapshot(snapshot_count);
}

void Simulation::amr_step(int ilevel, real_t dt, int icount) {
    if (ilevel > grid_.nlevelmax) return;
    if (grid_.count_grids_at_level(ilevel) == 0) return;

    real_t dx = p::boxlen / static_cast<real_t>(p::nx * (1 << (ilevel - 1)));

    // Refinement generation (beginning of step)
    if (ilevel == p::levelmin || icount > 1) {
        for (int i = ilevel; i < grid_.nlevelmax; ++i) {
            updater_.refine_fine(i);
        }
    }

    if (grid_.count_grids_at_level(ilevel) > 0) {
#ifdef MHD
        mhd_.set_unew(ilevel);
#else
        hydro_.set_unew(ilevel);
#endif
        rt_.set_unew(ilevel);
    }

    int nsub = 2;
    if (ilevel <= (int)nsubcycle_.size()) nsub = nsubcycle_[ilevel - 1];

    if (ilevel < grid_.nlevelmax && grid_.count_grids_at_level(ilevel + 1) > 0) {
        for (int isub = 1; isub <= nsub; ++isub) {
            amr_step(ilevel + 1, dt / static_cast<real_t>(nsub), isub);
        }
    }

    if (grid_.count_grids_at_level(ilevel) > 0) {
#ifdef MHD
        mhd_.godunov_fine(ilevel, dt, dx);
        mhd_.set_uold(ilevel);
#else
        hydro_.godunov_fine(ilevel, dt, dx);
        hydro_.set_uold(ilevel);
#endif
        rt_.set_uold(ilevel);
    }

    // Restriction
    if (ilevel == 1) grid_.restrict_coarse();
    else grid_.restrict_fine(ilevel); // restrict current level to father

    // Flagging for next step
    if (grid_.count_grids_at_level(ilevel) > 0 && ilevel < grid_.nlevelmax) {
        updater_.mark_cells(ilevel);
    }
}

void Simulation::dump_snapshot(int iout) {
    std::string dir = "output_" + std::string(5 - std::to_string(iout).length(), '0') + std::to_string(iout);
    mkdir(dir.c_str(), 0777);
    
    SnapshotInfo info;
    info.t = t_;
    info.nstep = nstep_;
    info.nstep_coarse = nstep_; // For now, assume coarse step = step
    info.noutput = tout_.size();
    info.iout = iout;
    info.gamma = grid_.gamma;
    info.tout = tout_;

    int myid = MpiManager::instance().rank() + 1;
    std::string snap_str = std::string(5 - std::to_string(iout).length(), '0') + std::to_string(iout);

    RamsesWriter amr_writer(dir + "/amr_" + snap_str + ".out00001");
    amr_writer.write_amr(grid_, info);

    RamsesWriter hydro_writer(dir + "/hydro_" + snap_str + ".out00001");
    hydro_writer.write_hydro(grid_, info);

    RamsesWriter info_writer(dir + "/info_" + snap_str + ".txt");
    info_writer.write_header(info);

    RamsesWriter desc_writer(dir + "/hydro_file_descriptor.txt");
    desc_writer.write_hydro_descriptor(grid_, info);

    RamsesWriter head_writer(dir + "/header_" + snap_str + ".txt");
    head_writer.write_extra_headers(info);

    std::cout << "[Simulation] Finished dump_snapshot " << iout << std::endl;
}

real_t Simulation::compute_total_mass() { return 0.0; }
real_t Simulation::compute_total_energy() { return 0.0; }
real_t Simulation::compute_potential_energy() { return 0.0; }

} // namespace ramses
