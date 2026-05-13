#include "ramses/Simulation.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/RamsesWriter.hpp"
#include "ramses/SolverFactory.hpp"
#ifdef RAMSES_USE_MPI
#include <mpi.h>
#endif
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <sys/stat.h>
#include <cmath>
#include <memory>
#include <chrono>

namespace ramses {

namespace p = ramses::params;

Simulation::Simulation() : grid_(), 
                           updater_(grid_, config_), 
                           load_balancer_(grid_, config_),
                           cosmo_() {
    hydro_ = create_hydro_solver(grid_, config_);
    cooling_ = create_cooling_solver(grid_, config_);
    turb_ = create_turbulence_solver(grid_, config_);
    sink_ = create_sink_solver(grid_, config_);
    star_ = create_star_solver(grid_, config_);
    feedback_ = create_feedback_solver(grid_, config_);
    clump_finder_ = create_clump_finder(grid_, config_);
    light_cone_ = create_light_cone(grid_, config_);
    mhd_ = create_mhd_solver(grid_, config_);
    rt_ = create_rt_solver(grid_, config_);
    poisson_ = create_poisson_solver(grid_, config_);
    particles_ = create_particle_solver(grid_, config_);
    initializer_ = create_initializer(grid_, config_);
}

void Simulation::initialize(const std::string& nml_path) {
    if (!config_.parse(nml_path)) return;
    
    p::nx = config_.get_int("amr_params", "nx", 1);
    p::ny = config_.get_int("amr_params", "ny", 1);
    p::nz = config_.get_int("amr_params", "nz", 1);
    p::boxlen = config_.get_double("amr_params", "boxlen", 1.0);

    int ngridmax = config_.get_int("amr_params", "ngridmax", 1000);
    ngridmax = config_.get_int("amr_params", "ngridtot", ngridmax);
    
    nener_ = config_.get_int("hydro_params", "nener", 0);
#ifdef RAMSES_NENER
    if (RAMSES_NENER > 0) nener_ = RAMSES_NENER;
#endif

    int nmetals = config_.get_int("hydro_params", "nmetals", 0);
    int npassive = config_.get_int("hydro_params", "npassive", nmetals);

#ifdef RAMSES_NPSCAL
    if (RAMSES_NPSCAL > 0 && npassive == 0) npassive = (int)RAMSES_NPSCAL;
#endif
#ifdef RAMSES_NMETALS
    if (RAMSES_NMETALS > 0 && npassive == 0) npassive = (int)RAMSES_NMETALS;
#endif

    int nvar = NDIM + 2 + nener_ + npassive;
#ifdef MHD
    nvar = 11 + nener_ + npassive;
#endif
#ifdef RT
    rt_->initialize();
    int nGroups = rt_->get_nGroups();
    int nIons = rt_->get_nIons();
    nvar += nIons + nGroups * (1 + NDIM);
#endif
    
    int levelmin = config_.get_int("amr_params", "levelmin", 1);
    int levelmax = config_.get_int("amr_params", "levelmax", 1);
    params::levelmin = levelmin; params::nlevelmax = levelmax;
    
    int nparttot = config_.get_int("amr_params", "nparttot", 0);
    
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    grid_.gamma = gamma;
    
    // Barotropic EOS parameters
    p::barotropic_eos = config_.get_bool("cooling_params", "barotropic_eos", false);
    p::barotropic_eos_form = config_.get("cooling_params", "barotropic_eos_form", "isothermal");
    p::polytrope_rho = config_.get_double("cooling_params", "polytrope_rho", 1e-15);
    p::polytrope_index = config_.get_double("cooling_params", "polytrope_index", 1.4);
    p::T_eos = config_.get_double("cooling_params", "T_eos", 10.0);
    p::mu_gas = config_.get_double("cooling_params", "mu_gas", 1.0);

    // Tracer parameters
    p::tracer = config_.get_bool("run_params", "tracer", false);
    p::MC_tracer = config_.get_bool("tracer_params", "mc_tracer", false);
    p::tracer_feed = config_.get_int("tracer_params", "tracer_feed", 0);
    p::tracer_feed_fmt = config_.get("tracer_params", "tracer_feed_fmt", "inplace");
    p::tracer_mass = config_.get_double("tracer_params", "tracer_mass", 0.0);
    p::tracer_first_balance_part_per_cell = config_.get_int("tracer_params", "tracer_first_balance_part_per_cell", 0);
    p::tracer_first_balance_levelmin = config_.get_int("tracer_params", "tracer_first_balance_levelmin", 0);

    // Physical units
    p::units_density = config_.get_double("units_params", "units_density", 1.0);
    p::units_time = config_.get_double("units_params", "units_time", 1.0);
    p::units_length = config_.get_double("units_params", "units_length", 1.0);
    p::units_velocity = p::units_length / p::units_time;
    p::units_mass = p::units_density * std::pow(p::units_length, 3);
    p::units_energy = p::units_mass * std::pow(p::units_velocity, 2);
    p::units_pressure = p::units_density * std::pow(p::units_velocity, 2);
    p::units_number_density = p::units_density / (p::mu_gas * constants::mH);

    p::scale_l = p::units_length;
    p::scale_t = p::units_time;
    p::scale_d = p::units_density;
    p::scale_v = p::units_velocity;
    p::scale_nH = p::units_number_density;
    p::scale_T2 = p::units_pressure / p::units_density; // Simplified scale for T/mu

    int ncpu = MpiManager::instance().size();
    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar, ncpu, levelmax);
    if (nparttot > 0) grid_.resize_particles(nparttot);
    hydro_->set_nener(nener_); hydro_->set_nvar_hydro(nvar);
    grid_.set_interpol_hook([this](const real_t u1[7][64], real_t u2[8][64]){ this->hydro_->interpol_hydro(u1, u2); });

    initializer_->apply_all();
    initializer_->init_tracers();

    nsubcycle_.assign(33, 1);
    std::string nsub_s = config_.get("run_params", "nsubcycle", "");
    if (!nsub_s.empty()) {
        std::stringstream ss(nsub_s); std::string item; int l = 1;
        while (std::getline(ss, item, ',')) {
            size_t star = item.find('*');
            if (star != std::string::npos) {
                int count = std::stoi(item.substr(0, star));
                int val = std::stoi(item.substr(star + 1));
                for(int i=0; i<count; ++i) if(l+i < 33) nsubcycle_[l+i] = val;
                l += count;
            } else {
                if(l < 33) nsubcycle_[l++] = std::stoi(item);
            }
        }
    }

    nexpand_.assign(33, 1);
    std::string nexp_s = config_.get("amr_params", "nexpand", "");
    if (!nexp_s.empty()) {
        std::stringstream ss(nexp_s); std::string item; int l = 1;
        while (std::getline(ss, item, ',')) {
            size_t star = item.find('*');
            if (star != std::string::npos) {
                int count = std::stoi(item.substr(0, star));
                int val = std::stoi(item.substr(star + 1));
                for(int i=0; i<count; ++i) if(l+i < 33) nexpand_[l+i] = val;
                l += count;
            } else {
                if(l < 33) nexpand_[l++] = std::stoi(item);
            }
        }
    }

    real_t ed = config_.get_double("refine_params", "err_grad_d", -1.0);
    real_t ep = config_.get_double("refine_params", "err_grad_p", -1.0);
    real_t ev = config_.get_double("refine_params", "err_grad_v", -1.0);
    real_t eb2 = config_.get_double("refine_params", "err_grad_b2", -1.0);

    for (int ipass = 1; ipass <= 2; ++ipass) {
        for (int il = 1; il <= levelmax; ++il) {
            initializer_->apply_all(); // Ensure current level is initialized
            updater_.flag_fine(il, ed, ep, ev, eb2, {}, nexpand_[std::min(il, 32)]);
            updater_.make_grid_fine(il);
            updater_.remove_grid_fine(il);
        }
        grid_.synchronize_level_counts();
        if (MpiManager::instance().size() > 1) {
            load_balancer_.calculate_hilbert_keys();
            load_balancer_.balance();
        }
        particles_->relink();
    }

    nstep_ = 0;
    initializer_->apply_all();
    if (config_.get_bool("run_params", "turb", false)) turb_->init();
    if (config_.get_bool("run_params", "sink", false)) sink_->init();
    for (int il = levelmax - 1; il >= 1; --il) updater_.restrict_fine(il);
    particles_->relink();

    tend_ = config_.get_double("run_params", "tend", 1e10);
    p::tend = tend_;
    nstepmax_ = config_.get_int("run_params", "nstepmax", 1000000);
    ncontrol_ = config_.get_int("run_params", "ncontrol", 1);
    
    if (MpiManager::instance().rank() == 0 && config_.get_bool("run_params", "verbose", false)) {
        std::cout << "[Simulation] nstepmax=" << nstepmax_ << " tend=" << tend_ << " ncontrol=" << ncontrol_ << std::endl;
    }

    std::string tout_s = config_.get("output_params", "tout", "");
    if (!tout_s.empty()) {
        tout_.clear(); std::replace(tout_s.begin(), tout_s.end(), ',', ' ');
        std::stringstream ss(tout_s); double t; while(ss >> t) tout_.push_back(t);
    }
    if (tout_.empty()) tout_.push_back(tend_);

    if (tend_ > tout_.back()) {
        tend_ = tout_.back();
    }
    
    bool do_cosmo = config_.get_bool("run_params", "cosmo", false);
    if (do_cosmo) {
        real_t omega_m = config_.get_double("cosmo_params", "omega_m", 0.3);
        real_t omega_l = config_.get_double("cosmo_params", "omega_l", 0.7);
        real_t omega_k = config_.get_double("cosmo_params", "omega_k", 0.0);
        real_t aexp_ini = config_.get_double("cosmo_params", "aexp_ini", 1e-3);
        cosmo_.solve_friedman(omega_m, omega_l, omega_k, aexp_ini);
        
        t_ = cosmo_.get_tau(aexp_ini);
        aexp_ = aexp_ini;
        real_t texp;
        cosmo_.get_cosmo_params(t_, aexp_, hexp_, texp);
    }
}

void Simulation::run() {
    real_t courant = config_.get_double("hydro_params", "courant_factor", 0.8);
    int snapshot_count = 1;
    int iout = 0;
    
    dump_snapshot(snapshot_count++);
    bool verbose = config_.get_bool("run_params", "verbose", false);

    if (MpiManager::instance().rank() == 0) {
        if (verbose) {
            std::cout << "[Simulation] Starting simulation with NDIM=" << NDIM << std::endl;
            std::cout << "[Simulation] Boxlen=" << p::boxlen << " Levelmin=" << p::levelmin << " Levelmax=" << p::nlevelmax << std::endl;
            std::cout << "[Simulation] Gamma=" << grid_.gamma << " Courant=" << courant << std::endl;
        }
        std::cout << "Initial mesh structure" << std::endl;
        for (int il = 1; il <= grid_.nlevelmax; ++il) {
            int ngrids = grid_.count_grids_at_level(il);
            if (ngrids > 0) {
                std::cout << " Level " << std::setw(2) << il << " has " << std::setw(10) << ngrids << " grids" << std::endl;
            }
        }
        std::cout << "Starting time integration" << std::endl;
    }

    auto t_start_loop = std::chrono::high_resolution_clock::now();
    while (t_ < tout_.back() && nstep_ < nstepmax_) {
        auto t_step_start = std::chrono::high_resolution_clock::now();
        real_t min_dt = 1e10;
        // Level 0 (ncoarse)
        real_t dx0 = p::boxlen / (real_t)p::nx;
        min_dt = std::min(min_dt, hydro_->compute_courant_step(0, dx0, grid_.gamma, courant));

        for (int il = 1; il <= grid_.nlevelmax; ++il) {
            if (grid_.count_grids_at_level(il) == 0) continue;
            real_t dx = p::boxlen / (real_t)(p::nx * (1 << il));
            min_dt = std::min(min_dt, hydro_->compute_courant_step(il, dx, grid_.gamma, courant));
        }

        real_t global_min_dt = min_dt;
#ifdef RAMSES_USE_MPI
        if (MpiManager::instance().size() > 1) {
            MPI_Allreduce(&min_dt, &global_min_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        }
#endif
        real_t dt = global_min_dt;
        if (iout < (int)tout_.size()) dt = std::min(dt, tout_[iout] - t_);

        amr_step(0, dt); t_ += dt; p::t = t_; nstep_++;

        if (config_.get_bool("run_params", "cosmo", false)) {
            real_t texp;
            cosmo_.get_cosmo_params(t_, aexp_, hexp_, texp);
        }

        particles_->relink();

        if (iout < (int)tout_.size() && t_ >= tout_[iout] - 1e-10 * dt) {
            dump_snapshot(snapshot_count++); iout++;
        }
        if (nstep_ % ncontrol_ == 0) {
            auto t_step_end = std::chrono::high_resolution_clock::now();
            double duration = std::chrono::duration<double>(t_step_end - t_step_start).count();
            double total_time = std::chrono::duration<double>(t_step_end - t_start_loop).count();
            if (MpiManager::instance().rank() == 0) {
                printf(" Step=%d t=%12.5e dt=%10.3e time_elapsed=%8.2fs total_time=%8.2fs\n", nstep_, t_, dt, duration, total_time);
            }
        }
    }
    dump_snapshot(snapshot_count++);

    if (MpiManager::instance().rank() == 0) {
        auto t_end_loop = std::chrono::high_resolution_clock::now();
        double total_duration = std::chrono::duration<double>(t_end_loop - t_start_loop).count();
        std::cout << "Simulation finished in " << std::fixed << std::setprecision(2) << total_duration << " seconds." << std::endl;
        std::cout << "Total steps: " << nstep_ << " Final time: " << std::scientific << std::setprecision(6) << t_ << std::endl;
    }
}

void Simulation::amr_step(int ilevel, real_t dt, int icount) {
    if (ilevel > grid_.nlevelmax) return;
    bool grids_exist = (ilevel == 0) || (grid_.count_grids_at_level(ilevel) > 0);
    if (!grids_exist) return;

    if (config_.get_bool("run_params", "verbose", false) && MpiManager::instance().rank() == 0) {
        std::cout << " Entering amr_step(" << icount << ") for level " << ilevel << std::endl;
    }

    // Dynamic Refinement (Phase 25)
    if (ilevel < grid_.nlevelmax) {
        real_t ed = config_.get_double("refine_params", "err_grad_d", -1.0);
        real_t ep = config_.get_double("refine_params", "err_grad_p", -1.0);
        real_t ev = config_.get_double("refine_params", "err_grad_v", -1.0);
        real_t eb2 = config_.get_double("refine_params", "err_grad_b2", -1.0);
        int nexp = config_.get_int("amr_params", "nexpand", 1);
        
        updater_.flag_fine(ilevel + 1, ed, ep, ev, eb2, {}, nexp);
        updater_.make_grid_fine(ilevel + 1);
        updater_.remove_grid_fine(ilevel + 1);
        grid_.synchronize_level_counts();
    }

    real_t dx = p::boxlen / (real_t)(p::nx * (1 << ilevel));
    bool do_poisson = config_.get_bool("run_params", "poisson", false);

#ifdef MHD
    mhd_->set_unew(ilevel);
#else
    hydro_->set_unew(ilevel);
#endif

    if (do_poisson) {
        hydro_->synchro_hydro_fine(ilevel, -0.5 * dt);

        // Calculate total mass and mean density for Poisson stability
        real_t total_mass = 0, total_vol = 0;
        int n2d_val = (1 << NDIM);
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            real_t dx_c = p::boxlen / (real_t)p::nx;
            real_t dV = std::pow(dx_c, NDIM);
            if (grid_.son[i-1] == 0) { total_mass += grid_.uold(i, 1) * dV; total_vol += dV; }
        }
        for (int il = 1; il <= grid_.nlevelmax; ++il) {
            int my_id = MpiManager::instance().rank() + 1;
            int ig = grid_.get_headl(my_id, il);
            real_t dx_l = p::boxlen / (real_t)(p::nx * (1 << il));
            real_t dV = std::pow(dx_l, NDIM);
            while (ig > 0) {
                for (int ic = 1; ic <= n2d_val; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                    if (grid_.son[idc] == 0) { total_mass += grid_.uold(idc + 1, 1) * dV; total_vol += dV; }
                }
                ig = grid_.next[ig - 1];
            }
        }
        real_t global_mass = total_mass, global_vol = total_vol;
#ifdef RAMSES_USE_MPI
        if (MpiManager::instance().size() > 1) {
            MPI_Allreduce(&total_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&total_vol, &global_vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
#endif
        real_t rho_mean = global_mass / std::max(global_vol, 1e-10);
        
        // Update density and aggregate for multi-grid
        for (int il = ilevel; il >= 0; --il) rho_fine(il);
        for (int il = ilevel; il >= 1; --il) {
            int my_id = MpiManager::instance().rank() + 1;
            int ig = grid_.get_headl(my_id, il);
            while (ig > 0) {
                int id_p = grid_.father[ig - 1];
                if (id_p > 0) {
                    real_t sum = 0;
                    for (int ic = 1; ic <= n2d_val; ++ic) {
                        size_t idx = (size_t)grid_.ncoarse + (size_t)(ic - 1) * grid_.ngridmax + ig - 1;
                        sum += grid_.rho[idx];
                    }
                    grid_.rho[id_p - 1] += sum;
                }
                ig = grid_.next[ig - 1];
            }
        }
        
        real_t omega_m = config_.get_double("cosmo_params", "omega_m", 0.3);
        poisson_->solve(ilevel, aexp_, omega_m, rho_mean);
        poisson_->compute_force(ilevel);

        hydro_->synchro_hydro_fine(ilevel, 0.5 * dt);
        particles_->move_fine(ilevel, 0.5 * dt);
    }

#ifdef MHD
    mhd_->godunov_fine(ilevel, dt, dx);
#else
    hydro_->godunov_fine(ilevel, dt, dx);
#endif

    if (config_.get_bool("run_params", "turb", false)) turb_->apply_forcing(ilevel, dt);

    if (config_.get_bool("run_params", "sink", false)) {
        sink_->create_sinks(ilevel);
        sink_->grow_sinks(ilevel, dt);
        sink_->synchronize_sinks();
    }

    if (config_.get_bool("run_params", "star", false)) {
        star_->form_stars(ilevel, dt);
    }

    if (config_.get_bool("run_params", "star", false) && config_.get_double("feedback_params", "eta_sn", 0.0) > 0) {
        feedback_->thermal_feedback(ilevel, dt);
    }

    if (config_.get_bool("run_params", "cosmo", false) && config_.get_bool("run_params", "lightcone", false)) {
        light_cone_->update(t_, dt);
    }

    if (do_poisson) {
        hydro_->add_gravity_source_terms(ilevel, dt);
        particles_->move_fine(ilevel, 0.5 * dt);
        particles_->exchange_particles();
    }

    cooling_->apply_cooling(ilevel, dt);

#ifdef RT
    rt_->godunov_fine(ilevel, dt, dx);
    rt_->apply_source_terms(ilevel, dt);
#endif

    int nsub = (ilevel < (int)nsubcycle_.size()) ? nsubcycle_[ilevel] : 1;
    if (ilevel < grid_.nlevelmax) {
        for (int i = 1; i <= nsub; ++i) amr_step(ilevel + 1, dt / nsub, i);
    }

    // RESTRICT FIRST: Update parents from children
    if (ilevel < grid_.nlevelmax) {
        updater_.restrict_fine(ilevel);
    }

    // THEN SET UOLD: Finalize the step for this level's cells.
#ifdef MHD
    mhd_->set_uold(ilevel);
#else
    hydro_->set_uold(ilevel);
#endif

    if (do_poisson) hydro_->synchro_hydro_fine(ilevel, 0.5 * dt);

#ifdef RT
    rt_->set_uold(ilevel);
#endif
}

void Simulation::rho_fine(int ilevel) {
    int myid = MpiManager::instance().rank() + 1, n2d_val = (1 << NDIM);
    if (ilevel == 0) {
        for (int i = 1; i <= grid_.ncoarse; ++i) grid_.rho[i - 1] = 0.0;
    } else {
        int igrid = grid_.get_headl(myid, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= n2d_val; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                grid_.rho[idc - 1] = 0.0;
            }
            igrid = grid_.next[igrid - 1];
        }
    }

    particles_->assign_mass_fine(ilevel);

    if (config_.get_bool("run_params", "hydro", true)) {
        if (ilevel == 0) {
            for (int i = 1; i <= grid_.ncoarse; ++i) if (grid_.son[i - 1] == 0) grid_.rho[i - 1] += grid_.uold(i, 1);
        } else {
            int igrid = grid_.get_headl(myid, ilevel);
            while (igrid > 0) {
                for (int ic = 1; ic <= n2d_val; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    if (grid_.son[idc - 1] == 0) grid_.rho[idc - 1] += grid_.uold(idc, 1);
                }
                igrid = grid_.next[igrid - 1];
            }
        }
    }

#ifdef RAMSES_USE_MPI
    if (ilevel == 0 && MpiManager::instance().size() > 1) {
        std::vector<real_t> global_rho(grid_.ncoarse);
        MPI_Allreduce(grid_.rho.data(), global_rho.data(), grid_.ncoarse, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        std::copy(global_rho.begin(), global_rho.end(), grid_.rho.begin());
    }
#endif
}

void Simulation::dump_snapshot(int iout) {
    std::stringstream ssd; ssd << "output_" << std::setfill('0') << std::setw(5) << iout;
    std::string dir = ssd.str(); mkdir(dir.c_str(), 0777);
    
    SnapshotInfo info; 
    info.t = t_; 
    info.nstep = nstep_; 
    info.iout = iout; 
    info.gamma = grid_.gamma; 
    info.nener = nener_;
    info.noutput = (int)tout_.size();
    info.aexp = aexp_;
    info.hexp = hexp_;
    
    auto get_path = [&](const std::string& prefix, const std::string& ext, bool use_rank) -> std::string {
        std::stringstream ss; 
        ss << dir << "/" << prefix << "_" << std::setfill('0') << std::setw(5) << iout << ext;
        if (use_rank) ss << std::setfill('0') << std::setw(5) << MpiManager::instance().rank() + 1;
        return ss.str();
    };

    RamsesWriter(get_path("amr", ".out", true)).write_amr(grid_, info);
    RamsesWriter(get_path("hydro", ".out", true)).write_hydro(grid_, info);
    RamsesWriter(get_path("grav", ".out", true)).write_grav(grid_, info);
    
    // Header for particles (required by visu_ramses even if npart=0)
    RamsesWriter(get_path("header", ".txt", false)).write_header_file(grid_, info);
    
    // Descriptors (once per rank or just rank 0?)
    RamsesWriter(dir + "/hydro_file_descriptor.txt").write_hydro_descriptor(grid_, info);
    RamsesWriter(dir + "/part_file_descriptor.txt").write_particles_descriptor(grid_, info);

    if (MpiManager::instance().rank() == 0) {
        std::stringstream ss_info;
        ss_info << dir << "/info_" << std::setfill('0') << std::setw(5) << iout << ".txt";
        RamsesWriter(ss_info.str()).write_header(grid_, info);
    }
}

} // namespace ramses
