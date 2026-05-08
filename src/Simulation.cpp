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

namespace ramses {

namespace p = ramses::params;

Simulation::Simulation() : grid_(), 
                           updater_(grid_, config_), 
                           load_balancer_(grid_, config_),
                           cosmo_() {
    hydro_ = create_hydro_solver(grid_, config_);
    cooling_ = create_cooling_solver(grid_, config_);
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
    int npassive = config_.get_int("hydro_params", "npassive", 0);
#ifdef RAMSES_NMETALS
    npassive = RAMSES_NMETALS;
#endif
    npassive = config_.get_int("hydro_params", "nmetals", npassive);

    int nvar = 5 + nener_ + npassive;
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
    
    bool do_poisson = config_.get_bool("run_params", "poisson", false);

    int ncpu = MpiManager::instance().size();
    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar, ncpu, levelmax);
    if (nparttot > 0) grid_.resize_particles(nparttot);
    hydro_->set_nener(nener_);
    hydro_->set_nvar_hydro(5 + nener_ + npassive);
    
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
        initializer_->apply_all();
        for (int il = 1; il < levelmax; ++il) {
            updater_.flag_fine(il, ed, ep, ev, eb2, {}, nexpand_[il]);
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
    for (int il = levelmax - 1; il >= 1; --il) updater_.restrict_fine(il);
    particles_->relink();

    tend_ = config_.get_double("run_params", "tend", 1e10);
    nstepmax_ = config_.get_int("run_params", "ncontrol", 1000000);
    ncontrol_ = 1;
    
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
    
    std::string tout_s = config_.get("output_params", "tout", "");
    if (!tout_s.empty()) {
        tout_.clear(); std::replace(tout_s.begin(), tout_s.end(), ',', ' ');
        std::stringstream ss(tout_s); double t; while(ss >> t) tout_.push_back(t);
    }
    if (tout_.empty()) tout_.push_back(tend_);
}

void Simulation::run() {
    real_t courant = config_.get_double("hydro_params", "courant_factor", 0.8);
    int snapshot_count = 1;
    int iout = 0;
    
    dump_snapshot(snapshot_count++);

    while (t_ < tout_.back() && nstep_ < nstepmax_) {
        bool do_poisson = config_.get_bool("run_params", "poisson", false);
        if (do_poisson) {
            for (int il = grid_.nlevelmax; il >= 1; --il) rho_fine(il);
            
            for (int il = grid_.nlevelmax; il >= 2; --il) {
                int myid = MpiManager::instance().rank() + 1, n2d_val = (1 << NDIM);
                int ig = grid_.get_headl(myid, il);
                while (ig > 0) {
                    int id_p = grid_.father[ig - 1];
                    if (id_p > 0) {
                        real_t sum = 0;
                        for (int ic = 1; ic <= n2d_val; ++ic) {
                            size_t idx = (size_t)grid_.ncoarse + (size_t)(ic - 1) * grid_.ngridmax + ig - 1;
                            sum += grid_.rho[idx];
                        }
                        grid_.rho[id_p - 1] += sum / (real_t)n2d_val;
                    }
                    ig = grid_.next[ig - 1];
                }
            }
            
            real_t omega_m = config_.get_double("cosmo_params", "omega_m", 0.3);
            real_t rho_tot = config_.get_double("poisson_params", "rho_tot", 0.0);
            if (config_.get_bool("run_params", "cosmo", false)) rho_tot = 1.0;
            
            for (int il = 1; il <= grid_.nlevelmax; ++il) {
                poisson_->solve(il, aexp_, omega_m, rho_tot);
                poisson_->compute_force(il);
            }
        }

        real_t min_dt = 1e10;
        for (int il = 1; il <= grid_.nlevelmax; ++il) {
            if (grid_.count_grids_at_level(il) == 0 && il > 1) continue;
            real_t dx = p::boxlen / (real_t)(p::nx * (1 << (il - 1)));
            
            if (il == 1) {
                for (int i = 1; i <= grid_.ncoarse; ++i) {
                    if (grid_.son[i - 1] > 0) continue;
                    real_t d = std::max(grid_.uold(i, 1), 1e-10), v2 = 0;
                    for (int j = 1; j <= NDIM; ++j) { real_t v = grid_.uold(i, 1 + j) / d; v2 += v * v; }
                    real_t p = std::max((grid_.uold(i, 5) - 0.5 * d * v2) * (grid_.gamma - 1.0), d * 1e-10);
                    real_t c = std::sqrt(grid_.gamma * p / d);
                    min_dt = std::min(min_dt, courant * dx / (std::sqrt(v2) + c));
                    for (int j = 1; j <= NDIM; ++j) {
                        real_t acc = std::abs(grid_.f(i, j));
                        if (acc > 0) min_dt = std::min(min_dt, courant * std::sqrt(dx / acc));
                    }
                }
            } else {
                int myid = MpiManager::instance().rank() + 1;
                int ig = grid_.get_headl(myid, il);
                int n2d_val = (1 << NDIM);
                while (ig > 0) {
                    for (int ic = 1; ic <= n2d_val; ++ic) {
                        int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
                        if (grid_.son[idc - 1] > 0) continue;

                        real_t d = std::max(grid_.uold(idc, 1), 1e-10), v2 = 0;
                        for(int j=1; j<=NDIM; ++j) { real_t v = grid_.uold(idc, 1+j)/d; v2 += v*v; }
                        real_t p = std::max((grid_.uold(idc, 5) - 0.5*d*v2)*(grid_.gamma-1.0), d * 1e-10);
                        real_t c = std::sqrt(grid_.gamma*p/d);
                        min_dt = std::min(min_dt, courant * dx / (std::sqrt(v2) + c));
                        for (int j = 1; j <= NDIM; ++j) {
                            real_t acc = std::abs(grid_.f(idc, j));
                            if (acc > 0) min_dt = std::min(min_dt, courant * std::sqrt(dx / acc));
                        }
                    }
                    ig = grid_.next[ig - 1];
                }
            }
        }

        real_t global_min_dt = min_dt;
#ifdef RAMSES_USE_MPI
        if (MpiManager::instance().size() > 1) {
            MPI_Allreduce(&min_dt, &global_min_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        }
#endif
        real_t dt = global_min_dt;
        if (iout < (int)tout_.size()) dt = std::min(dt, tout_[iout] - t_);

        amr_step(1, dt); t_ += dt; nstep_++;

        if (config_.get_bool("run_params", "cosmo", false)) {
            real_t texp;
            cosmo_.get_cosmo_params(t_, aexp_, hexp_, texp);
        }

        particles_->relink();

        if (iout < (int)tout_.size() && t_ >= tout_[iout] - 1e-12 * t_) {
            dump_snapshot(snapshot_count++); iout++;
        }
        if (nstep_ % ncontrol_ == 0) printf(" Step=%d t=%12.5e dt=%10.3e\n", nstep_, t_, dt);
    }
    dump_snapshot(snapshot_count);
}

void Simulation::amr_step(int ilevel, real_t dt, int icount) {
    if (ilevel > grid_.nlevelmax) return;
    bool grids_exist = (ilevel == 1) || (grid_.count_grids_at_level(ilevel) > 0);
    if (!grids_exist) return;

    real_t dx = p::boxlen / (real_t)(p::nx * (1 << (ilevel - 1)));
    bool do_poisson = config_.get_bool("run_params", "poisson", false);

    if (do_poisson) {
        hydro_->synchro_hydro_fine(ilevel, -0.5 * dt);
        particles_->move_fine(ilevel, 0.5 * dt);
    }

#ifdef MHD
    mhd_->godunov_fine(ilevel, dt, dx);
#else
    hydro_->godunov_fine(ilevel, dt, dx);
#endif

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
    for (int i = 1; i <= nsub; ++i) amr_step(ilevel + 1, dt / nsub, i);

    updater_.restrict_fine(ilevel);
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
    if (ilevel == 1) {
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
        if (ilevel == 1) {
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
    if (ilevel == 1 && MpiManager::instance().size() > 1) {
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
    info.noutput = 100; 
    info.iout = iout; 
    info.gamma = grid_.gamma; 
    info.nener = nener_;
    info.aexp = aexp_;
    info.hexp = hexp_;
    
    if (config_.get_bool("run_params", "cosmo", false)) {
        info.omega_m = config_.get_double("cosmo_params", "omega_m", 0.3);
        info.omega_l = config_.get_double("cosmo_params", "omega_l", 0.7);
        info.omega_k = config_.get_double("cosmo_params", "omega_k", 0.0);
        info.h0 = config_.get_double("cosmo_params", "h0", 70.0);
    }
    auto get_path = [&](const std::string& prefix, const std::string& ext) -> std::string {
        std::stringstream ss; ss << dir << "/" << prefix << "_" << std::setfill('0') << std::setw(5) << iout << ext << std::setfill('0') << std::setw(5) << 1;
        return ss.str();
    };
    RamsesWriter(get_path("amr", ".out")).write_amr(grid_, info);
    RamsesWriter(get_path("hydro", ".out")).write_hydro(grid_, info);
    if (grid_.npart > 0) {
        RamsesWriter(get_path("part", ".out")).write_particles(grid_, info);
        RamsesWriter(dir + "/part_file_descriptor.txt").write_particles_descriptor(grid_, info);
    }
#ifdef RT
    int nGroups = rt_->get_nGroups();
    if (nGroups > 0) {
        RamsesWriter(get_path("rt", ".out")).write_rt(grid_, info, nGroups, rt_->get_c_speed());
        RamsesWriter(dir + "/rt_file_descriptor.txt").write_rt_descriptor(grid_, info, nGroups);
    }
#endif
    std::stringstream ssinfo; ssinfo << dir << "/info_" << std::setfill('0') << std::setw(5) << iout << ".txt";
    RamsesWriter(ssinfo.str()).write_header(grid_, info);
    RamsesWriter(dir + "/hydro_file_descriptor.txt").write_hydro_descriptor(grid_, info);
    std::stringstream ssh; ssh << dir << "/header_" << std::setfill('0') << std::setw(5) << iout << ".txt";
    RamsesWriter(ssh.str()).write_header_file(grid_, info);
}

} // namespace ramses
