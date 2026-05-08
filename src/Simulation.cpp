#include "ramses/Simulation.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/RamsesWriter.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <sys/stat.h>

namespace ramses {

namespace p = ramses::params;

void Simulation::initialize(const std::string& nml_path) {
    if (!config_.parse(nml_path)) return;
    
    p::nx = config_.get_int("amr_params", "nx", 1);
    p::ny = config_.get_int("amr_params", "ny", 1);
    p::nz = config_.get_int("amr_params", "nz", 1);
    p::boxlen = config_.get_double("amr_params", "boxlen", 1.0);

    int ngridmax = config_.get_int("amr_params", "ngridmax", 1000);
    ngridmax = config_.get_int("amr_params", "ngridtot", ngridmax);
    
    nener_ = config_.get_int("hydro_params", "nener", 0);
    int nvar = 5 + nener_;
#ifdef MHD
    nvar = 11 + nener_;
#endif
#ifdef RT
    rt_.initialize();
    int nGroups = rt_.get_nGroups();
    int nIons = 3; 
    if (config_.get_bool("rt_params", "isH2", false)) nIons = 6;
    nvar += nIons + nGroups * (1 + NDIM);
#endif
    
    int levelmin = config_.get_int("amr_params", "levelmin", 1);
    int levelmax = config_.get_int("amr_params", "levelmax", 1);
    params::levelmin = levelmin; params::nlevelmax = levelmax;
    
    int nparttot = config_.get_int("amr_params", "nparttot", 0);
    
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    grid_.gamma = gamma;
    
    bool do_poisson = config_.get_bool("run_params", "poisson", false);
    std::cout << "[Simulation] Poisson enabled: " << (do_poisson ? "YES" : "NO") << std::endl;

    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar, 1, levelmax);
    if (nparttot > 0) grid_.resize_particles(nparttot);
    hydro_.set_nener(nener_);
    hydro_.set_nvar_hydro(5 + nener_);
    
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

    for (int ipass = 1; ipass <= levelmax; ++ipass) {
        initializer_.apply_all();
        for (int il = 1; il < levelmax; ++il) {
            updater_.flag_fine(il, ed, ep, ev, eb2, {}, nexpand_[il]);
            updater_.make_grid_fine(il);
            updater_.remove_grid_fine(il);
        }
        particles_.relink();
    }

    for(int il=1; il<=levelmax; ++il) std::cout << "[Simulation] Level " << il << " ngrid=" << grid_.count_grids_at_level(il) << std::endl;

    nstep_ = 0;
    initializer_.apply_all();
    for (int il = levelmax - 1; il >= 1; --il) updater_.restrict_fine(il);
    
    particles_.relink();

    tend_ = config_.get_double("run_params", "tend", 1e10);
    nstepmax_ = config_.get_int("run_params", "ncontrol", 1000000);
    ncontrol_ = 1;
    
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
    
    // Dump initial state
    dump_snapshot(snapshot_count++);

    while (t_ < tout_.back() && nstep_ < nstepmax_) {
        bool do_poisson = config_.get_bool("run_params", "poisson", false);
        if (do_poisson) {
            for (int il = grid_.nlevelmax; il >= 1; --il) rho_fine(il);
            // Restrict rho
            for (int il = grid_.nlevelmax; il >= 2; --il) {
                int myid = 1, n2d_val = (1 << NDIM);
                int ig = grid_.get_headl(myid, il - 1);
                while (ig > 0) {
                    int id_p = grid_.father[ig - 1];
                    if (id_p > 0) {
                        real_t sum = 0;
                        for (int ic = 1; ic <= n2d_val; ++ic) sum += grid_.rho[grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1];
                        grid_.rho[id_p - 1] += sum / std::pow(2.0, NDIM);
                    }
                    ig = grid_.next[ig - 1];
                }
            }
            for (int il = 1; il <= grid_.nlevelmax; ++il) {
                poisson_.solve(il);
                poisson_.compute_force(il);
            }
        }

        real_t min_dt = 1e10;
        for (int il = 1; il <= grid_.nlevelmax; ++il) {
            bool cells_exist = (il == 1) || (grid_.count_grids_at_level(il - 1) > 0);
            if (!cells_exist) continue;
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
                int ig = grid_.get_headl(1, il - 1);
                while(ig > 0) {
                    for(int ic=1; ic<=(1<<NDIM); ++ic) {
                        int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig - 1;
                        if (grid_.son[idc] > 0) continue;
                        real_t d = std::max(grid_.uold(idc + 1, 1), 1e-10), v2 = 0;
                        for(int j=1; j<=NDIM; ++j) { real_t v = grid_.uold(idc + 1, 1+j)/d; v2 += v*v; }
                        real_t p = std::max((grid_.uold(idc + 1, 5) - 0.5*d*v2)*(grid_.gamma-1.0), d * 1e-10);
                        real_t c = std::sqrt(grid_.gamma*p/d);
                        min_dt = std::min(min_dt, courant * dx / (std::sqrt(v2) + c));
                        for (int j = 1; j <= NDIM; ++j) {
                            real_t acc = std::abs(grid_.f(idc + 1, j));
                            if (acc > 0) min_dt = std::min(min_dt, courant * std::sqrt(dx / acc));
                        }
                    }
                    ig = grid_.next[ig - 1];
                }
            }
        }
        
        real_t dt = min_dt;
        if (iout < (int)tout_.size()) dt = std::min(dt, tout_[iout] - t_);
        
        amr_step(1, dt); t_ += dt; nstep_++;
        particles_.relink();
        
        if (iout < (int)tout_.size() && t_ >= tout_[iout] - 1e-12 * t_) {
            dump_snapshot(snapshot_count++); iout++;
        }
        if (nstep_ % ncontrol_ == 0) printf(" Step=%d t=%12.5e dt=%10.3e\n", nstep_, t_, dt);
    }
    dump_snapshot(snapshot_count);
}

void Simulation::amr_step(int ilevel, real_t dt, int icount) {
    if (ilevel > grid_.nlevelmax) return;
    bool cells_exist = (ilevel == 1) || (grid_.count_grids_at_level(ilevel - 1) > 0);
    if (!cells_exist) return;
    
    real_t dx = p::boxlen / (real_t)(p::nx * (1 << (ilevel - 1)));
    bool do_poisson = config_.get_bool("run_params", "poisson", false);
    
    if (do_poisson) {
        hydro_.synchro_hydro_fine(ilevel, -0.5 * dt);
        particles_.move_fine(ilevel, 0.5 * dt);
    }

#ifdef MHD
    mhd_.godunov_fine(ilevel, dt, dx);
#else
    hydro_.godunov_fine(ilevel, dt, dx);
#endif

    if (do_poisson) {
        hydro_.add_gravity_source_terms(ilevel, dt);
        particles_.move_fine(ilevel, 0.5 * dt);
        particles_.exchange_particles();
    }

#ifdef RT
    rt_.godunov_fine(ilevel, dt, dx);
    rt_.apply_source_terms(ilevel, dt);
#endif
    
    int nsub = (ilevel < (int)nsubcycle_.size()) ? nsubcycle_[ilevel] : 1;
    for (int i = 1; i <= nsub; ++i) amr_step(ilevel + 1, dt / nsub, i);
    
    updater_.restrict_fine(ilevel);
#ifdef MHD
    mhd_.set_uold(ilevel);
#else
    hydro_.set_uold(ilevel);
#endif

    if (do_poisson) hydro_.synchro_hydro_fine(ilevel, 0.5 * dt);

#ifdef RT
    rt_.set_uold(ilevel);
#endif
}

void Simulation::rho_fine(int ilevel) {
    int myid = 1, n2d_val = (1 << NDIM);
    if (ilevel == 1) {
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            if (grid_.son[i - 1] == 0) grid_.rho[i - 1] = grid_.uold(i, 1);
            else grid_.rho[i - 1] = 0.0;
        }
    } else {
        int igrid = grid_.get_headl(myid, ilevel - 1);
        while (igrid > 0) {
            for (int ic = 1; ic <= n2d_val; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                if (grid_.son[idc - 1] == 0) grid_.rho[idc - 1] = grid_.uold(idc, 1);
                else grid_.rho[idc - 1] = 0.0;
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    particles_.assign_mass_fine(ilevel);
}

void Simulation::dump_snapshot(int iout) {
    std::stringstream ssd; ssd << "output_" << std::setfill('0') << std::setw(5) << iout;
    std::string dir = ssd.str(); mkdir(dir.c_str(), 0777);
    SnapshotInfo info; info.t = t_; info.nstep = nstep_; info.noutput = 100; info.iout = iout; info.gamma = grid_.gamma; info.nener = nener_;
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
    int nGroups = rt_.get_nGroups();
    if (nGroups > 0) {
        RamsesWriter(get_path("rt", ".out")).write_rt(grid_, info, nGroups, rt_.get_c_speed());
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
