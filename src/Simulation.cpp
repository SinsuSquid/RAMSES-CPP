#include "ramses/Simulation.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/RamsesWriter.hpp"
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <algorithm>
#include <sstream>
#include <iomanip>

namespace ramses {

namespace p = ramses::params;

void Simulation::initialize(const std::string& nml_path) {
    if (!config_.parse(nml_path)) return;
    
    p::nx = config_.get_int("amr_params", "nx", 1);
    p::ny = config_.get_int("amr_params", "ny", 1);
    p::nz = config_.get_int("amr_params", "nz", 1);
    p::boxlen = config_.get_double("amr_params", "boxlen", 1.0);
    int ngridmax = config_.get_int("amr_params", "ngridmax", 1000);
    nener_ = config_.get_int("hydro_params", "nener", 0);
    hydro_.set_nener(nener_); mhd_.set_nener(nener_);

    int levelmin = config_.get_int("amr_params", "levelmin", 1);
    int levelmax = config_.get_int("amr_params", "levelmax", 1);
    params::levelmin = levelmin; params::nlevelmax = levelmax;
    int ncpu = 1; 
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    grid_.gamma = gamma;
    
    // Calculate nvar: 5 (hydro) + nener (non-thermal energy) + npassive
    int npassive = 0;
    std::string var_s = config_.get("init_params", "var_region", "");
    if (!var_s.empty()) {
        std::stringstream ss(var_s); std::string item;
        while (std::getline(ss, item, ',')) npassive++;
    }
    int nreg = config_.get_int("init_params", "nregion", 1);
    if (nreg > 0) npassive /= nreg;

    int nvar = 5 + nener_ + npassive;
#ifdef MHD
    nvar = 11 + nener_ + npassive;
#endif
    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar, ncpu, levelmax);
    grid_.nvar = nvar;
    updater_.set_interpol_hook([this](const real_t u1[7][20], real_t u2[8][20]){
#ifdef MHD
        mhd_.interpol_mhd(u1, u2);
#else
        hydro_.interpol_hydro(u1, u2);
#endif
    });
    
    real_t ed = config_.get_double("refine_params", "err_grad_d", 0.05);
    real_t ep = config_.get_double("refine_params", "err_grad_p", 0.0);
    real_t ev = config_.get_double("refine_params", "err_grad_u", 0.0);
    real_t eb2 = config_.get_double("refine_params", "err_grad_b2", ep);
    std::vector<real_t> evar;
    std::string evar_s = config_.get("refine_params", "err_grad_var", "");
    if (!evar_s.empty()) {
        std::stringstream ss(evar_s); std::string item;
        while (std::getline(ss, item, ',')) evar.push_back(std::stod(item));
    }
    if (evar.empty() && npassive > 0) evar.assign(npassive, 0.05);

    std::vector<std::string> nexp_s = config_.get_string_array("amr_params", "nexpand");
    nexpand_.assign(grid_.nlevelmax, 1);
    int curr_l = 0;
    for (const auto& s : nexp_s) {
        if (s.find('*') != std::string::npos) {
            int count = std::stoi(s.substr(0, s.find('*'))), val = std::stoi(s.substr(s.find('*') + 1));
            for (int i = 0; i < count && curr_l < grid_.nlevelmax; ++i) nexpand_[curr_l++] = val;
        } else if (!s.empty()) { if (curr_l < grid_.nlevelmax) nexpand_[curr_l++] = std::stoi(s); }
    }

    // Multi-pass initial refinement loop
    for (int ipass = 1; ipass <= levelmax; ++ipass) {
        initializer_.apply_all();
        for (int il = 1; il < levelmax; ++il) {
            updater_.flag_fine(il, ed, ep, ev, eb2, evar);
            updater_.make_grid_fine(il);
        }
    }
    nstep_ = 0; // Reset step counter
    initializer_.apply_all();
    for (int il = levelmax - 1; il >= 1; --il) updater_.restrict_fine(il);

    grid_.nboundary = config_.get_int("boundary_params", "nboundary", 0);
    if (grid_.nboundary > 0) {
        auto parse_vec = [](const std::string& s) {
            std::vector<int> v; std::string sc = s; std::replace(sc.begin(), sc.end(), ',', ' ');
            std::stringstream ss(sc); int val; while (ss >> val) v.push_back(val); return v;
        };
        grid_.ibound_min = parse_vec(config_.get("boundary_params", "ibound_min", ""));
        grid_.ibound_max = parse_vec(config_.get("boundary_params", "ibound_max", ""));
        grid_.jbound_min = parse_vec(config_.get("boundary_params", "jbound_min", ""));
        grid_.jbound_max = parse_vec(config_.get("boundary_params", "jbound_max", ""));
        grid_.kbound_min = parse_vec(config_.get("boundary_params", "kbound_min", ""));
        grid_.kbound_max = parse_vec(config_.get("boundary_params", "kbound_max", ""));
        grid_.bound_type = parse_vec(config_.get("boundary_params", "bound_type", ""));
    }

    t_ = 0.0; tend_ = config_.get_double("output_params", "tend", 1.0);
    nstepmax_ = config_.get_int("run_params", "nstepmax", 1000000);
    ncontrol_ = config_.get_int("run_params", "ncontrol", 1);
    std::string tout_s = config_.get("output_params", "tout", "");
    if (!tout_s.empty()) {
        std::replace(tout_s.begin(), tout_s.end(), ',', ' ');
        std::stringstream ss(tout_s); double val; while (ss >> val) tout_.push_back(val);
        if (!tout_.empty()) tend_ = tout_.back();
    }

    std::vector<std::string> nsub_s = config_.get_string_array("run_params", "nsubcycle");
    nsubcycle_.assign(grid_.nlevelmax, 1);
    curr_l = 0;
    for (const auto& s : nsub_s) {
        if (s.find('*') != std::string::npos) {
            int count = std::stoi(s.substr(0, s.find('*'))), val = std::stoi(s.substr(s.find('*') + 1));
            for (int i = 0; i < count && curr_l < grid_.nlevelmax; ++i) nsubcycle_[curr_l++] = val;
        } else if (!s.empty()) { if (curr_l < grid_.nlevelmax) nsubcycle_[curr_l++] = std::stoi(s); }
    }
}

void Simulation::run() {
    int iout = 0, snapshot_count = 1;
    real_t courant = config_.get_double("hydro_params", "courant_factor", 0.8);
    
    while (t_ < tend_ && nstep_ < nstepmax_) {
        real_t dt = 1e10;
        // Search for min density to estimate dt
        real_t min_dt = 1e10;
        for (int il = 1; il <= grid_.nlevelmax; ++il) {
            int myid = 1; int ig = grid_.get_headl(myid, il);
            real_t dx = p::boxlen / (real_t)(p::nx * (1 << (il - 1)));
            while(ig > 0) {
                for(int ic=1; ic<=(1<<NDIM); ++ic) {
                    int idc = grid_.ncoarse + (ic-1)*grid_.ngridmax + ig - 1;
                    real_t d = std::max(grid_.uold(idc+1, 1), 1e-10);
                    real_t v2 = 0; for(int j=1; j<=NDIM; ++j) { real_t v = grid_.uold(idc+1, 1+j)/d; v2 += v*v; }
                    real_t p = (grid_.uold(idc+1, 5) - 0.5*d*v2)*(grid_.gamma-1.0); // Simplified
                    real_t c = std::sqrt(std::max(grid_.gamma*p/d, 1e-10));
                    real_t v_max = std::sqrt(v2) + c;
                    min_dt = std::min(min_dt, courant * dx / v_max);
                }
                ig = grid_.next[ig - 1];
            }
        }
        dt = min_dt;

        if (iout < (int)tout_.size() && t_ + dt > tout_[iout]) {
            dt = tout_[iout] - t_;
        }
        
        // Global refinement sweep
        real_t ed = config_.get_double("refine_params", "err_grad_d", 0.05);
        real_t ep = config_.get_double("refine_params", "err_grad_p", 0.0);
        real_t ev = config_.get_double("refine_params", "err_grad_u", 0.0);
        real_t eb2 = config_.get_double("refine_params", "err_grad_b2", ep);
        std::vector<real_t> evar;
        std::string evar_s = config_.get("refine_params", "err_grad_var", "");
        if (!evar_s.empty()) {
            std::stringstream ss(evar_s); std::string item;
            while (std::getline(ss, item, ',')) evar.push_back(std::stod(item));
        }
        int npassive = grid_.nvar - (11 + nener_);
#ifndef MHD
        npassive = grid_.nvar - (5 + nener_);
#endif
        if (evar.empty() && npassive > 0) evar.assign(npassive, 0.05);

        for (int il = 1; il < grid_.nlevelmax; ++il) {
             updater_.flag_fine(il, ed, ep, ev, eb2, evar);
             updater_.make_grid_fine(il);
             updater_.remove_grid_fine(il);
        }

        amr_step(1, dt); t_ += dt; nstep_++;
        if (iout < (int)tout_.size() && t_ >= tout_[iout]) {
            dump_snapshot(snapshot_count++); iout++;
        }
        if (nstep_ % ncontrol_ == 0) printf(" Step=%d t=%12.5e dt=%10.3e\n", nstep_, t_, dt);
    }
    dump_snapshot(snapshot_count);
}

void Simulation::amr_step(int ilevel, real_t dt, int icount) {
    if (ilevel > grid_.nlevelmax || grid_.count_grids_at_level(ilevel) == 0) return;
    
    real_t dx = p::boxlen / (real_t)(p::nx * (1 << (ilevel - 1)));
#ifdef MHD
    mhd_.godunov_fine(ilevel, dt, dx);
#else
    hydro_.godunov_fine(ilevel, dt, dx);
#endif
    int nsub = (ilevel <= (int)nsubcycle_.size()) ? nsubcycle_[ilevel - 1] : 1;
    for (int i = 1; i <= nsub; ++i) amr_step(ilevel + 1, dt / nsub, i);
    updater_.restrict_fine(ilevel);
    
#ifdef MHD
    mhd_.set_uold(ilevel);
#else
    hydro_.set_uold(ilevel);
#endif
}

void Simulation::dump_snapshot(int iout) {
    std::stringstream ssd; ssd << "output_" << std::setfill('0') << std::setw(5) << iout;
    std::string dir = ssd.str(); mkdir(dir.c_str(), 0777);
    SnapshotInfo info; info.t = t_; info.nstep = nstep_; info.noutput = 100; info.iout = iout; info.gamma = grid_.gamma; info.nener = nener_;
    int myid = 1;
    auto get_path = [&](const std::string& prefix, const std::string& ext) -> std::string {
        std::stringstream ss; ss << dir << "/" << prefix << "_" << std::setfill('0') << std::setw(5) << iout << ext << std::setfill('0') << std::setw(5) << myid;
        return ss.str();
    };
    RamsesWriter(get_path("amr", ".out")).write_amr(grid_, info);
    RamsesWriter(get_path("hydro", ".out")).write_hydro(grid_, info);
    
    std::stringstream ssinfo; ssinfo << dir << "/info_" << std::setfill('0') << std::setw(5) << iout << ".txt";
    RamsesWriter(ssinfo.str()).write_header(grid_, info);
    
    RamsesWriter(dir + "/hydro_file_descriptor.txt").write_hydro_descriptor(grid_, info);
    
    std::stringstream ssh; ssh << dir << "/header_" << std::setfill('0') << std::setw(5) << iout << ".txt";
    RamsesWriter(ssh.str()).write_header_file(grid_, info);
}

} // namespace ramses
