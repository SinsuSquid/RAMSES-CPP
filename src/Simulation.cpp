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

namespace ramses {
namespace p = params;

void Simulation::initialize(const std::string& nml_path) {
    config_.parse(nml_path);
    p::nx = config_.get_int("amr_params", "nx", 1);
    p::ny = config_.get_int("amr_params", "ny", 1);
    p::nz = config_.get_int("amr_params", "nz", 1);
    p::boxlen = config_.get_double("amr_params", "boxlen", 1.0);
    int ngridmax = config_.get_int("amr_params", "ngridmax", 1000);
    p::ngridmax = ngridmax;
    nener_ = config_.get_int("hydro_params", "nener", 0);
#ifdef RAMSES_NENER
    if (nener_ == 0) nener_ = RAMSES_NENER;
#endif
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
    
    for (int il = 1; il < levelmax; ++il) {
        initializer_.apply_all();
        real_t ed = config_.get_double("refine_params", "err_grad_d", 0.05);
        real_t ep = config_.get_double("refine_params", "err_grad_p", 0.0);
        real_t ev = config_.get_double("refine_params", "err_grad_u", 0.0);
        real_t eb2 = config_.get_double("refine_params", "err_grad_b2", ep);
        updater_.flag_fine(il, ed, ep, ev, eb2);
        updater_.make_grid_fine(il);
    }
    nstep_ = 0; // Reset step counter
    initializer_.apply_all();

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
    int curr_l = 0;
    for (const std::string& s : nsub_s) {
        size_t star = s.find('*');
        if (star != std::string::npos) {
            int count = std::stoi(s.substr(0, star));
            int val = std::stoi(s.substr(star + 1));
            for (int i = 0; i < count && curr_l < grid_.nlevelmax; ++i) nsubcycle_[curr_l++] = val;
        } else if (!s.empty()) {
            if (curr_l < grid_.nlevelmax) nsubcycle_[curr_l++] = std::stoi(s);
        }
    }
}

void Simulation::run() {
    std::cout << "[Simulation] Starting main loop..." << std::endl;
    dump_snapshot(1);
    int iout = 0, snapshot_count = 2;
    while (t_ < tend_ && nstep_ < nstepmax_) {
        nstep_++; real_t dt = 1e30;
        for (int il = 1; il <= grid_.nlevelmax; ++il) {
            if (grid_.count_grids_at_level(il) == 0) continue;
            real_t dx = p::boxlen / (real_t)(p::nx * (1 << (il - 1)));
#ifdef MHD
            real_t ldt = mhd_.compute_courant_step(il, dx, grid_.gamma, 0.8);
#else
            real_t ldt = hydro_.compute_courant_step(il, dx, grid_.gamma, 0.8);
#endif
            dt = std::min(dt, ldt);
        }
        if (iout < (int)tout_.size() && t_ + dt >= tout_[iout]) {
            dt = tout_[iout] - t_;
        }
        
        // Global refinement sweep (Single Pass)
        real_t ed = config_.get_double("refine_params", "err_grad_d", 0.05);
        real_t ep = config_.get_double("refine_params", "err_grad_p", 0.0);
        real_t ev = config_.get_double("refine_params", "err_grad_u", 0.0);
        real_t eb2 = config_.get_double("refine_params", "err_grad_b2", ep);
        
        for (int il = 1; il < grid_.nlevelmax; ++il) updater_.flag_fine(il, ed, ep, ev, eb2);
        for (int il = 1; il < grid_.nlevelmax; ++il) {
            updater_.make_grid_fine(il);
            updater_.remove_grid_fine(il);
        }

        amr_step(1, dt); t_ += dt;
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
    RamsesWriter writer_amr(get_path("amr", ".out")); writer_amr.write_amr(grid_, info);
    RamsesWriter writer_hydro(get_path("hydro", ".out")); writer_hydro.write_hydro(grid_, info);
    std::stringstream ssi; ssi << dir << "/info_" << std::setfill('0') << std::setw(5) << iout << ".txt";
    RamsesWriter writer_header(ssi.str()); writer_header.write_header(grid_, info);
    std::stringstream ssh; ssh << dir << "/header_" << std::setfill('0') << std::setw(5) << iout << ".txt";
    RamsesWriter writer_header_file(ssh.str()); writer_header_file.write_header_file(grid_, info);
    RamsesWriter writer_desc(dir + "/hydro_file_descriptor.txt"); writer_desc.write_hydro_descriptor(grid_, info);
    std::ofstream fpd(dir + "/part_file_descriptor.txt"); fpd << "# total 0" << std::endl; fpd.close();
}

} // namespace ramses
