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
    int nvar = 5 + nener_;
#ifdef MHD
    nvar = 8 + nener_;
#endif
    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar, ncpu, levelmax);
    updater_.set_interpol_hook([this](const real_t u1[7][20], real_t u2[8][20]){
#ifdef MHD
        mhd_.interpol_mhd(u1, u2);
#else
        hydro_.interpol_hydro(u1, u2);
#endif
    });
    
    for (int il = 1; il < levelmax; ++il) {
        initializer_.apply_all();
        if (il < levelmin) {
            if (il == 1) { for(int i=0; i<grid_.ncoarse; ++i) grid_.flag1[i] = 1; }
            else {
                int myid = 1;
                int ig = grid_.get_headl(myid, il);
                while(ig > 0) {
                    for(int ic=1; ic<=constants::twotondim; ++ic) {
                        int idc = grid_.ncoarse + (ic-1)*grid_.ngridmax + ig - 1;
                        grid_.flag1[idc] = 1;
                    }
                    ig = grid_.next[ig-1];
                }
            }
        } else {
            updater_.flag_fine(il, config_.get_double("refine_params", "err_grad_d", 0.05), 0, 0);
        }
        updater_.make_grid_fine(il);
    }
    initializer_.apply_all();

    grid_.nboundary = config_.get_int("boundary_params", "nboundary", 0);
    if (grid_.nboundary > 0) {
        auto parse_vec = [](const std::string& s) {
            std::vector<int> v; std::string sc = s; std::replace(sc.begin(), sc.end(), ',', ' ');
            std::stringstream ss(sc); int val; while (ss >> val) v.push_back(val); return v;
        };
        grid_.ibound_min = parse_vec(config_.get("boundary_params", "ibound_min", ""));
        grid_.ibound_max = parse_vec(config_.get("boundary_params", "ibound_max", ""));
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
}

void Simulation::run() {
    std::cout << "[Simulation] Starting main loop..." << std::endl;
    dump_snapshot(1);
    int iout = 0, snapshot_count = 2;
    while (t_ < tend_ && nstep_ < nstepmax_) {
        nstep_++; real_t dt = 1e30;
        for (int il = p::levelmin; il <= grid_.nlevelmax; ++il) {
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
            dt = tout_[iout] - t_; amr_step(p::levelmin, dt); t_ += dt;
            dump_snapshot(snapshot_count++); iout++;
        } else {
            amr_step(p::levelmin, dt); t_ += dt;
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
    int nsub = 1; 
    for (int i = 1; i <= nsub; ++i) amr_step(ilevel + 1, dt / nsub, i);
    updater_.restrict_fine(ilevel);
    if (icount == nsub) {
        updater_.flag_fine(ilevel, 0.05, 0, 0);
        updater_.make_grid_fine(ilevel);
    }
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
    RamsesWriter writer_header(ssi.str()); writer_header.write_header(info);
    std::stringstream ssh; ssh << dir << "/header_" << std::setfill('0') << std::setw(5) << iout << ".txt";
    RamsesWriter writer_header_file(ssh.str()); writer_header_file.write_header_file(grid_, info);
    RamsesWriter writer_desc(dir + "/hydro_file_descriptor.txt"); writer_desc.write_hydro_descriptor(grid_, info);
    std::ofstream fpd(dir + "/part_file_descriptor.txt"); fpd << "# total 0" << std::endl; fpd.close();
}

} // namespace ramses
