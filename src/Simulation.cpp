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

namespace ramses {
namespace p = params;

void Simulation::initialize(const std::string& nml_path) {
    config_.parse(nml_path);
    
    p::nx = config_.get_int("amr_params", "nx", 1);
    p::ny = config_.get_int("amr_params", "ny", 1);
    p::nz = config_.get_int("amr_params", "nz", 1);
    
    int ngridmax = config_.get_int("amr_params", "ngridmax", 0);
    if (ngridmax == 0) ngridmax = config_.get_int("amr_params", "ngridtot", 1000000);
    
    nener_ = config_.get_int("hydro_params", "nener", 0);
    if (nener_ == 0) {
        for (int i = 1; i <= 10; ++i) {
            if (!config_.get("init_params", "prad_region(" + std::to_string(i) + ",1)", "").empty())
                nener_ = std::max(nener_, i);
        }
    }
    
    int nvar_default = 8 + nener_;
#ifdef MHD
    nvar_default = 8 + nener_;
#endif
#ifdef RT
    int nGroups = config_.get_int("rt_params", "nGroups", 0);
    nvar_default += nGroups * (1 + NDIM);
#endif
    int nvar = config_.get_int("hydro_params", "nvar", nvar_default);
    
    int ncpu = MpiManager::instance().size();
#ifdef MHD
    int nvar_all = nvar + 3;
    int nlevelmax = config_.get_int("amr_params", "levelmax", 10);
    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar_all, ncpu, nlevelmax);
#else
    int nlevelmax = config_.get_int("amr_params", "levelmax", 10);
    grid_.allocate(p::nx, p::ny, p::nz, ngridmax, nvar, ncpu, nlevelmax);
#endif
    
    noutput_ = config_.get_int("output_params", "noutput", 0);
    ncontrol_ = config_.get_int("run_params", "ncontrol", 1);

    p::units_length = config_.get_double("units_params", "units_length", 1.0);
    p::units_density = config_.get_double("units_params", "units_density", 1.0);
    p::units_time = config_.get_double("units_params", "units_time", 1.0);
    p::units_velocity = p::units_length / p::units_time;
    p::units_mass = p::units_density * std::pow(p::units_length, 3);
    p::units_energy = p::units_mass * std::pow(p::units_velocity, 2);
    p::units_pressure = p::units_density * std::pow(p::units_velocity, 2);
    p::units_number_density = p::units_density / (1.67e-24);
    
    grid_.gamma = config_.get_double("hydro_params", "gamma", 1.4);
    grid_.err_grad_d = config_.get_double("refine_params", "err_grad_d", 0.05);
    grid_.err_grad_u = config_.get_double("refine_params", "err_grad_u", 0.05);
    grid_.err_grad_p = config_.get_double("refine_params", "err_grad_p", 0.05);

    std::string tout_str = config_.get("output_params", "tout", "");
    if (!tout_str.empty()) {
        std::stringstream ss(tout_str); double val;
        while (ss >> val) { tout_.push_back(val); if (ss.peek() == ',' || ss.peek() == ' ') ss.ignore(); }
    }
    tend_ = config_.get_double("output_params", "tend", 1e10);
    if (!tout_.empty()) tend_ = tout_.back();

    // 1. Initial partition of coarse grid
    balancer_.calculate_hilbert_keys();
    balancer_.balance();

    // 2. Initial coarse level refinement
    updater_.mark_all(0);
    updater_.refine_coarse();
    
    // 3. Dynamic partition of Level 1
    balancer_.calculate_hilbert_keys();
    balancer_.balance();

    // 4. Partitioned refinement up to levelmin
    int levelmin = config_.get_int("amr_params", "levelmin", 1);
    for (int il = 1; il < levelmin; ++il) {
        updater_.mark_all(il);
        updater_.refine_fine(il);
    }
    
    Initializer init(grid_, config_);
    init.apply_all();
    rt_.initialize();
}

void Simulation::run() {
    std::cout << "[Simulation] Starting main loop..." << std::endl;
    real_t gamma = grid_.gamma;
    dump_snapshot(1);
    int iout = 0, snapshot_count = 2;

    while (t_ < tend_) {
        if (nstep_ % ncontrol_ == 0) {
            std::cout << "Mesh structure" << std::endl;
            for (int i = 1; i <= grid_.nlevelmax; ++i) {
                int ngrids = grid_.count_grids_at_level(i);
                if (ngrids > 0) printf(" Level %2d has %10d grids\n", i, ngrids);
            }
            printf(" Main step=%7d t=%12.5e dt=%10.3e\n", nstep_, t_, 0.0);
        }
        nstep_++;
        real_t dx = p::boxlen / static_cast<real_t>(p::nx), dt;
#ifdef MHD
        dt = mhd_.compute_courant_step(1, dx, gamma, config_.get_double("hydro_params", "courant_factor", 0.8));
#else
        dt = hydro_.compute_courant_step(1, dx, gamma, config_.get_double("hydro_params", "courant_factor", 0.8));
#endif
        if (iout < (int)tout_.size()) {
            if (t_ + dt >= tout_[iout]) {
                dt = tout_[iout] - t_;
                amr_step(1, dt);
                t_ += dt;
                dump_snapshot(snapshot_count++);
                iout++;
                continue;
            }
        }
        if (t_ + dt > tend_) dt = tend_ - t_;
        amr_step(1, dt);
        t_ += dt;

        // Dynamic Load Balancing
        if (MpiManager::instance().size() > 1) {
            balancer_.calculate_hilbert_keys();
            balancer_.balance();
        }

        if (nstep_ >= nstepmax_) break;
    }
    if (iout < (int)tout_.size() || t_ >= tend_) dump_snapshot(snapshot_count);
}

void Simulation::amr_step(int ilevel, real_t dt) {
    if (ilevel > grid_.nlevelmax) return;
    if (grid_.count_grids_at_level(ilevel) == 0) return;
    real_t dx = p::boxlen / static_cast<real_t>(p::nx * (1 << (ilevel - 1)));

#ifdef MHD
    mhd_.set_unew(ilevel);
#else
    hydro_.set_unew(ilevel);
#endif
    rt_.set_unew(ilevel);

    if (ilevel < grid_.nlevelmax) {
        if (grid_.count_grids_at_level(ilevel + 1) > 0) {
            amr_step(ilevel + 1, dt / 2.0);
            amr_step(ilevel + 1, dt / 2.0);
        }
    }

    if (config_.get_bool("run_params", "poisson", false)) {
        poisson_.compute_force(ilevel);
        hydro_.add_gravity_source_terms(ilevel, dt);
    }
    
#ifdef MHD
    mhd_.godunov_fine(ilevel, dt, dx);
#else
    hydro_.godunov_fine(ilevel, dt, dx);
#endif
    cooling_.apply_cooling(ilevel, dt);
    rt_.godunov_fine(ilevel, dt, dx);
    rt_.apply_source_terms(ilevel, dt);

#ifdef MHD
    mhd_.set_uold(ilevel);
#else
    hydro_.set_uold(ilevel);
#endif
    rt_.set_uold(ilevel);

    if (ilevel < grid_.nlevelmax) {
        updater_.mark_cells(ilevel);
        updater_.refine_fine(ilevel);
    }

    real_t mind, maxv, maxdb = 0;
#ifdef MHD
    mhd_.get_diagnostics(ilevel, dx, mind, maxv, maxdb);
    printf(" Fine step=%7d t=%12.5e dt=%10.3e level=%2d min_d=%9.2e max_v=%9.2e max_div_b=%9.2e\n", nstep_, t_ + dt, dt, ilevel, mind, maxv, maxdb);
#else
    hydro_.get_diagnostics(ilevel, dx, mind, maxv);
    printf(" Fine step=%7d t=%12.5e dt=%10.3e level=%2d min_d=%9.2e max_v=%9.2e\n", nstep_, t_ + dt, dt, ilevel, mind, maxv);
#endif
}

void Simulation::dump_snapshot(int iout) {
    std::stringstream ss; ss << "output_" << std::setw(5) << std::setfill('0') << iout;
    std::string dir = ss.str(); mkdir(dir.c_str(), 0777);
    std::string nchar = ss.str().substr(7);
    
    SnapshotInfo info;
    info.t = t_; info.nstep = nstep_; info.noutput = noutput_;
    info.tout = tout_; info.iout = iout; info.gamma = grid_.gamma;
    info.nener = nener_;

    RamsesWriter amr_w(dir + "/amr_" + nchar + ".out00001"); amr_w.write_amr(grid_, info);
    RamsesWriter hydro_w(dir + "/hydro_" + nchar + ".out00001"); hydro_w.write_hydro(grid_, info);
    if (config_.get_bool("run_params", "poisson", false)) {
        RamsesWriter grav_w(dir + "/grav_" + nchar + ".out00001"); grav_w.write_grav(grid_, info);
    }

    std::ofstream infof(dir + "/info_" + nchar + ".txt");
    if (infof.is_open()) {
        infof << std::scientific << std::setprecision(15);
        infof << "ncpu         = " << grid_.ncpu << "\nndim         = " << NDIM << "\nnx           = " << p::nx << "\nny           = " << p::ny << "\nnz           = " << p::nz << "\n";
        infof << "levelmin     = 1\nlevelmax     = " << grid_.nlevelmax << "\nngridmax     = " << grid_.ngridmax << "\nboxlen       = " << p::boxlen << "\ntime         = " << t_ << "\n";
        infof << "unit_l       = " << p::units_length << "\nunit_d       = " << p::units_density << "\nunit_t       = " << p::units_time << "\n";
        infof.close();
    }

    std::ofstream desc(dir + "/hydro_file_descriptor.txt");
    if (desc.is_open()) {
        int iv = 1; desc << iv++ << ", density, double\n";
        for (int i = 1; i <= NDIM; ++i) desc << iv++ << ", velocity_" << (i==1?'x':(i==2?'y':'z')) << ", double\n";
#ifdef MHD
        desc << iv++ << ", B_x_left, double\n" << iv++ << ", B_y_left, double\n" << iv++ << ", B_z_left, double\n";
        desc << iv++ << ", B_x_right, double\n" << iv++ << ", B_y_right, double\n" << iv++ << ", B_z_right, double\n";
        for (int i = 1; i <= nener_; ++i) desc << iv++ << ", non_thermal_pressure_" << std::setw(2) << std::setfill('0') << i << ", double\n";
        desc << iv++ << ", pressure, double\n";
        int nvp = grid_.nvar; int npscal = (nvp > 8 + nener_) ? nvp - (8 + nener_) : 0;
        for (int i = 1; i <= npscal; ++i) desc << iv++ << ", scalar_" << std::setw(2) << std::setfill('0') << i << ", double\n";
#else
        desc << iv++ << ", pressure, double\n";
        for (int i = 1; i <= nener_; ++i) desc << iv++ << ", non_thermal_pressure_" << std::setw(2) << std::setfill('0') << i << ", double\n";
        int nvp = grid_.nvar; int npscal = (nvp > 5 + nener_) ? nvp - (5 + nener_) : 0;
        if (npscal < 0) npscal = 0; // Guard
        for (int i = 1; i <= npscal; ++i) desc << iv++ << ", scalar_" << std::setw(2) << std::setfill('0') << i << ", double\n";
#endif
        desc.close();
    }
}

} // namespace ramses
