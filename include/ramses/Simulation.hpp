#ifndef RAMSES_SIMULATION_HPP
#define RAMSES_SIMULATION_HPP

#include "AmrGrid.hpp"
#include "HydroSolver.hpp"
#include "MhdSolver.hpp"
#include "PoissonSolver.hpp"
#include "TreeUpdater.hpp"
#include "Config.hpp"
#include "ParticleSystem.hpp"
#include "ParticleSolver.hpp"
#include "LoadBalancer.hpp"

namespace ramses {

/**
 * @brief Main simulation driver for RAMSES-CPP.
 */
class Simulation {
public:
    Simulation() : hydro_(grid_), mhd_(grid_, config_), poisson_(grid_, config_), updater_(grid_), ptcl_solver_(grid_, ps_), balancer_(grid_, ps_) {}

    void initialize(const std::string& nml_path);
    void run();

private:
    void amr_step(int ilevel, real_t dt);
    void dump_snapshot(int iout);

    AmrGrid grid_;
    HydroSolver hydro_;
    MhdSolver mhd_;
    PoissonSolver poisson_;
    TreeUpdater updater_;
    Config config_;
    ParticleSystem ps_;
    ParticleSolver ptcl_solver_;
    LoadBalancer balancer_;

    real_t t_ = 0.0;
    real_t tend_ = 1.0;
    int nstep_ = 0;
    int nstepmax_ = 10;
    int noutput_ = 0;
    int nener_ = 0;
    std::vector<double> tout_;
};

} // namespace ramses

#endif // RAMSES_SIMULATION_HPP
