#ifndef RAMSES_SIMULATION_HPP
#define RAMSES_SIMULATION_HPP

#include "AmrGrid.hpp"
#include "HydroSolver.hpp"
#include "MhdSolver.hpp"
#include "PoissonSolver.hpp"
#include "TreeUpdater.hpp"
#include "Config.hpp"
#include "Initializer.hpp"
#include "RtSolver.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Main simulation driver for RAMSES-CPP.
 */
class Simulation {
public:
    Simulation() : grid_(), 
                   hydro_(grid_, config_), 
                   mhd_(grid_, config_), 
                   rt_(grid_, config_),
                   poisson_(grid_, config_), 
                   updater_(grid_, config_), 
                   initializer_(grid_, config_) {}

    void initialize(const std::string& nml_path);
    void run();

private:
    void amr_step(int ilevel, real_t dt, int icount = 1);
    void dump_snapshot(int iout);
    
    Config config_;
    AmrGrid grid_;
    HydroSolver hydro_;
    MhdSolver mhd_;
    RtSolver rt_;
    PoissonSolver poisson_;
    TreeUpdater updater_;
    Initializer initializer_;

    real_t t_ = 0.0;
    real_t tend_ = 1.0;
    int nstep_ = 0;
    int nstepmax_ = 10;
    int ncontrol_ = 1;
    int noutput_ = 0;
    int nener_ = 0;
    std::vector<int> nsubcycle_;
    std::vector<int> nexpand_;
    std::vector<double> tout_;
};

} // namespace ramses

#endif // RAMSES_SIMULATION_HPP
