#ifndef RAMSES_SIMULATION_HPP
#define RAMSES_SIMULATION_HPP

#include "AmrGrid.hpp"
#include "HydroSolver.hpp"
#include "CoolingSolver.hpp"
#include "TurbulenceSolver.hpp"
#include "MhdSolver.hpp"
#include "PoissonSolver.hpp"
#include "TreeUpdater.hpp"
#include "Config.hpp"
#include "Initializer.hpp"
#include "RtSolver.hpp"
#include "Cosmology.hpp"
#include "ParticleSolver.hpp"
#include "LoadBalancer.hpp"
#include <vector>
#include <memory>

namespace ramses {

/**
 * @brief Main simulation driver for RAMSES-CPP.
 */
class Simulation {
public:
    Simulation();

    void initialize(const std::string& nml_path);
    void run();

private:
    void amr_step(int ilevel, real_t dt, int icount = 1);
    void dump_snapshot(int iout);
    void rho_fine(int ilevel);

    Config config_;
    AmrGrid grid_;
    std::unique_ptr<HydroSolver> hydro_;
    std::unique_ptr<CoolingSolver> cooling_;
    std::unique_ptr<TurbulenceSolver> turb_;
    std::unique_ptr<MhdSolver> mhd_;
    std::unique_ptr<RtSolver> rt_;
    std::unique_ptr<PoissonSolver> poisson_;
    TreeUpdater updater_;
    std::unique_ptr<Initializer> initializer_;
    std::unique_ptr<ParticleSolver> particles_;
    LoadBalancer load_balancer_;
    Cosmology cosmo_;

    real_t t_ = 0.0;
    real_t tend_ = 1.0;
    real_t aexp_ = 1.0;
    real_t hexp_ = 0.0;

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
