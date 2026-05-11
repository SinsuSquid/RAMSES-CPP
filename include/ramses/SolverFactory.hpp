#ifndef RAMSES_SOLVER_FACTORY_HPP
#define RAMSES_SOLVER_FACTORY_HPP

#include "HydroSolver.hpp"
#include "RhdSolver.hpp"
#include "CoolingSolver.hpp"
#include "TurbulenceSolver.hpp"
#include "SinkSolver.hpp"
#include "MhdSolver.hpp"
#include "RtSolver.hpp"
#include "PoissonSolver.hpp"
#include "ParticleSolver.hpp"
#include "Initializer.hpp"
#include <memory>

#include "StarSolver.hpp"
#include "FeedbackSolver.hpp"
#include "ClumpFinder.hpp"
#include "LightCone.hpp"

namespace ramses {

std::unique_ptr<HydroSolver> create_hydro_solver(AmrGrid& grid, Config& config);
std::unique_ptr<RhdSolver> create_rhd_solver(AmrGrid& grid, Config& config);
std::unique_ptr<TurbulenceSolver> create_turbulence_solver(AmrGrid& grid, Config& config);
std::unique_ptr<SinkSolver> create_sink_solver(AmrGrid& grid, Config& config);
std::unique_ptr<StarSolver> create_star_solver(AmrGrid& grid, Config& config);
std::unique_ptr<FeedbackSolver> create_feedback_solver(AmrGrid& grid, Config& config);
std::unique_ptr<ClumpFinder> create_clump_finder(AmrGrid& grid, Config& config);
std::unique_ptr<LightCone> create_light_cone(AmrGrid& grid, Config& config);
std::unique_ptr<CoolingSolver> create_cooling_solver(AmrGrid& grid, Config& config);
std::unique_ptr<MhdSolver> create_mhd_solver(AmrGrid& grid, Config& config);
std::unique_ptr<RtSolver> create_rt_solver(AmrGrid& grid, Config& config);
std::unique_ptr<PoissonSolver> create_poisson_solver(AmrGrid& grid, Config& config);
std::unique_ptr<ParticleSolver> create_particle_solver(AmrGrid& grid, Config& config);
std::unique_ptr<Initializer> create_initializer(AmrGrid& grid, Config& config);

} // namespace ramses

#endif
