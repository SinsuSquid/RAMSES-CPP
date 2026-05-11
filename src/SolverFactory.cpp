#include "ramses/SolverFactory.hpp"

namespace ramses {

#ifndef RAMSES_HAS_PATCH_HYDRO
std::unique_ptr<HydroSolver> create_hydro_solver(AmrGrid& grid, Config& config) {
    return std::make_unique<HydroSolver>(grid, config);
}
#endif

#ifndef RAMSES_HAS_PATCH_RHD
std::unique_ptr<RhdSolver> create_rhd_solver(AmrGrid& grid, Config& config) {
    return std::make_unique<RhdSolver>(grid, config);
}
#endif

#ifndef RAMSES_HAS_PATCH_TURBULENCE
std::unique_ptr<TurbulenceSolver> create_turbulence_solver(AmrGrid& grid, Config& config) {
    return std::make_unique<TurbulenceSolver>(grid, config);
}
#endif

#ifndef RAMSES_HAS_PATCH_SINK
std::unique_ptr<SinkSolver> create_sink_solver(AmrGrid& grid, Config& config) {
    return std::make_unique<SinkSolver>(grid, config);
}
#endif

#ifndef RAMSES_HAS_PATCH_STAR
std::unique_ptr<StarSolver> create_star_solver(AmrGrid& grid, Config& config) {
    return std::make_unique<StarSolver>(grid, config);
}
#endif

#ifndef RAMSES_HAS_PATCH_COOLING
std::unique_ptr<CoolingSolver> create_cooling_solver(AmrGrid& grid, Config& config) {
    return std::make_unique<CoolingSolver>(grid, config);
}
#endif

#ifndef RAMSES_HAS_PATCH_MHD
std::unique_ptr<MhdSolver> create_mhd_solver(AmrGrid& grid, Config& config) {
    return std::make_unique<MhdSolver>(grid, config);
}
#endif

#ifndef RAMSES_HAS_PATCH_RT
std::unique_ptr<RtSolver> create_rt_solver(AmrGrid& grid, Config& config) {
    return std::make_unique<RtSolver>(grid, config);
}
#endif

#ifndef RAMSES_HAS_PATCH_POISSON
std::unique_ptr<PoissonSolver> create_poisson_solver(AmrGrid& grid, Config& config) {
    return std::make_unique<PoissonSolver>(grid, config);
}
#endif

#ifndef RAMSES_HAS_PATCH_PARTICLE
std::unique_ptr<ParticleSolver> create_particle_solver(AmrGrid& grid, Config& config) {
    return std::make_unique<ParticleSolver>(grid, config);
}
#endif

#ifndef RAMSES_HAS_PATCH_INITIALIZER
std::unique_ptr<Initializer> create_initializer(AmrGrid& grid, Config& config) {
    return std::make_unique<Initializer>(grid, config);
}
#endif

} // namespace ramses
