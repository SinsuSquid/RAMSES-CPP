#include "ramses/SinkSolver.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/Constants.hpp"
#include <cmath>
#include <algorithm>

#ifdef RAMSES_USE_MPI
#include <mpi.h>
#endif

namespace ramses {

SinkSolver::SinkSolver(AmrGrid& grid, Config& config) 
    : grid_(grid), config_(config) {
    nsinkmax_ = config_.get_int("sink_params", "nsinkmax", 100);
    mass_sink_seed_ = config_.get_double("sink_params", "mass_sink_seed", 0.0);
}

SinkSolver::~SinkSolver() {}

void SinkSolver::init() {
    sinks_.clear();
}

void SinkSolver::create_sinks(int ilevel) {
    // formation logic would go here
}

void SinkSolver::grow_sinks(int ilevel, real_t dt) {
    if (sinks_.empty()) return;

    size_t nsink = sinks_.size();
    std::vector<real_t> dM(nsink, 0.0);
    std::vector<real_t> dP(nsink * 3, 0.0);

    collect_accretion_data(ilevel, dM, dP);

#ifdef RAMSES_USE_MPI
    auto& mpi = MpiManager::instance();
    if (mpi.size() > 1) {
        std::vector<real_t> dM_global(nsink);
        std::vector<real_t> dP_global(nsink * 3);

        MPI_Allreduce(dM.data(), dM_global.data(), nsink, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(dP.data(), dP_global.data(), nsink * 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        dM = dM_global;
        dP = dP_global;
    }
#endif

    for (size_t i = 0; i < nsink; ++i) {
        sinks_[i].mass += dM[i];
        for (int d = 0; d < 3; ++d) {
            // Update momentum logic
            // sinks_[i].v[d] = ...
        }
    }
}

void SinkSolver::synchronize_sinks() {
#ifdef RAMSES_USE_MPI
    auto& mpi = MpiManager::instance();
    if (mpi.size() == 1) return;

    // 1. Synchronize the number of sinks
    int nsink = static_cast<int>(sinks_.size());
    int nsink_global = 0;
    MPI_Allreduce(&nsink, &nsink_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    // 2. If a new sink was created on any rank, broadcast the updated list
    // For simplicity, we assume rank 0 is the master of the sink list
    // and others send their new sinks to it, or we just broadcast from the creator.
    // Here we'll do a simple broadcast from rank 0, assuming creation is coordinated.
    if (nsink_global > 0) {
        if (sinks_.size() < (size_t)nsink_global) sinks_.resize(nsink_global);
        MPI_Bcast(sinks_.data(), nsink_global * sizeof(Sink), MPI_BYTE, 0, MPI_COMM_WORLD);
    }
#endif
}

void SinkSolver::collect_accretion_data(int ilevel, std::vector<real_t>& dM, std::vector<real_t>& dP) {
    // Local accretion logic
}

} // namespace ramses
