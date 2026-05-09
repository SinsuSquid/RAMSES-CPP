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

    // 1. Local accretion
    int n2d_val = (1 << NDIM);
    int myid = MpiManager::instance().rank() + 1;
    real_t r_acc = config_.get_double("sink_params", "r_acc", 0.0);

    auto process_cell = [&](int idc) {
        if (grid_.son.at(idc - 1) > 0) return;
        real_t xc[3];
        grid_.get_cell_center(idc, xc);
        
        for (size_t i = 0; i < nsink; ++i) {
            real_t d2 = 0;
            for (int d = 0; d < NDIM; ++d) d2 += std::pow(xc[d] - sinks_[i].x[d], 2);
            if (d2 < r_acc * r_acc) {
                real_t rho = grid_.uold(idc, 1);
                real_t dm = std::min(rho * 0.1 * dt, rho - 1e-10); // Cap accretion to keep density positive
                if (dm > 0) {
                    grid_.unew(idc, 1) -= dm;
                    dM[i] += dm;
                    for (int d = 0; d < NDIM; ++d) {
                        real_t mom = grid_.uold(idc, 2 + d);
                        real_t dp = mom * (dm / rho);
                        grid_.unew(idc, 2 + d) -= dp;
                        dP[i * 3 + d] += dp;
                    }
                    // Update total energy to reflect mass/momentum loss
                    int iener = NDIM + 2;
                    real_t v2 = 0; for(int d=0; d<NDIM; ++d) v2 += std::pow(grid_.uold(idc, 2+d)/rho, 2);
                    grid_.unew(idc, iener) -= 0.5 * dm * v2;
                }
            }
        }
    };

    if (ilevel == 1) {
        for (int idc = 1; idc <= grid_.ncoarse; ++idc) process_cell(idc);
    }
    int ig = grid_.get_headl(myid, ilevel);
    while (ig > 0) {
        for (int ic = 1; ic <= n2d_val; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
            process_cell(idc);
        }
        ig = grid_.next[ig - 1];
    }

#ifdef RAMSES_USE_MPI
    auto& mpi = MpiManager::instance();
    if (mpi.size() > 1) {
        std::vector<real_t> dM_global(nsink);
        std::vector<real_t> dP_global(nsink * 3);
        MPI_Allreduce(dM.data(), dM_global.data(), nsink, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(dP.data(), dP_global.data(), nsink * 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dM = dM_global; dP = dP_global;
    }
#endif

    for (size_t i = 0; i < nsink; ++i) {
        sinks_[i].mass += dM[i];
        if (sinks_[i].mass > 0) {
            for (int d = 0; d < NDIM; ++d) {
                sinks_[i].v[d] = (sinks_[i].v[d] * (sinks_[i].mass - dM[i]) + dP[i * 3 + d]) / sinks_[i].mass;
            }
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
