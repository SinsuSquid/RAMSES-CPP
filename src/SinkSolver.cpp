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
    if (!config_.get_bool("sink_params", "create_sinks", false)) return;
    
    // In a real port, we'd call the clump finder here.
    // For now, we'll implement a simple density-based seed formation if a clump finder isn't available.
    real_t rho_crit = config_.get_double("sink_params", "rho_crit", 1e6);
    int n2d_val = (1 << NDIM);
    int myid = MpiManager::instance().rank() + 1;
    real_t boxlen = config_.get_double("amr_params", "boxlen", 1.0);

    auto check_and_create = [&](int idc) {
        if (grid_.son[idc - 1] > 0) return;
        if (grid_.uold(idc, 1) > rho_crit) {
            real_t xc[3];
            grid_.get_cell_center(idc, xc);
            
            // Check if there is already a sink nearby
            real_t dx = boxlen / (1 << (ilevel - 1));
            bool too_close = false;
            for (const auto& s : sinks_) {
                real_t d2 = 0;
                for (int d=0; d<NDIM; ++d) d2 += std::pow(xc[d] - s.x[d], 2);
                if (d2 < 4 * dx * dx) { too_close = true; break; }
            }

            if (!too_close) {
                Sink new_sink;
                new_sink.id = static_cast<int>(sinks_.size() + 1);
                for(int d=0; d<NDIM; ++d) {
                    new_sink.x[d] = xc[d];
                    new_sink.v[d] = grid_.uold(idc, 2+d) / grid_.uold(idc, 1);
                }
                new_sink.mass = mass_sink_seed_;
                new_sink.level = ilevel;
                new_sink.t_creation = 0.0; // Should use simulation time
                sinks_.push_back(new_sink);
            }
        }
    };

    if (ilevel == 1) {
        for (int idc = 1; idc <= grid_.ncoarse; ++idc) check_and_create(idc);
    }
    int ig = grid_.get_headl(myid, ilevel);
    while (ig > 0) {
        for (int ic = 1; ic <= n2d_val; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
            check_and_create(idc);
        }
        ig = grid_.next[ig - 1];
    }
}

void SinkSolver::grow_sinks(int ilevel, real_t dt) {
    if (sinks_.empty()) return;

    size_t nsink = sinks_.size();
    std::vector<real_t> dM(nsink, 0.0);
    std::vector<real_t> dP(nsink * 3, 0.0);
    std::vector<real_t> dL(nsink * 3, 0.0);

    // 1. Local accretion
    int n2d_val = (1 << NDIM);
    int myid = MpiManager::instance().rank() + 1;
    real_t r_acc = config_.get_double("sink_params", "r_acc", 0.1);

    auto process_cell = [&](int idc) {
        if (grid_.son[idc - 1] > 0) return;
        real_t xc[3];
        grid_.get_cell_center(idc, xc);
        
        for (size_t i = 0; i < nsink; ++i) {
            if (sinks_[i].merged) continue;
            real_t d2 = 0;
            for (int d = 0; d < NDIM; ++d) d2 += std::pow(xc[d] - sinks_[i].x[d], 2);
            if (d2 < r_acc * r_acc) {
                real_t rho = grid_.uold(idc, 1);
                real_t dm = std::min(rho * 0.1 * dt, rho - 1e-10); 
                if (dm > 0) {
                    grid_.unew(idc, 1) -= dm;
                    dM[i] += dm;
                    for (int d = 0; d < NDIM; ++d) {
                        real_t mom = grid_.uold(idc, 2 + d);
                        real_t dp = mom * (dm / rho);
                        grid_.unew(idc, 2 + d) -= dp;
                        dP[i * 3 + d] += dp;
                        
                        // Angular momentum: L = r x p
                        int d1 = (d + 1) % 3;
                        int d2_idx = (d + 2) % 3;
                        dL[i * 3 + d] += (xc[d1] - sinks_[i].x[d1]) * dp; // This is a simplification
                    }
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
        std::vector<real_t> dL_global(nsink * 3);
        MPI_Allreduce(dM.data(), dM_global.data(), nsink, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(dP.data(), dP_global.data(), nsink * 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(dL.data(), dL_global.data(), nsink * 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dM = dM_global; dP = dP_global; dL = dL_global;
    }
#endif

    for (size_t i = 0; i < nsink; ++i) {
        if (sinks_[i].merged) continue;
        real_t m_old = sinks_[i].mass;
        sinks_[i].mass += dM[i];
        if (sinks_[i].mass > 0) {
            for (int d = 0; d < NDIM; ++d) {
                sinks_[i].v[d] = (sinks_[i].v[d] * m_old + dP[i * 3 + d]) / sinks_[i].mass;
                sinks_[i].ang_mom[d] += dL[i * 3 + d];
            }
        }
    }
}

void SinkSolver::merge_sinks() {
    if (sinks_.size() < 2) return;
    real_t boxlen = config_.get_double("amr_params", "boxlen", 1.0);
    real_t dx_min = boxlen / (1 << (grid_.nlevelmax - 1));
    real_t merge_dist2 = 4 * dx_min * dx_min;

    for (size_t i = 0; i < sinks_.size(); ++i) {
        if (sinks_[i].merged) continue;
        for (size_t j = i + 1; j < sinks_.size(); ++j) {
            if (sinks_[j].merged) continue;
            
            real_t d2 = 0;
            for (int d = 0; d < NDIM; ++d) d2 += std::pow(sinks_[i].x[d] - sinks_[j].x[d], 2);
            
            if (d2 < merge_dist2) {
                // Merge j into i
                real_t m_tot = sinks_[i].mass + sinks_[j].mass;
                if (m_tot > 0) {
                    for (int d = 0; d < NDIM; ++d) {
                        sinks_[i].x[d] = (sinks_[i].x[d] * sinks_[i].mass + sinks_[j].x[d] * sinks_[j].mass) / m_tot;
                        sinks_[i].v[d] = (sinks_[i].v[d] * sinks_[i].mass + sinks_[j].v[d] * sinks_[j].mass) / m_tot;
                        sinks_[i].ang_mom[d] += sinks_[j].ang_mom[d];
                    }
                }
                sinks_[i].mass = m_tot;
                sinks_[j].merged = true;
            }
        }
    }

    // Remove merged sinks
    sinks_.erase(std::remove_if(sinks_.begin(), sinks_.end(), [](const Sink& s) { return s.merged; }), sinks_.end());
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
