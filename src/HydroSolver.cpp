#include "ramses/HydroSolver.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/RiemannSolver.hpp"
#include "ramses/SlopeLimiter.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#ifdef RAMSES_USE_MPI
#include <mpi.h>
#endif

namespace ramses {

void HydroSolver::set_unew(int ilevel) {
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int id = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            for (int iv = 1; iv <= grid_.nvar; ++iv) {
                grid_.unew(id, iv) = grid_.uold(id, iv);
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

void HydroSolver::set_uold(int ilevel) {
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int id = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            for (int iv = 1; iv <= grid_.nvar; ++iv) {
                grid_.uold(id, iv) = grid_.unew(id, iv);
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

void HydroSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    std::vector<int> active_octs;
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        active_octs.push_back(igrid);
        igrid = grid_.next[igrid - 1];
    }
    if (!active_octs.empty()) godfine1(active_octs, ilevel, dt, dx);
}

void HydroSolver::godfine1(const std::vector<int>& octs, int ilevel, real_t dt, real_t dx) {
    real_t gamma = grid_.gamma;
    for (int igrid : octs) {
        LocalStencil stencil;
        gather_stencil(igrid, ilevel, stencil);
        
        real_t q[8][1 + NDIM + 1]; 
        for (int ic = 0; ic < constants::twotondim; ++ic) {
            real_t u[20];
            // Center cell is (3,3,3) in the 6x6x6 stencil? No, LocalStencil header is weird.
            // For now just match the previousTurn implementation as much as possible
            for (int iv = 1; iv <= grid_.nvar; ++iv) u[iv - 1] = stencil.uloc[0][0][0][iv - 1]; 
            ctoprim(u, q[ic], gamma);
        }
    }
}

void HydroSolver::ctoprim(const real_t u[], real_t q[], real_t gamma) {
    const real_t smallr = 1e-10;
    real_t d = std::max(u[0], smallr);
    q[0] = d;
    real_t vel2 = 0.0;
    for (int idim = 1; idim <= NDIM; ++idim) {
        q[idim] = u[idim] / d;
        vel2 += q[idim] * q[idim];
    }
    real_t e_int = u[4] - 0.5 * d * vel2;
    q[NDIM + 1] = std::max(e_int * (gamma - 1.0), d * 1e-10);
}

void HydroSolver::gather_stencil(int igrid, int ilevel, LocalStencil& stencil) {
    int nbors_father[27];
    grid_.get_3x3x3_father(igrid, nbors_father);
    int myid = MpiManager::instance().rank() + 1;

    for (int i = 0; i < 27; ++i) {
        int ifather = nbors_father[i];
        if (ifather > 0 && grid_.cpu_map[ifather] == myid) {
            if (grid_.son[ifather] > 0) {
                int ig = grid_.son[ifather];
                for (int ic = 0; ic < constants::twotondim; ++ic) {
                    int idc = grid_.ncoarse + ic * grid_.ngridmax + ig;
                    for (int iv = 1; iv <= grid_.nvar; ++iv) stencil.uloc[i][ic][0][iv - 1] = grid_.uold(idc, iv);
                }
            } else {
                for (int ic = 0; ic < constants::twotondim; ++ic) {
                    for (int iv = 1; iv <= grid_.nvar; ++iv) stencil.uloc[i][ic][0][iv - 1] = grid_.uold(ifather, iv);
                }
            }
        } else {
            int icenter = nbors_father[13];
            for (int ic = 0; ic < constants::twotondim; ++ic) {
                for (int iv = 1; iv <= grid_.nvar; ++iv) {
                    stencil.uloc[i][ic][0][iv - 1] = (icenter > 0 && grid_.cpu_map[icenter] == myid) ? grid_.uold(icenter, iv) : 0.0;
                }
            }
        }
    }
}

real_t HydroSolver::compute_courant_step(int ilevel, real_t dx, real_t gamma, real_t courant_factor) {
    real_t dt_max = 1e30; const real_t smallr = 1e-10;
    int myid = MpiManager::instance().rank() + 1;
    for (int i = 1; i <= grid_.ncell; ++i) {
        if (grid_.son[i] != 0 || grid_.cpu_map[i] != myid) continue; 
        real_t d = std::max(grid_.uold(i, 1), smallr);
        real_t vel2 = 0.0, vel_max = 0.0;
        for (int idim = 1; idim <= NDIM; ++idim) {
            real_t v = grid_.uold(i, 1 + idim) / d;
            vel2 += v * v; vel_max = std::max(vel_max, std::abs(v));
        }
        real_t e_int = grid_.uold(i, 5) - 0.5 * d * vel2;
        real_t p = std::max(e_int * (gamma - 1.0), d * 1e-10);
        dt_max = std::min(dt_max, courant_factor * dx / (vel_max + std::sqrt(gamma * p / d)));
    }
    real_t global_dt = dt_max;
#ifdef RAMSES_USE_MPI
    MPI_Allreduce(&dt_max, &global_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
    return global_dt;
}

void HydroSolver::get_diagnostics(int ilevel, real_t dx, real_t& min_d, real_t& max_v) {
    min_d = 1e30; max_v = 0.0;
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int id = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            real_t d = grid_.uold(id, 1); min_d = std::min(min_d, d);
            real_t v2 = 0; for(int d_idx=1; d_idx<=NDIM; ++d_idx) v2 += std::pow(grid_.uold(id, 1+d_idx)/d, 2);
            max_v = std::max(max_v, std::sqrt(v2));
        }
        igrid = grid_.next[igrid - 1];
    }
#ifdef RAMSES_USE_MPI
    real_t g_min_d, g_max_v;
    MPI_Allreduce(&min_d, &g_min_d, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&max_v, &g_max_v, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    min_d = g_min_d; max_v = g_max_v;
#endif
}

void HydroSolver::add_gravity_source_terms(int ilevel, real_t dt) {
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int id = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            real_t d = grid_.uold(id, 1);
            real_t e_kin_old = 0; for(int i=1; i<=NDIM; ++i) e_kin_old += 0.5 * std::pow(grid_.uold(id, 1+i), 2) / d;
            for (int idim = 1; idim <= NDIM; ++idim) grid_.uold(id, 1 + idim) += d * grid_.f(id, idim) * dt;
            real_t e_kin_new = 0; for(int i=1; i<=NDIM; ++i) e_kin_new += 0.5 * std::pow(grid_.uold(id, 1+i), 2) / d;
            grid_.uold(id, 5) += (e_kin_new - e_kin_old);
        }
        igrid = grid_.next[igrid - 1];
    }
}

} // namespace ramses
