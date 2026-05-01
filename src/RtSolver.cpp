#include "ramses/RtSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/MpiManager.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace ramses {

void RtSolver::initialize() {
    nGroups = config_.get_int("rt_params", "nGroups", 0);
    if (nGroups > 0) {
        load_hll_eigenvalues();
        chem_ = std::make_unique<RtChemistry>(nGroups);
    }
}

void RtSolver::load_hll_eigenvalues() {
    std::string path = config_.get("rt_params", "hll_evals_file", "legacy/rt/hll_evals.list");
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "[RtSolver] Error: Could not open HLL eigenvalues file: " << path << std::endl;
        return;
    }

    int n; file >> n;
    lambda1.assign(101, std::vector<real_t>(101, 0.0));
    lambda4.assign(101, std::vector<real_t>(101, 0.0));

    for (int i = 0; i < 101; ++i) {
        for (int j = 0; j < 101; ++j) {
            int ii, jj;
            real_t d1, d2;
            file >> ii >> jj >> lambda1[ii][jj] >> d1 >> d2 >> lambda4[ii][jj];
        }
    }
    std::cout << "[RtSolver] Loaded HLL eigenvalues from " << path << std::endl;
}

void RtSolver::godunov_fine(int ilevel, real_t dt, real_t dx) {
    if (nGroups <= 0) return;
    std::vector<int> active_octs;
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        active_octs.push_back(igrid);
        igrid = grid_.next[igrid - 1];
    }
    if (!active_octs.empty()) rt_godfine1(active_octs, ilevel, dt, dx);
}

void RtSolver::apply_source_terms(int ilevel, real_t dt) {
    if (nGroups <= 0 || !chem_) return;
    int myid = MpiManager::instance().rank() + 1;
    int nvar_hydro = 5; 
#ifdef MHD
    nvar_hydro = 8;
#endif
    int iIons = nvar_hydro + 1; // Assuming ions start after hydro

    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int id = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son[id] != 0) continue;

            real_t nH = grid_.uold(id, 1) * params::units_density / 1.67e-24;
            real_t T2 = (grid_.uold(id, 5) - 0.5 * grid_.uold(id, 2)*grid_.uold(id, 2)/grid_.uold(id, 1)) * (grid_.gamma - 1.0) / grid_.uold(id, 1) * params::units_pressure / params::units_density * 1.67e-24 / 1.38e-16;

            real_t xion[3] = {grid_.uold(id, iIons)/grid_.uold(id, 1), 0, 0};
            real_t Np[10], Fp[10][3];
            for (int ig = 0; ig < nGroups; ++ig) {
                Np[ig] = grid_.uold(id, iIons + 3 + ig * (1+NDIM));
                for(int d=0; d<NDIM; ++d) Fp[ig][d] = grid_.uold(id, iIons + 3 + ig * (1+NDIM) + 1 + d);
            }

            chem_->solve_chemistry(T2, xion, Np, Fp, nH, dt * params::units_time, 1.0);

            // Update state
            grid_.uold(id, iIons) = xion[0] * grid_.uold(id, 1);
            for (int ig = 0; ig < nGroups; ++ig) {
                grid_.uold(id, iIons + 3 + ig * (1+NDIM)) = Np[ig];
                for(int d=0; d<NDIM; ++d) grid_.uold(id, iIons + 3 + ig * (1+NDIM) + 1 + d) = Fp[ig][d];
            }
            // Update energy
            real_t mu = 1.0 / (1.0 + xion[0]);
            real_t p_new = nH * 1.38e-16 * T2 / mu / params::units_pressure;
            grid_.uold(id, 5) = p_new / (grid_.gamma - 1.0) + 0.5 * grid_.uold(id, 2)*grid_.uold(id, 2)/grid_.uold(id, 1);
        }
        igrid = grid_.next[igrid - 1];
    }
}

void RtSolver::rt_godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx) {
    int nvar_hydro = 5; 
#ifdef MHD
    nvar_hydro = 8;
#endif
    int iIons = nvar_hydro + 1;

    for (int igrid : ind_grid) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son[idc] != 0) continue;
            for (int ig = 0; ig < nGroups; ++ig) {
                int idx = iIons + 3 + ig * (1 + NDIM);
                grid_.unew(idc, idx) = std::max(smallNp, grid_.uold(idc, idx));
            }
        }
    }
}

void RtSolver::set_unew(int ilevel) {
    if (nGroups <= 0) return;
    int nvar_hydro = 5;
#ifdef MHD
    nvar_hydro = 8;
#endif
    int nrtvar = nGroups * (1 + NDIM) + 3; // +3 for ions
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            for (int iv = 1; iv <= nrtvar; ++iv) {
                grid_.unew(idc, nvar_hydro + iv) = grid_.uold(idc, nvar_hydro + iv);
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

void RtSolver::set_uold(int ilevel) {
    if (nGroups <= 0) return;
    int nvar_hydro = 5;
#ifdef MHD
    nvar_hydro = 8;
#endif
    int nrtvar = nGroups * (1 + NDIM) + 3;
    int myid = MpiManager::instance().rank() + 1;
    int igrid = grid_.headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            for (int iv = 1; iv <= nrtvar; ++iv) {
                grid_.uold(idc, nvar_hydro + iv) = grid_.unew(idc, nvar_hydro + iv);
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

} // namespace ramses
