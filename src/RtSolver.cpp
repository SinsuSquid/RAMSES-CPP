#include "ramses/RtSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace ramses {

void RtSolver::initialize() {
    nGroups = config_.get_int("rt_params", "nGroups", 0);
    if (nGroups > 0) {
        load_hll_eigenvalues();
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
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            active_octs.push_back(igrid);
            igrid = grid_.next[igrid - 1];
        }
    }
    
    if (!active_octs.empty()) {
        rt_godfine1(active_octs, ilevel, dt, dx);
    }
}

void RtSolver::rt_godfine1(const std::vector<int>& ind_grid, int ilevel, real_t dt, real_t dx) {
    real_t dt_dx = dt / dx;
    int nvar_hydro = 5; // Default hydro variables
#ifdef MHD
    nvar_hydro = 8;
#endif
    int nrtvar = nGroups * (1 + NDIM);

    for (int igrid : ind_grid) {
        // Simplified stencil and flux calculation for RT (1st order Godunov for now)
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            if (grid_.son[idc] != 0) continue;

            for (int ig = 0; ig < nGroups; ++ig) {
                int iNp = nvar_hydro + ig * (1 + NDIM) + 1;
                
                // Simple upwind or HLL flux update (TBD: full implementation)
                // For now, this is a placeholder for the full RT flux logic
                grid_.unew(idc, iNp) = std::max(smallNp, grid_.uold(idc, iNp));
            }
        }
    }
}

void RtSolver::set_unew(int ilevel) {
    int nrtvar = nGroups * (1 + NDIM);
    int nvar_hydro = 5;
#ifdef MHD
    nvar_hydro = 8;
#endif

    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
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
}

void RtSolver::set_uold(int ilevel) {
    int nrtvar = nGroups * (1 + NDIM);
    int nvar_hydro = 5;
#ifdef MHD
    nvar_hydro = 8;
#endif

    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
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
}

} // namespace ramses
