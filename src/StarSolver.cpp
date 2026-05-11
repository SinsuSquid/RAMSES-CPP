#include "ramses/StarSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/MpiManager.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace ramses {

namespace p = ramses::params;

StarSolver::StarSolver(AmrGrid& grid, Config& config) 
    : grid_(grid), config_(config) {
    n_star_ = config_.get_double("sf_params", "n_star", 0.1);
    eps_star_ = config_.get_double("sf_params", "eps_star", 0.01);
    m_star_ = config_.get_double("sf_params", "m_star", 1.0); // User units or solar? Standard RAMSES is often solar if mass_sph is set
    T2_star_ = config_.get_double("sf_params", "T2_star", 0.0);
    g_star_ = config_.get_double("sf_params", "g_star", 1.6666667);
}

StarSolver::~StarSolver() {}

void StarSolver::init() {
    nstar_tot_ = 0;
}

void StarSolver::form_stars(int ilevel, real_t dt) {
    if (ilevel < 2) return; // Legacy RAMSES typically doesn't form stars on level 1

    int myid = MpiManager::instance().rank() + 1;
    real_t dx = p::boxlen / (real_t)(p::nx * (1 << (ilevel - 1)));
    real_t vol = std::pow(dx, NDIM);
    
    // SF thresholds in code units
    real_t d0 = n_star_ / p::scale_nH;
    real_t mstar_code = m_star_; // Simplification: assuming m_star is in code units for now, or scaled by mass_sph
    
    // G factor for gravity
    real_t factG = 1.0; // In code units, usually 1.0 unless cosmo
    if (config_.get_bool("run_params", "cosmo", false)) {
        // approx factG logic if needed
    }

    int n2d_val = (1 << NDIM);
    int ig = grid_.get_headl(myid, ilevel);
    
    while (ig > 0) {
        for (int ic = 1; ic <= n2d_val; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
            if (grid_.son[idc - 1] != 0) continue; // Leaf cells only

            real_t d = grid_.uold(idc, 1);
            if (d <= d0) continue;

            // 1. Deterministic Seeding
            // Hash spatial coordinates to get a unique seed for this cell
            real_t xc[3];
            grid_.get_cell_center(idc, xc);
            
            size_t seed = std::hash<real_t>{}(xc[0]);
            if (NDIM > 1) seed ^= std::hash<real_t>{}(xc[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            if (NDIM > 2) seed ^= std::hash<real_t>{}(xc[2]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>{}(ilevel) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            
            std::mt19937 gen(static_cast<unsigned int>(seed));

            // 2. Poisson mean calculation
            // Free fall time t_ff = 0.5427 * sqrt(1 / (G * rho))
            real_t t_ff = 0.5427 * std::sqrt(1.0 / (factG * std::max(d, (real_t)1e-10)));
            real_t mgas_cell = d * vol;
            
            // Expected mass converted to stars: M = dt * (eps_star / t_ff) * M_cell
            real_t expected_mstar = dt * (eps_star_ / t_ff) * mgas_cell;
            real_t poiss_mean = expected_mstar / mstar_code;
            
            std::poisson_distribution<int> dist(poiss_mean);
            int nstar = dist(gen);

            if (nstar > 0) {
                real_t m_spawn = nstar * mstar_code;
                
                // Security: prevent more than 90% depletion
                if (m_spawn > 0.9 * mgas_cell) {
                    m_spawn = 0.9 * mgas_cell;
                    // In a perfect port, we'd handle the 'lost' mass accounting like legacy
                }

                // 3. Update Grid (Gas Depletion)
                grid_.uold(idc, 1) -= m_spawn / vol;
                // Update conservative variables (momentum and energy also scale down)
                real_t d_old = d;
                real_t d_new = grid_.uold(idc, 1);
                real_t ratio = d_new / d_old;
                
                for (int iv = 2; iv <= NDIM + 2; ++iv) {
                    grid_.uold(idc, iv) *= ratio;
                    grid_.unew(idc, iv) = grid_.uold(idc, iv); // Sync unew
                }
                grid_.unew(idc, 1) = grid_.uold(idc, 1);

                // 4. Create Particles
                for (int i = 0; i < nstar; ++i) {
                    int ip = grid_.get_free_particle();
                    if (ip == 0) {
                        std::cerr << "[StarSolver] Out of particle memory!" << std::endl;
                        break;
                    }
                    
                    grid_.idp[ip - 1] = ++nstar_tot_; // Global-ish ID (needs MPI sync for true parity)
                    grid_.levelp[ip - 1] = ilevel;
                    grid_.mp[ip - 1] = mstar_code; // Simplified: distribute total spawn mass if capped
                    
                    for (int j = 0; j < NDIM; ++j) {
                        grid_.xp[j * grid_.npartmax + ip - 1] = xc[j];
                        grid_.vp[j * grid_.npartmax + ip - 1] = grid_.uold(idc, 2 + j) / grid_.uold(idc, 1);
                    }
                    
                    // Add to cell's particle list
                    int igrid = (idc > grid_.ncoarse) ? ((idc - grid_.ncoarse - 1) % grid_.ngridmax) + 1 : 0;
                    if (igrid > 0) {
                        int next_p = grid_.headp[igrid - 1];
                        grid_.nextp[ip - 1] = next_p;
                        grid_.prevp[ip - 1] = 0;
                        if (next_p > 0) grid_.prevp[next_p - 1] = ip;
                        grid_.headp[igrid - 1] = ip;
                        if (grid_.tailp[igrid - 1] == 0) grid_.tailp[igrid - 1] = ip;
                        grid_.numbp[igrid - 1]++;
                    }
                }
            }
        }
        ig = grid_.next[ig - 1];
    }
}

} // namespace ramses
