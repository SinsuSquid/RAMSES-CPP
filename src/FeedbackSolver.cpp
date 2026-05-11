#include "ramses/FeedbackSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/MpiManager.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace ramses {

namespace p = ramses::params;

FeedbackSolver::FeedbackSolver(AmrGrid& grid, Config& config) 
    : grid_(grid), config_(config) {
    metal_ = config_.get_int("hydro_params", "nmetals", 0) > 0;
    yield_ = config_.get_double("feedback_params", "yield", 0.1);
    eta_sn_ = config_.get_double("feedback_params", "eta_sn", 0.1);
    f_w_ = config_.get_double("feedback_params", "f_w", 0.0);
    t_sne_ = config_.get_double("feedback_params", "t_sne", 1e7); // in years
    delayed_cooling_ = config_.get_bool("feedback_params", "delayed_cooling", false);
}

FeedbackSolver::~FeedbackSolver() {}

void FeedbackSolver::init() {}

void FeedbackSolver::thermal_feedback(int ilevel, real_t dt) {
    if (eta_sn_ <= 0) return;

    int myid = MpiManager::instance().rank() + 1;
    int ig = grid_.get_headl(myid, ilevel);
    
    real_t dx = p::boxlen / (real_t)(p::nx * (1 << (ilevel - 1)));
    real_t vol = std::pow(dx, NDIM);

    while (ig > 0) {
        int ip = grid_.headp[ig - 1];
        while (ip > 0) {
            int next_p = grid_.nextp[ip - 1];
            
            // Check if star particle
            if (grid_.family[ip - 1] == 2) { // FAM_STAR usually 2
                real_t age = (grid_.gamma > 0) ? (p::units_time * (p::tend - grid_.tp[ip - 1])) : 0; // Simplified age logic
                // In reality, we need to compare t_now - t_birth
                
                // For now, let's implement a simplified check:
                // If the star was formed in the PREVIOUS step, it explodes now.
                // In RAMSES, feedback is often instantaneous or delayed.
                
                // Let's assume instantaneous for the port verification if needed, 
                // or check the age window.
                
                // ... (Logic to determine if SN happens in this dt) ...
            }
            ip = next_p;
        }
        ig = grid_.next[ig - 1];
    }
}

void FeedbackSolver::sink_feedback(int ilevel, real_t dt) {
    // Placeholder for Phase 28/29
}

void FeedbackSolver::feedbk(int ilevel, real_t dt, const std::vector<int>& ind_part) {
    // Core logic for energy injection
}

} // namespace ramses
