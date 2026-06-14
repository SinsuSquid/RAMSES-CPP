#ifndef RAMSES_FEEDBACK_SOLVER_HPP
#define RAMSES_FEEDBACK_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"

namespace ramses {

/**
 * @brief Handles supernova feedback from star particles and sinks.
 */
class FeedbackSolver {
public:
    FeedbackSolver(AmrGrid& grid, Config& config);
    ~FeedbackSolver();

    void init();
    void thermal_feedback(int ilevel, real_t dt);
    
    // Sink feedback (Phase 27)
    void sink_feedback(int ilevel, real_t dt);

private:
    void feedbk(int ilevel, real_t dt, const std::vector<int>& ind_part);

    AmrGrid& grid_;
    Config& config_;

    // Parameters
    bool metal_;
    real_t yield_;
    real_t eta_sn_;
    real_t f_w_;
    real_t mass_gmc_;
    real_t t_sne_;
    real_t kappa_IR_;
    real_t z_ave_;
    bool delayed_cooling_;
    bool use_proper_time_;

    // Sink specific
    real_t sn_p_ref_;
    real_t sn_e_ref_;
    real_t Tsat_;
};

} // namespace ramses

#endif // RAMSES_FEEDBACK_SOLVER_HPP
