#ifndef RAMSES_SINK_SOLVER_HPP
#define RAMSES_SINK_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <vector>
#include <string>

namespace ramses {

struct Sink {
    int id;
    real_t x[3];
    real_t v[3];
    real_t mass;
    real_t ang_mom[3];
    real_t acc_rate;
    int level;
    real_t t_creation;
    bool merged = false;
};

/**
 * @brief Handles sink particle creation, accretion, and MPI synchronization.
 */
class SinkSolver {
public:
    SinkSolver(AmrGrid& grid, Config& config);
    ~SinkSolver();

    void init();
    void create_sinks(int ilevel);
    void grow_sinks(int ilevel, real_t dt);
    void merge_sinks();
    void synchronize_sinks();

private:
    AmrGrid& grid_;
    Config& config_;
    std::vector<Sink> sinks_;
    int nsinkmax_ = 100;
    real_t mass_sink_seed_ = 0.0;

    void collect_accretion_data(int ilevel, std::vector<real_t>& dM, std::vector<real_t>& dP);
};

} // namespace ramses

#endif // RAMSES_SINK_SOLVER_HPP
