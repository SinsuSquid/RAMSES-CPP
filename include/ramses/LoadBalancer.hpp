#ifndef RAMSES_LOAD_BALANCER_HPP
#define RAMSES_LOAD_BALANCER_HPP

#include "AmrGrid.hpp"
#include "ParticleSystem.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Handles MPI domain decomposition and load balancing.
 */
class LoadBalancer {
public:
    LoadBalancer(AmrGrid& grid, ParticleSystem& ps) : grid_(grid), ps_(ps) {}

    void balance();

    void calculate_hilbert_keys();

    struct OctPacket {
        int ilevel;
        real_t xg[3];
        real_t father_x[3]; // Center of father cell
        int nbor[6];
        real_t uold[8 * 30]; 
        real_t unew[8 * 30];
        qdp_t hilbert_keys[8];
    };

private:
    void compute_new_cpu_map(std::vector<int>& cpu_map_new);
    void move_grids(const std::vector<int>& cpu_map_new);
    
    void remove_grid_from_list(int igrid, int ilevel, int icpu);

    AmrGrid& grid_;
    ParticleSystem& ps_;
};

} // namespace ramses

#endif // RAMSES_LOAD_BALANCER_HPP
