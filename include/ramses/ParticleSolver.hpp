#ifndef RAMSES_PARTICLE_SOLVER_HPP
#define RAMSES_PARTICLE_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"

namespace ramses {

/**
 * @brief Handles particle dynamics (movement, force interpolation).
 */
class ParticleSolver {
public:
    ParticleSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}

    /**
     * @brief Assigns particle mass to the grid density field (Cloud-In-Cell).
     */
    void assign_mass(int ilevel);

    /**
     * @brief Pushes particles using current force field.
     */
    void move_fine(int ilevel, real_t dt);

    /**
     * @brief Re-links particles to their current leaf cells.
     */
    void relink();

    /**
     * @brief Assigns particle mass to the density grid (CIC).
     */
    void assign_mass_fine(int ilevel);

    /**
     * @brief Performs MPI particle exchange between ranks.
     */
    void exchange_particles();

    struct ParticlePacket {
        real_t xp[3];
        real_t vp[3];
        real_t mp;
        int idp;
        int levelp;
    };

    private:
    void move_particles(const std::vector<int>& ind_part, real_t dt);
    int find_cell_by_coords(real_t x, real_t y, real_t z, int* level = nullptr);

    AmrGrid& grid_;
    Config& config_;
    };

} // namespace ramses

#endif // RAMSES_PARTICLE_SOLVER_HPP
